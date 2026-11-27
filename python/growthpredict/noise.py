"""Observation noise for the parametric bootstrap and its dispersion estimates.

Python port of ``AddErrorStructure.m``, ``estimatePearsonDispersion.m`` and
``getMeanVarianceRatio.m``.
"""

from __future__ import annotations

import numpy as np

from .errors import GrowthPredictError
from .utils import Struct, get_rng


def AddErrorStructure(yi, M, dist1, factor1, d, rng=None):
    """Draw M noisy incidence curves around a model curve.

    Port of ``AddErrorStructure.m``.

    Parameters
    ----------
    yi : sequence
        CUMULATIVE model curve; the mean of observation t is ``yi[t]-yi[t-1]``.
    M : int
        Number of realizations (columns of the result).
    dist1 : int
        0 normal, sd = factor1
        1 Poisson
        2 negative binomial, VAR = factor1*mean
        3 negative binomial, VAR = mean + factor1*mean
        4 negative binomial, VAR = mean + factor1*mean^2
        5 negative binomial, VAR = mean + factor1*mean^d
    factor1, d : float
        Noise parameters (see ``dist1``).
    rng : numpy.random.Generator, optional
        Random stream; default is the toolbox stream (:func:`growthpredict.rng`).

    Returns
    -------
    curves : ndarray, len(yi)-by-M
        The first point of every realization is ``yi[0]``.

    Notes
    -----
    A zero model mean is replaced by 0.001 for the negative binomial cases, a
    negative binomial with no over-dispersion (VAR <= mean) is drawn as
    Poisson, and non-finite draws stop with an error instead of entering the
    bootstrap data.
    """
    rng = get_rng(rng)

    yi = np.asarray(yi, dtype=float).ravel()
    n = yi.size
    M = int(M)

    curves = np.zeros((n, M))
    if n == 0:
        return curves

    curves[0, :] = yi[0]

    if n > 1:
        signed = np.diff(yi)[:, None]  # normal case uses the signed increment
        means = np.abs(np.diff(yi))  # model mean of observations 2..n
        size = (n - 1, M)

        with np.errstate(all="ignore"):

            if dist1 == 0:  # normal

                draws = rng.normal(loc=np.broadcast_to(signed, size), scale=factor1)

            elif dist1 == 1:  # Poisson

                draws = _poisson(rng, np.broadcast_to(means[:, None], size))

            elif dist1 in (2, 3, 4, 5):

                mean1 = np.maximum(means, 0.001)
                if dist1 == 2:  # negative binomial, VAR = factor1*mean
                    var1 = mean1 * factor1
                elif dist1 == 3:  # negative binomial, VAR = mean + alpha*mean
                    var1 = mean1 + mean1 * factor1
                elif dist1 == 4:  # negative binomial, VAR = mean + alpha*mean^2
                    var1 = mean1 + factor1 * mean1**2
                else:  # negative binomial, VAR = mean + alpha*mean^d
                    var1 = mean1 + factor1 * mean1**d

                draws = _negbin(rng, mean1, var1, M)

            else:
                raise GrowthPredictError(
                    "GrowthPredict:AddErrorStructure:UnsupportedDist",
                    "dist1 = %g is not supported; use 0 to 5." % dist1,
                )

        curves[1:, :] = draws

    if np.any(~np.isfinite(curves)):
        raise GrowthPredictError(
            "GrowthPredict:AddErrorStructure:InvalidDraws",
            "Non-finite draws were generated (dist1 = %g, factor1 = %g). "
            "Check the fitted curve and the noise parameters." % (dist1, factor1),
        )

    return curves


def _poisson(rng, lam):
    """Poisson draws as floats; NaN where the mean is not a valid rate."""
    lam = np.asarray(lam, dtype=float)
    ok = np.isfinite(lam) & (lam >= 0) & (lam < 1e18)
    out = np.full(lam.shape, np.nan)
    if np.any(ok):
        out[ok] = rng.poisson(lam[ok])
    return out


def _negbin(rng, mean1, var1, M):
    """Negative binomial draws with the given mean and variance, one row per time point.

    Without over-dispersion (VAR <= mean, up to rounding) the limit is Poisson.
    A variance that gives no valid distribution (e.g. infinite) leaves NaN,
    which the caller reports.
    """
    n = mean1.size
    out = np.full((n, M), np.nan)
    with np.errstate(all="ignore"):
        over = var1 > mean1 * (1 + 1e-10)  # False also for NaN
        p1 = mean1 / var1
        r1 = mean1 * p1 / (1 - p1)
    finite = np.isfinite(mean1)

    rows_po = finite & ~over
    if np.any(rows_po):
        out[rows_po, :] = _poisson(rng, np.broadcast_to(mean1[rows_po][:, None], (int(rows_po.sum()), M)))

    rows_nb = finite & over & np.isfinite(r1) & (r1 > 0) & (p1 > 0) & (p1 < 1)
    if np.any(rows_nb):
        out[rows_nb, :] = rng.negative_binomial(
            r1[rows_nb][:, None], p1[rows_nb][:, None], size=(int(rows_nb.sum()), M)
        )
    return out


def estimatePearsonDispersion(ydata, fitcurve, numModelParams, nSeries=1, fixI0=0, minMean=0.5):
    """Variance-to-mean ratio of the data around the fitted curve.

    Port of ``estimatePearsonDispersion.m``. Pearson (quasi-Poisson) dispersion
    estimate for the noise model VAR = phi * mean, used to scale the negative
    binomial bootstrap noise (``dist1 = 2``) when the model is fitted by least
    squares (``method1 = 0``):

        phi = sum( (y_i - mu_i)^2 / mu_i ) / (n_used - numModelParams)

    Parameters
    ----------
    ydata, fitcurve : sequences
        Stacked observations and fitted curve (same length; one block of equal
        length per fitted series).
    numModelParams : int
        Estimated quantities of the CURVE (model parameters and free initial
        values), not counting the Gaussian variance of the least-squares fit;
        the runners pass ``numparams - 1``.
    nSeries : int
        Number of fitted series stacked in ``ydata`` (default 1).
    fixI0 : 0 or 1
        1 -> the first point of each series is fixed to the data by
        construction (zero residual) and is excluded.
    minMean : float
        Points with fitted mean below this are excluded (default 0.5).

    Returns
    -------
    phi : float
        The estimate; NaN when fewer than ``numModelParams + 1`` points are usable.
    info : Struct
        ``nUsed, dof, minMean, excludedSmallMean, excludedFixed``.
    """
    if nSeries is None:
        nSeries = 1
    if fixI0 is None:
        fixI0 = 0
    if minMean is None:
        minMean = 0.5

    y = np.asarray(ydata, dtype=float).ravel()
    mu = np.asarray(fitcurve, dtype=float).ravel()
    if y.size != mu.size:
        raise GrowthPredictError(
            "estimatePearsonDispersion:size",
            "ydata (%d) and fitcurve (%d) must have the same number of elements." % (y.size, mu.size),
        )
    if y.size % nSeries != 0:
        raise GrowthPredictError(
            "estimatePearsonDispersion:series",
            "%d observations cannot be split into %d series of equal length." % (y.size, nSeries),
        )
    nT = y.size // nSeries

    use = np.isfinite(y) & np.isfinite(mu)
    fixedMask = np.zeros(y.size, dtype=bool)
    if fixI0 and nT > 0:
        fixedMask[::nT] = True  # first point of every series
    smallMask = mu < minMean
    use = use & ~fixedMask & ~smallMask

    dof = int(np.count_nonzero(use)) - numModelParams
    if dof < 1:
        phi = float("nan")
    else:
        phi = float(np.sum((y[use] - mu[use]) ** 2 / mu[use]) / dof)

    info = Struct(
        nUsed=int(np.count_nonzero(use)),
        dof=dof,
        minMean=minMean,
        excludedSmallMean=int(np.count_nonzero(smallMask & ~fixedMask)),
        excludedFixed=int(np.count_nonzero(fixedMask)),
    )
    return phi, info


def getMeanVarianceRatio(data1, binsize1, flag1):
    """Variance-to-mean ratio of the raw data in bins (``getMeanVarianceRatio.m``).

    ``flag1 = 1`` uses consecutive bins of ``binsize1`` points, ``flag1 = 2`` a
    moving window. Returns ``(ratios, timebins1)``; each row of ``ratios`` is
    ``[variance/mean, mean, variance]``. The runners no longer use it (see
    :func:`estimatePearsonDispersion`); it is kept for earlier scripts.
    """
    data1 = np.asarray(data1, dtype=float).ravel()
    binsize1 = int(binsize1)
    ratios = []
    timebins1 = []

    if flag1 == 1:  # no moving window
        starts = range(0, data1.size - binsize1 + 1, binsize1)
    elif flag1 == 2:  # moving window
        starts = range(0, data1.size - binsize1 + 1, 1)
    else:
        starts = []

    for i in starts:
        chunk = data1[i : i + binsize1]
        mean1 = chunk.mean()
        var1 = chunk.var(ddof=1) if chunk.size > 1 else 0.0
        with np.errstate(all="ignore"):
            ratio1 = np.float64(var1) / np.float64(mean1)
        ratios.append([ratio1, mean1, var1])
        if flag1 == 1:
            timebins1.append(i + (binsize1 + 1) / 2)
        else:
            timebins1.append(i + 1)

    return np.array(ratios, dtype=float).reshape(-1, 3), np.array(timebins1, dtype=float)
