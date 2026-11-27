"""Predictive quantiles, forecast performance metrics and bound diagnostics.

Python port of ``computeQuantiles.m``, ``computeWIS.m``,
``computeforecastperformance.m`` and ``boundHitFractions.m``.

Matrices of curves have one row per time point and one column per trajectory,
as in MATLAB.
"""

from __future__ import annotations

import warnings

import numpy as np

from .errors import GrowthPredictError
from .utils import quantile

ALPHAQUANTILES = np.array(
    [0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500,
     0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 0.975, 0.990]
)  # fmt: skip

# Central prediction intervals of the weighted interval score: 98%, 95%, 90%, ..., 10%.
WIS_ALPHAS = np.concatenate(([0.02, 0.05], np.arange(1, 10) / 10.0))


def _clip0(q):
    """MATLAB's ``(q>=0).*q``: negative limits become zero."""
    return np.where(q >= 0, q, 0.0)


def _row_quantiles(curves, probs):
    """Quantiles of each row (time point) across the columns (trajectories)."""
    curves = np.atleast_2d(np.asarray(curves, dtype=float))
    probs = np.atleast_1d(probs)
    if curves.shape[0] == 0:
        return np.zeros((0, probs.size))
    return np.asarray(quantile(curves, probs, axis=1)).T  # rows = time points


def computeQuantiles(data1, curvesforecasts2, forecastingperiod):
    """The 23 predictive quantiles of the calibration and forecast periods.

    Port of ``computeQuantiles.m``. Returns ``(quantilesc1, quantilesf1,
    alphaquantiles)``; negative quantiles are set to zero.
    """
    curves = np.atleast_2d(np.asarray(curvesforecasts2, dtype=float))
    tf2 = np.asarray(data1).shape[0]
    forecastingperiod = int(forecastingperiod)

    quantilesc1 = _clip0(_row_quantiles(curves[:tf2, :], ALPHAQUANTILES))

    quantilesf1 = np.zeros((max(forecastingperiod, 0), ALPHAQUANTILES.size))
    if forecastingperiod > 0:
        quantilesf1 = _clip0(_row_quantiles(curves[tf2 : tf2 + forecastingperiod, :], ALPHAQUANTILES))

    return quantilesc1, quantilesf1, ALPHAQUANTILES.copy()


def _wis_rows(y, curves):
    """Weighted interval score of each observation ``y[i]`` against row i of ``curves``."""
    alphas = WIS_ALPHAS
    K = alphas.size
    w0 = 0.5

    Lt = _clip0(_row_quantiles(curves, alphas / 2))
    Ut = _clip0(_row_quantiles(curves, 1 - alphas / 2))
    m = _row_quantiles(curves, [0.5])[:, 0]

    y = np.asarray(y, dtype=float).ravel()
    yc = y[:, None]
    IS = (Ut - Lt) + (2 / alphas) * (Lt - yc) * (yc < Lt) + (2 / alphas) * (yc - Ut) * (yc > Ut)
    sum1 = np.sum((alphas / 2) * IS, axis=1)

    return (1 / (K + 0.5)) * (w0 * np.abs(y - m) + sum1)


def computeWIS(data1, datalatest, curvesforecasts2, forecastingperiod):
    """Weighted interval score of the calibration and forecast periods.

    Port of ``computeWIS.m``.

    Returns
    -------
    WISC : float
        Mean WIS over the calibration period.
    WISFS : ndarray, h-by-2
        Rows ``[horizon, mean WIS of the first `horizon` forecast points]``;
        empty (0-by-2) without a forecast or when the data do not cover it.
    quantilesc1, quantilesf1 : ndarray
        The 23 predictive quantiles of both periods.
    """
    curves = np.atleast_2d(np.asarray(curvesforecasts2, dtype=float))
    datalatest = np.atleast_2d(np.asarray(datalatest, dtype=float))
    tf2 = np.asarray(data1).shape[0]
    forecastingperiod = int(forecastingperiod)

    datac = datalatest[:tf2, 1]

    WISC = float(np.mean(_wis_rows(datac, curves[:tf2, :])))
    quantilesc1 = _clip0(_row_quantiles(curves[:tf2, :], ALPHAQUANTILES))

    WISFS = np.zeros((0, 2))
    quantilesf1 = np.zeros((max(forecastingperiod, 0), ALPHAQUANTILES.size))

    if forecastingperiod > 0:

        # Forecasting period
        if datalatest.shape[0] < tf2 + forecastingperiod:
            warnings.warn("forecasting period is too long to evaluate with available data")
            return WISC, WISFS, quantilesc1, quantilesf1

        dataf = datalatest[tf2 : tf2 + forecastingperiod, 1]
        WISF = _wis_rows(dataf, curves[tf2 : tf2 + forecastingperiod, :])

        horizons = np.arange(1, forecastingperiod + 1)
        WISFS = np.column_stack([horizons, np.cumsum(WISF) / horizons])

        quantilesf1 = _clip0(_row_quantiles(curves[tf2 : tf2 + forecastingperiod, :], ALPHAQUANTILES))

    return WISC, WISFS, quantilesc1, quantilesf1


def computeforecastperformance(data1, datalatest, curvesforecasts1, curvesforecasts2, forecastingperiod):
    """Point errors, coverage of the 95% prediction interval and mean interval score.

    Port of ``computeforecastperformance.m``.

    Parameters
    ----------
    data1 : array, n-by-2
        Calibration window (only its length is used).
    datalatest : array, N-by-2
        Observations from the start of the window to the end of the series.
    curvesforecasts1 : array
        Model trajectories of the bootstrap refits (no observation noise);
        MAE, MSE and RMSE use their median.
    curvesforecasts2 : array
        Trajectories with observation noise; coverage and interval score use them.
    forecastingperiod : int

    Returns
    -------
    RMSECS, MSECS, MAECS, PICS, MISCS : ndarray, 1-by-1
        Calibration metrics.
    RMSEFS, MSEFS, MAEFS, PIFS, MISFS : ndarray, h-by-2
        Rows ``[horizon, metric over the first `horizon` forecast points]``;
        empty (0-by-2) without a forecast or when the data do not cover it.
    """
    c1 = np.atleast_2d(np.asarray(curvesforecasts1, dtype=float))
    c2 = np.atleast_2d(np.asarray(curvesforecasts2, dtype=float))
    datalatest = np.atleast_2d(np.asarray(datalatest, dtype=float))
    tf2 = np.asarray(data1).shape[0]
    forecastingperiod = int(forecastingperiod)

    # **** compute performance metrics ******

    # calibration period

    ycmean = np.median(c1[:tf2, :], axis=1)

    datac = datalatest[:tf2, 1]

    MSEc = np.mean((datac - ycmean) ** 2)
    RMSEc = np.sqrt(MSEc)
    MAEc = np.mean(np.abs(datac - ycmean))

    # 95% prediction coverage
    lims = _clip0(_row_quantiles(c2[:tf2, :], [0.025, 0.975]))
    Lt, Ut = lims[:, 0], lims[:, 1]

    PIc = 100 * (np.count_nonzero((datac >= Lt) & (datac <= Ut)) / datac.size)

    MISc = np.mean((Ut - Lt) + (2 / 0.05) * (Lt - datac) * (datac < Lt) + (2 / 0.05) * (datac - Ut) * (datac > Ut))

    RMSECS = np.array([[RMSEc]])
    MSECS = np.array([[MSEc]])
    MAECS = np.array([[MAEc]])
    PICS = np.array([[PIc]])
    MISCS = np.array([[MISc]])

    # forecasting period

    empty = np.zeros((0, 2))
    RMSEFS, MSEFS, MAEFS, PIFS, MISFS = empty, empty.copy(), empty.copy(), empty.copy(), empty.copy()

    if forecastingperiod > 0:

        if datalatest.shape[0] < tf2 + forecastingperiod:
            warnings.warn("forecasting period is too long to evaluate with available data")
            return RMSECS, MSECS, MAECS, PICS, MISCS, RMSEFS, MSEFS, MAEFS, PIFS, MISFS

        rows = slice(tf2, tf2 + forecastingperiod)

        yfmean = np.median(c1[rows, :], axis=1)
        dataf = datalatest[rows, 1]

        lims = _clip0(_row_quantiles(c2[rows, :], [0.025, 0.975]))
        Lt, Ut = lims[:, 0], lims[:, 1]

        horizons = np.arange(1, forecastingperiod + 1)

        sq = (dataf - yfmean) ** 2
        ab = np.abs(dataf - yfmean)
        inside = ((dataf >= Lt) & (dataf <= Ut)).astype(float)
        score = (Ut - Lt) + (2 / 0.05) * (Lt - dataf) * (dataf < Lt) + (2 / 0.05) * (dataf - Ut) * (dataf > Ut)

        MSEf = np.cumsum(sq) / horizons
        MAEf = np.cumsum(ab) / horizons
        PIf = 100 * np.cumsum(inside) / horizons
        MISf = np.cumsum(score) / horizons

        RMSEFS = np.column_stack([horizons, np.sqrt(MSEf)])
        MSEFS = np.column_stack([horizons, MSEf])
        MAEFS = np.column_stack([horizons, MAEf])
        PIFS = np.column_stack([horizons, PIf])
        MISFS = np.column_stack([horizons, MISf])

    return RMSECS, MSECS, MAECS, PICS, MISCS, RMSEFS, MSEFS, MAEFS, PIFS, MISFS


def boundHitFractions(Phatss, LB, UB, relTol=None):
    """Share of bootstrap parameter estimates sitting on the optimizer box.

    Port of ``boundHitFractions.m``.

    Parameters
    ----------
    Phatss : array, M-by-d
        Bootstrap estimates (rows = draws, columns = parameters).
    LB, UB : array
        Optimizer box: a vector of length d or more (one common box) or an
        M-by-d (or wider) matrix (the per-draw box).
    relTol : float
        Relative tolerance (default 1e-3). A draw P counts as hitting LB when
        ``P <= LB + relTol*|LB|`` and UB when ``P >= UB - relTol*|UB|``.

    Returns
    -------
    fracLB, fracUB : ndarray, length d
        Fractions in [0, 1]; NaN for pinned columns (LB == UB) and for columns
        whose draws are all NaN.
    hitLB, hitUB : ndarray of bool, M-by-d
        The individual hits.
    """
    if relTol is None or np.size(relTol) == 0:
        relTol = 1e-3

    Phatss = np.asarray(Phatss, dtype=float)
    if Phatss.ndim == 1:
        Phatss = Phatss.reshape(-1, 1)
    M, d = Phatss.shape

    def _box(B):
        B = np.asarray(B, dtype=float)
        # Broadcast a common box (given as a vector) to M-by-d. A column whose length
        # already equals M is a per-draw box with d == 1 and is kept as is.
        if B.ndim <= 1:
            B = np.tile(B.reshape(1, -1), (M, 1))
        elif min(B.shape) == 1 and B.shape[0] != M:
            B = np.tile(B.reshape(1, -1), (M, 1))
        return B

    LB = _box(LB)
    UB = _box(UB)

    if LB.shape[0] != M or UB.shape[0] != M or LB.shape[1] < d or UB.shape[1] < d:
        raise GrowthPredictError(
            "boundHitFractions:sizeMismatch",
            "LB/UB must be 1-by-d (or longer) or M-by-d (or wider) to match Phatss (%d-by-%d)." % (M, d),
        )

    # Trim a box longer than Phatss (the box also carries X0 / error-model entries).
    LB = LB[:, :d].copy()
    UB = UB[:, :d].copy()

    # Guard against flipped boxes.
    flipMask = LB > UB
    if np.any(flipMask):
        tmp = LB[flipMask].copy()
        LB[flipMask] = UB[flipMask]
        UB[flipMask] = tmp

    tolLB = relTol * np.abs(LB)
    tolUB = relTol * np.abs(UB)

    with np.errstate(invalid="ignore"):
        hitLB = Phatss <= LB + tolLB
        hitUB = Phatss >= UB - tolUB

    valid = ~np.isnan(Phatss)
    hitLB = hitLB & valid
    hitUB = hitUB & valid

    nValid = valid.sum(axis=0)
    fracLB = hitLB.sum(axis=0) / np.maximum(nValid, 1)
    fracUB = hitUB.sum(axis=0) / np.maximum(nValid, 1)

    # Pinned columns (LB==UB in every draw) and all-NaN columns carry no information.
    pinned = np.all(LB == UB, axis=0) if M > 0 else np.ones(d, dtype=bool)
    mask = pinned | (nValid == 0)
    fracLB = fracLB.astype(float)
    fracUB = fracUB.astype(float)
    fracLB[mask] = np.nan
    fracUB[mask] = np.nan

    return fracLB, fracUB, hitLB, hitUB
