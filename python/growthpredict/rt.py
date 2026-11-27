"""Effective reproduction number from an incidence curve.

Python port of ``get_Rt.m``, ``selectRtRows.m`` and the generation-interval
option files ``options_Rt.m``, ``options_Rt_mpox.m`` and ``options_Rt_cholera.m``.

Indexing. MATLAB indices are one-based, Python's are zero-based. In MATLAB
``Rs(k)`` is the estimate at ``ts(k+1)``; here ``Rs[k]`` is the estimate at
``ts[k+1]``, so the pairing with ``ts[1:]`` is the same. The row indices
returned by :func:`selectRtRows` are zero-based: MATLAB's ``300:100:1900``
is ``299, 399, ..., 1899`` here, and ``ts[idxRows + 1]`` gives the same times.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.stats import expon, gamma

from .errors import GrowthPredictError

_DT1 = 0.01  # grid step of the renewal calculation


def _is_real_scalar(x):
    return np.ndim(x) == 0 and isinstance(x, (int, float, np.integer, np.floating)) and not isinstance(x, bool)


def _colon(a, step, b):
    """``a:step:b`` as MATLAB builds it (the last point never exceeds b)."""
    n = int(math.floor((b - a) / step + 1e-9)) + 1
    if n < 1:
        return np.zeros(0)
    return a + step * np.arange(n)


def get_Rt(timevect, incidence1, type_GId1, mean1, var1=None):
    """Reproduction number from an incidence curve with the renewal equation.

    Port of ``get_Rt.m``.

    Parameters
    ----------
    timevect : sequence
        Time labels of the incidence curve (increasing).
    incidence1 : sequence
        Incidence at those times.
    type_GId1 : int
        Generation-interval distribution: 1 = gamma, 2 = exponential,
        3 = fixed interval.
    mean1 : float
        Mean of the generation interval, in the units of ``timevect`` (for
        type 3: the fixed interval, rounded to the nearest 0.01).
    var1 : float
        Variance of the generation interval (gamma only).

    Returns
    -------
    Rs : ndarray
        ``len(ts) - 1`` values; ``Rs[k]`` is the estimate at ``ts[k+1]``, so
        plot it against ``ts[1:]``. NaN where the denominator is zero.
    ts : ndarray
        Grid ``timevect[0]:0.01:timevect[-1]``.
    giCoverage : ndarray
        Like ``Rs``: fraction of the generation-interval distribution that the
        incidence history available at ``ts[k+1]`` covers.

    Notes
    -----
    Incidence is interpolated linearly on the grid, and

        R(t) = I(t) / sum_s I(t-s) w(s),

    where w(s) is the probability that the generation interval lies in the
    grid step ending at s. The interval is measured from the start of the
    curve: s is elapsed time, so the result does not depend on the value of
    ``timevect[0]``. For a fixed interval g the estimate is ``I(t)/I(t-g)``.

    Only the incidence inside the curve is used. Until the history covers most
    of the generation-interval distribution the denominator is too small and
    the estimate too high; ``giCoverage`` shows when that is the case.
    """
    dt1 = _DT1

    if not _is_real_scalar(type_GId1) or type_GId1 not in (1, 2, 3):
        raise GrowthPredictError(
            "GrowthPredict:get_Rt:type",
            "type_GId1 must be 1 (gamma), 2 (exponential) or 3 (fixed interval).",
        )

    if not _is_real_scalar(mean1) or not math.isfinite(mean1) or mean1 <= 0:
        raise GrowthPredictError(
            "GrowthPredict:get_Rt:mean", "The mean of the generation interval must be a positive number."
        )

    if type_GId1 == 1 and (var1 is None or not _is_real_scalar(var1) or not math.isfinite(var1) or var1 <= 0):
        raise GrowthPredictError(
            "GrowthPredict:get_Rt:variance",
            "The variance of the gamma generation interval must be a positive number.",
        )

    timevect = np.asarray(timevect, dtype=float).ravel()
    incidence1 = np.asarray(incidence1, dtype=float).ravel()

    Rs = np.zeros(0)
    giCoverage = np.zeros(0)

    if timevect.size < 2:
        return Rs, timevect.copy(), giCoverage

    timevect2 = _colon(timevect[0], dt1, timevect[-1])

    incidence2 = np.interp(timevect2, timevect, incidence1)

    ts = timevect2

    n = incidence2.size

    if n < 2:
        return Rs, ts, giCoverage

    # Elapsed time since the start of the curve: the lag axis of the generation interval.
    tau = np.arange(n) * dt1

    if type_GId1 == 1:  # gamma

        b = var1 / mean1
        a = mean1 / b

        prob1 = np.diff(gamma.cdf(tau, a, scale=b))

    elif type_GId1 == 2:  # exponential

        prob1 = np.diff(expon.cdf(tau, scale=mean1))

    else:  # fixed interval

        lag1 = int(math.floor(mean1 / dt1 + 0.5))  # number of grid steps (MATLAB round)

        if lag1 < 1:
            raise GrowthPredictError(
                "GrowthPredict:get_Rt:lag",
                "A fixed generation interval must be at least %g (the grid step)." % dt1,
            )

        prob1 = np.zeros(n - 1)

        if lag1 <= n - 1:
            prob1[lag1 - 1] = 1

    # Renewal denominator by full convolution: for the estimate at grid point i,
    # c[i-1] = sum_{j=1}^{i} I[i-j] * prob1[j-1], where prob1[j-1] belongs to lag j*dt1.

    c = np.convolve(incidence2, prob1, mode="full")
    denominator = c[: n - 1]

    Rs = np.full(n - 1, np.nan)

    currentIncidence = incidence2[1:n]

    valid = denominator > 0

    Rs[valid] = currentIncidence[valid] / denominator[valid]

    giCoverage = np.cumsum(prob1)

    return Rs, ts, giCoverage


def selectRtRows(giCoverage, minCoverage=None, stepsPerUnit=None):
    """Reproduction-number estimates backed by enough incidence history.

    Port of ``selectRtRows.m``.

    Parameters
    ----------
    giCoverage : sequence
        Third output of :func:`get_Rt`.
    minCoverage : float
        Smallest fraction accepted (default 0.95).
    stepsPerUnit : int
        Grid steps of :func:`get_Rt` per unit of time (default 100).

    Returns
    -------
    idxRows : ndarray of int (zero-based)
        Indices into the estimates for a table with one row per unit of time
        after the start of the curve, ending at the last time, restricted to
        reliable estimates. Estimate ``idxRows[k]`` belongs to time
        ``ts[idxRows[k] + 1]``. When no estimate is reliable, ``idxRows`` is
        the last index alone, so that a table is never empty; its coverage
        shows that it is not reliable.
    reliable : ndarray of bool
        ``giCoverage >= minCoverage``. Use it to draw the curve:
        ``plot(ts[1:][reliable], Rs[reliable])``.

    Why
    ---
    :func:`get_Rt` uses only the incidence inside the curve. Near its start
    the renewal denominator misses the infections that happened earlier, and
    the estimate is too high: with the mpox default of :func:`options_Rt`,
    about 2.5 times the true value one week in and 20% above it two weeks in
    when incidence falls by 20% per week.
    """
    if minCoverage is None or np.size(minCoverage) == 0:
        minCoverage = 0.95

    if stepsPerUnit is None or np.size(stepsPerUnit) == 0:
        stepsPerUnit = 100

    if not _is_real_scalar(minCoverage) or not math.isfinite(minCoverage) or minCoverage <= 0 or minCoverage > 1:
        raise GrowthPredictError(
            "GrowthPredict:selectRtRows:minCoverage", "minCoverage must be a number above 0 and at most 1."
        )

    if (
        not _is_real_scalar(stepsPerUnit)
        or not math.isfinite(stepsPerUnit)
        or stepsPerUnit < 1
        or stepsPerUnit != round(stepsPerUnit)
    ):
        raise GrowthPredictError(
            "GrowthPredict:selectRtRows:stepsPerUnit", "stepsPerUnit must be a positive integer."
        )
    stepsPerUnit = int(stepsPerUnit)

    giCoverage = np.asarray(giCoverage, dtype=float).ravel()

    n = giCoverage.size

    reliable = giCoverage >= minCoverage

    idxRows = np.zeros(0, dtype=int)

    if n == 0:
        return idxRows, reliable

    # One index per unit of time after the start, and the last index.
    idxAll = np.arange(stepsPerUnit - 1, n, stepsPerUnit, dtype=int)

    if idxAll.size == 0 or idxAll[-1] != n - 1:
        idxAll = np.append(idxAll, n - 1)

    idxRows = idxAll[reliable[idxAll]]

    if idxRows.size == 0:
        idxRows = np.array([n - 1], dtype=int)

    return idxRows, reliable


# <=======================================================================================>
# <========================== Reproduction Number Parameters =============================>
# <=======================================================================================>


def options_Rt():
    """Default generation-interval options, ``(type_GId1, mean_GI1, var_GI1)`` (``options_Rt.m``).

    ``type_GId1``: 1 = gamma, 2 = exponential, 3 = delta (fixed interval).
    ``mean_GI1`` and ``var_GI1`` are in the time units of the data (and their
    square).

    The default data file is the weekly U.S. mpox series, so the default is the
    serial interval of mpox estimated from U.S. case pairs by symptom onset:
    mean 8.5 days (95% CrI 7.3-9.9), SD 5.0 days (95% CrI 4.0-6.4), written in
    weeks.

        Madewell ZJ, Charniga K, Masters NB, et al. Serial interval and
        incubation period estimates of monkeypox virus infection in 12
        jurisdictions, United States, May-August 2022. Emerg Infect Dis.
        2023;29(4). doi:10.3201/eid2904.221622

    To use other settings, write a function that returns the three values and
    pass it to the plotting function.
    """
    type_GId1 = 1  # gamma distribution
    mean_GI1 = 8.5 / 7  # mean generation interval: 8.5 days, in weeks
    var_GI1 = (5.0 / 7) ** 2  # variance: SD 5.0 days, in weeks^2
    return type_GId1, mean_GI1, var_GI1


def options_Rt_mpox():
    """Generation-interval options for mpox, weekly data (``options_Rt_mpox.m``).

    The values of :func:`options_Rt`: mean 8.5 days, SD 5.0 days (Madewell et
    al., 2023). For daily data use ``mean_GI1 = 8.5`` and ``var_GI1 = 5.0**2``.
    """
    type_GId1 = 1
    mean_GI1 = 8.5 / 7  # weeks
    var_GI1 = (5.0 / 7) ** 2  # weeks^2
    return type_GId1, mean_GI1, var_GI1


def options_Rt_cholera():
    """Generation-interval options for cholera, weekly data (``options_Rt_cholera.m``).

    Serial interval commonly assumed for cholera: gamma distribution with mean
    5 days and SD 8 days (e.g. Camacho A, et al. Lancet Glob Health.
    2018;6(6):e680-e690. doi:10.1016/S2214-109X(18)30230-4, which shifts it by
    one day). :func:`get_Rt` has no shift, so the distribution is applied from
    lag zero; about 41% of its mass lies within one day, which pulls Rt toward
    1. For daily data use ``mean_GI1 = 5`` and ``var_GI1 = 8**2``.
    """
    type_GId1 = 1
    mean_GI1 = 5 / 7  # weeks
    var_GI1 = (8 / 7) ** 2  # weeks^2
    return type_GId1, mean_GI1, var_GI1


RT_OPTIONS = {
    "options_Rt": options_Rt,
    "options_Rt_mpox": options_Rt_mpox,
    "options_Rt_cholera": options_Rt_cholera,
}
