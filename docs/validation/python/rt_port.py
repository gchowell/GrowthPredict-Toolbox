"""Python ports of get_Rt.m: 'old' = repository HEAD b7829ae, 'new' = proposed.

Both return (Rs, ts) like the MATLAB function: ts is the 0.01 grid over the window and
Rs has one value fewer; Rs[k] is the estimate at ts[k+1].
"""
import numpy as np
from scipy.stats import gamma, expon
DT = 0.01


def _grid(timevect):
    n = int(np.floor((timevect[-1] - timevect[0]) / DT + 1e-9)) + 1     # MATLAB colon
    return timevect[0] + DT * np.arange(n)


def get_Rt_old(timevect, inc, gi_type, mean1, var1=0.0):
    ts = _grid(np.asarray(timevect, float))
    I = np.interp(ts, timevect, inc)
    n = len(I)
    if gi_type == 1:
        b = var1 / mean1; a = mean1 / b
        w = np.diff(gamma.cdf(ts, a, scale=b))              # absolute time labels
    elif gi_type == 2:
        w = np.diff(expon.cdf(ts, scale=mean1))
    else:
        idx = mean1 / DT + 1                                 # MATLAB 1-based, not rounded
        if idx != int(idx) or idx < 1:
            raise IndexError(f"Array indices must be positive integers (got {idx!r})")
        w = np.zeros(max(n, int(idx))); w[int(idx) - 1] = 1  # MATLAB grows the array
    return _renewal(I, w), ts


def get_Rt_new(timevect, inc, gi_type, mean1, var1=0.0):
    timevect = np.asarray(timevect, float)
    if gi_type not in (1, 2, 3):
        raise ValueError("type_GId1 must be 1, 2 or 3")
    if not (np.isfinite(mean1) and mean1 > 0):
        raise ValueError("mean1 must be positive")
    ts = _grid(timevect)
    I = np.interp(ts, timevect, inc)
    n = len(I)
    if n < 2:
        return np.array([]), ts
    tau = DT * np.arange(n)                                  # elapsed time since the window start
    if gi_type == 1:
        if not (np.isfinite(var1) and var1 > 0):
            raise ValueError("var1 must be positive")
        b = var1 / mean1; a = mean1 / b
        w = np.diff(gamma.cdf(tau, a, scale=b))
    elif gi_type == 2:
        w = np.diff(expon.cdf(tau, scale=mean1))
    else:
        lag = int(round(mean1 / DT))                         # number of grid steps
        if lag < 1:
            raise ValueError("fixed lag below the grid step")
        w = np.zeros(n - 1)
        if lag <= n - 1:
            w[lag - 1] = 1
    return _renewal(I, w), ts


def _renewal(I, w):
    n = len(I)
    c = np.convolve(I, w, 'full')[:n - 1]
    Rs = np.full(n - 1, np.nan)
    ok = c > 0
    Rs[ok] = I[1:][ok] / c[ok]
    return Rs


def gi_coverage(n_rs, gi_type, mean1, var1=0.0):
    """Third output of the new get_Rt.m: cumulative generation-interval mass covered at each Rs index."""
    tau = DT * np.arange(n_rs + 1)
    if gi_type == 1:
        b = var1 / mean1; a = mean1 / b
        w = np.diff(gamma.cdf(tau, a, scale=b))
    elif gi_type == 2:
        w = np.diff(expon.cdf(tau, scale=mean1))
    else:
        w = np.zeros(n_rs); lag = int(round(mean1 / DT))
        if 1 <= lag <= n_rs:
            w[lag - 1] = 1
    return np.cumsum(w)


def select_rt_rows(cov, min_cov=0.95, steps=100):
    """Port of selectRtRows.m. Returns (1-based row indices into Rs, boolean mask of reliable estimates)."""
    cov = np.asarray(cov, float)
    if not (np.isscalar(min_cov) and 0 < min_cov <= 1):
        raise ValueError("minCoverage must be in (0, 1]")
    if not (float(steps).is_integer() and steps >= 1):
        raise ValueError("stepsPerUnit must be a positive integer")
    n = len(cov)
    reliable = cov >= min_cov
    if n == 0:
        return np.array([], int), reliable
    idx = list(range(steps, n + 1, steps))
    if not idx or idx[-1] != n:
        idx.append(n)
    idx = np.array(idx, int)
    keep = idx[reliable[idx - 1]]
    if len(keep) == 0:
        keep = np.array([n], int)
    return keep, reliable
