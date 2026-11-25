"""Python port of getDoublingTimeCurve.m. prior = 0 reproduces repository HEAD b7829ae."""
import numpy as np


def doubling_times(curve, DT=1.0, prior=0.0):
    """curve: incidence at successive time points. Returns (doubling times, tds rows [time, cumulative])."""
    curve = np.asarray(curve, float)
    if len(curve) <= 1:
        return np.array([]), np.zeros((0, 2))
    tb = np.arange(len(curve)) * DT
    cum = prior + np.cumsum(curve)
    n = int(np.floor((tb[-1] - tb[0]) / 0.01 + 1e-9)) + 1
    XX = tb[0] + 0.01 * np.arange(n)
    x = np.interp(XX, tb, cum)
    C0 = x[0]
    rows = [(0.0, C0)]
    for _ in range(25):
        C0 = 2 * C0
        idx = np.where(x <= C0)[0]
        rows.append((XX[idx[-1]], x[idx[-1]]))
    tds = np.array(sorted(set(rows)))          # unique(tds,'rows')
    tds = tds[:-1]                             # drop the terminal row
    if len(tds) <= 1:
        return np.array([]), tds
    return np.diff(tds[:, 0])[:30], tds
