"""Python ports of forecastPanelGrid.m and adjustLastWindowStart.m, and of the code they replace.

'old' = repository at v1.1.1 (commit 64648e7): the grid computed inline in plotForecast_GrowthModels.m,
and plotting functions that use the requested last window start as it is.
"""


def factor(n):
    """MATLAB factor(n): prime factors in ascending order; factor(1) = 1."""
    if n < 2:
        return [n]
    out, f = [], 2
    while f * f <= n:
        while n % f == 0:
            out.append(f)
            n //= f
        f += 1
    if n > 1:
        out.append(n)
    return out


def panel_grid_old(n):
    """plotForecast_GrowthModels.m, v1.1.1: rows and columns of figure 400."""
    f = factor(n)
    if len(f) == 1:
        return f[0], 1
    if len(f) == 3:
        return f[0] * f[1], f[2]
    return f[0], f[1]


def panel_grid_new(n):
    """forecastPanelGrid.m."""
    if not (isinstance(n, int) and n >= 1):
        raise ValueError("nWindows must be a positive integer")
    f = factor(n)
    if len(f) == 1:
        return f[0], 1
    if len(f) == 2:
        return f[0], f[1]
    if len(f) == 3:
        return f[0] * f[1], f[2]
    rows = 1
    for d in range(2, n + 1):
        if d * d > n:
            break
        if n % d == 0:
            rows = d
    return rows, n // rows


def adjust_last_window_start(tend1, windowsize1, n_observations):
    """adjustLastWindowStart.m, which is the rule of the two runners. Returns (tend1, adjusted)."""
    if (tend1 + windowsize1 - 1) > n_observations:
        return n_observations - windowsize1 + 1, True
    return tend1, False


def window_mat_name(cadfilename1, flag1, fixI0, method1, dist1, i, tend1, windowsize1, forecastingperiod):
    """Name of the per-window MAT file, as the runners write it and the plotting functions read it."""
    return ("./output/Forecast-growthModel-%s-flag1-%d-fixI0-%d-method-%d-dist-%d-tstart-%d-tend-%d-"
            "calibrationperiod-%d-forecastingperiod-%d.mat"
            % (cadfilename1, flag1, fixI0, method1, dist1, i, tend1, windowsize1, forecastingperiod))
