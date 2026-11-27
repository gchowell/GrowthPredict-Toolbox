"""Parameter estimation: initial guesses and the multi-start fit of one window.

Python port of ``initialParams.m`` and ``fit_model.m``.

MATLAB minimizes with ``fmincon`` (algorithm ``sqp``) started from every point
of a ``MultiStart`` set. This port runs SciPy's SQP implementation (``SLSQP``)
from the same kind of start set. Everything around the optimizer is the same:
the parameter box, the seed, the unit-box / log10 optimizer coordinates, the
objective scaling, the bounded transform of failed solves, the
rescale-and-confirm pass, the single final evaluation and the retry with a
wider interval for ``r``.
"""

from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np
from scipy.optimize import Bounds, minimize

from .errors import GrowthPredictError, warn
from .models import growthIncidenceCurve
from .objective import growthObjective
from .utils import Struct, get_rng

__all__ = ["initialParams", "fit_model", "FitResult"]

_REALMIN = np.finfo(float).tiny


def initialParams(data1, flag1):
    """Initial guess ``[r, p, a, K, I0, alpha, d]`` for the model ``flag1`` (``initialParams.m``)."""
    data1 = np.asarray(data1, dtype=float).ravel()
    s = float(np.sum(data1))
    y1 = float(data1[0])

    if flag1 == -1:
        params0 = [1, 1, 1, s, y1, 1, 1]
    elif flag1 == 0:
        params0 = [1, 0.9, 1, s, y1, 1, 1]
    elif flag1 == 1:
        params0 = [0.1, 0.9, 1, s, y1, 1, 1]
    elif flag1 == 2:
        params0 = [0.1, 0.9, 1, s, y1, 1, 1]
    elif flag1 == 3:
        params0 = [0.1, 1, 1, s, y1, 1, 1]
    elif flag1 == 4:
        params0 = [0.1, 1, 0.9, s, y1, 1, 1]
    elif flag1 == 5:
        K_guess = s
        with np.errstate(all="ignore"):
            r_guess = 1 - y1 / K_guess  # r = 1 - C0/K
            a_guess = float(np.float64(r_guess) / np.log(np.float64(K_guess) / np.float64(y1)))  # r/log(K/C0)
        params0 = [r_guess + 0.001, 1, a_guess, s, y1, 1, 1]
    else:
        params0 = [0.1, 0.9, 1, s, y1, 1, 1]

    return np.array(params0, dtype=float)


class FitResult(NamedTuple):
    """The nine outputs of ``fit_model.m``, in the same order.

    Unpack it like the MATLAB call,
    ``P, residual, fitcurve, forecastcurve, timevect2, initialguess, fval, flagg, fitDiagnostics = fit_model(...)``,
    or use the names (``result.P``, ``result.fval`` ...).
    """

    P: np.ndarray
    residual: np.ndarray
    fitcurve: np.ndarray
    forecastcurve: np.ndarray
    timevect2: np.ndarray
    initialguess: np.ndarray
    fval: float
    flagg: float
    fitDiagnostics: Struct


# <============================================================================>
# <============================ Local functions ===============================>
# <============================================================================>

_T_BOUND = 1e6
_FAILGAP = 1e3


def _localBound(v):
    """Monotone, bounded transform of the scaled objective handed to the optimizer.

    v <= T       : g = v
    v >  T       : g = T + log1p(v - T)   (continuous, <= T + 710)
    v not finite : g = T + FAILGAP        (a failed solve)

    A failed trial is therefore worse than every finite objective value.
    """
    if not math.isfinite(v):
        return _T_BOUND + _FAILGAP
    if v > _T_BOUND:
        return _T_BOUND + math.log1p(v - _T_BOUND)
    return v


def _localUnbound(g):
    """Inverse of :func:`_localBound`; the failure value maps back to ``inf``."""
    if not math.isfinite(g) or g >= _T_BOUND + _FAILGAP:
        return math.inf
    if g > _T_BOUND:
        try:
            return _T_BOUND + math.expm1(g - _T_BOUND)
        except OverflowError:
            return math.inf
    return g


def _localToT(Z, tf):
    T = np.array(Z, dtype=float, copy=True)
    if np.any(tf):
        T[..., tf] = np.log10(np.maximum(T[..., tf], _REALMIN))
    return T


def _localFromT(T, tf):
    Z = np.array(T, dtype=float, copy=True)
    if np.any(tf):
        with np.errstate(over="ignore"):
            Z[..., tf] = 10.0 ** Z[..., tf]
    return Z


def _lhsdesign(n, p, rng, iterations=50):
    """Latin hypercube sample, best of ``iterations`` by the maximin criterion.

    Stands in for ``lhsdesign(n,p,'criterion','maximin','iterations',50)``:
    each column is a random permutation of the n strata with a uniform draw
    inside each stratum; among the candidate designs the one with the largest
    minimum distance between points is kept.
    """
    best = None
    best_score = -math.inf
    for _ in range(max(1, int(iterations))):
        X = np.empty((n, p))
        for j in range(p):
            X[:, j] = (rng.permutation(n) + rng.random(n)) / n
        if n < 2:
            return X
        diff = X[:, None, :] - X[None, :, :]
        d2 = np.einsum("ijk,ijk->ij", diff, diff)
        score = d2[np.triu_indices(n, 1)].min()
        if score > best_score:
            best_score = score
            best = X
    return best


def _localMakeStarts(nBox, nJit, LBt, UBt, t0, free, rng):
    """Start points in the optimizer's (possibly log10) coordinates, one per row.

    nBox points cover the box of the free coordinates (Latin hypercube),
    nJit points are jittered copies of the seed (sd = 5% of each range),
    the seed itself is always the last row.
    """
    nz = t0.size
    nF = free.size
    span = UBt - LBt

    S1 = np.tile(t0, (nBox, 1))
    if nBox > 0 and nF > 0:
        X = _lhsdesign(nBox, nF, rng)
        S1[:, free] = LBt[free] + X * span[free]

    S2 = np.tile(t0, (nJit, 1))
    if nJit > 0 and nF > 0:
        S2[:, free] = S2[:, free] + 0.05 * rng.standard_normal((nJit, nF)) * span[free]
        S2 = np.minimum(np.maximum(S2, LBt), UBt)

    starts = np.vstack([S1, S2, t0.reshape(1, nz)])

    # unique(starts,'rows','stable')
    keep = []
    seen = set()
    for i in range(starts.shape[0]):
        key = starts[i].tobytes()
        if key not in seen:
            seen.add(key)
            keep.append(i)
    return starts[keep].reshape(-1, nz)


class _EvaluationLimit(Exception):
    pass


# fmincon option names accepted in fitopts['optimOptions'] and what they set here.
_DEFAULT_OPTIM = {
    "StepTolerance": 1e-10,
    "FunctionTolerance": 1e-10,
    "OptimalityTolerance": 1e-8,
    "FiniteDifferenceStepSize": 1e-6,
    "MaxFunctionEvaluations": 20000,
    "MaxIterations": 20000,
}


def _local_solve(fun, u0, nFree, optim):
    """One bounded local minimization in the unit box (stand-in for ``fmincon``/sqp).

    Returns ``(u, g, exitflag)`` with MATLAB-like exit flags: 1 converged,
    0 ran out of iterations or evaluations, -2 no usable point.

    SLSQP is run in rounds. Each round minimizes the objective divided by its
    value at the start of the round, because SLSQP, unlike fmincon, is not
    invariant to the scale of the objective: from a start point where the
    objective is 1e5 times its value at the seed it stops at once. A round
    starts from the end point of the previous one and the rounds stop when the
    objective no longer decreases, which is what the step-size test of fmincon
    achieves.
    """
    counter = {"n": 0, "best_g": math.inf, "best_u": None}
    maxfev = int(optim["MaxFunctionEvaluations"])
    scale = [1.0]

    def counted(u):
        if counter["n"] >= maxfev:
            raise _EvaluationLimit()
        counter["n"] += 1
        g = fun(u)
        if g < counter["best_g"]:
            counter["best_g"] = g
            counter["best_u"] = np.array(u, dtype=float, copy=True)
        return g / scale[0]

    bounds = Bounds(np.zeros(nFree), np.ones(nFree), keep_feasible=True)
    max_rounds = 5
    iters_left = int(optim["MaxIterations"])
    options = {
        "ftol": float(optim["FunctionTolerance"]),
        "eps": float(optim["FiniteDifferenceStepSize"]),
        "disp": False,
    }

    u = np.minimum(np.maximum(np.asarray(u0, dtype=float).ravel(), 0.0), 1.0)
    exitflag = -2
    try:
        g_prev = counted(u) * scale[0]
        for _ in range(max_rounds):
            if not math.isfinite(g_prev) or iters_left <= 0:
                break
            scale[0] = max(1.0, abs(g_prev))
            options["maxiter"] = int(min(iters_left, 200))
            res = minimize(counted, u, method="SLSQP", bounds=bounds, options=options)
            iters_left -= int(getattr(res, "nit", 0))
            u_new = np.minimum(np.maximum(np.asarray(res.x, dtype=float), 0.0), 1.0)
            g_new = float(res.fun) * scale[0]
            if not math.isfinite(g_new):
                break
            if res.status == 9:
                exitflag = 0  # iteration limit
            elif res.success or res.status == 8:
                # status 8 ('positive directional derivative for linesearch') is how
                # SLSQP reports that no further descent is possible at the precision
                # of the finite-difference gradient.
                exitflag = 1
            else:
                exitflag = -2
            improved = g_prev - g_new > 1e-8 * max(1.0, abs(g_new))
            if g_new <= g_prev:
                u = u_new
                g_prev = g_new
            if exitflag < 0 or not improved:
                if exitflag == 0 and not improved:
                    exitflag = 1  # stalled at the precision of the objective
                break
    except _EvaluationLimit:
        exitflag = 0

    # Return the best point visited (SLSQP can end on a point slightly worse
    # than one it evaluated on the way).
    if counter["best_u"] is None or not math.isfinite(counter["best_g"]):
        return u, math.inf, -2
    u_best = np.minimum(np.maximum(counter["best_u"], 0.0), 1.0)
    return u_best, float(counter["best_g"]), exitflag


def _multistart(fqs, starts_q, nFree, optim):
    """Run the local solver from every start point (stand-in for ``MultiStart``).

    Returns ``(x, fval, exitflag, runs)``. ``runs`` holds one record per start;
    the overall exit flag follows MultiStart: 1 all runs converged, 2 some
    did, 0 none converged but a run stopped at its iteration limit, -8 no
    solution.
    """
    runs = []
    for k in range(starts_q.shape[0]):
        try:
            u, g, flag = _local_solve(fqs, starts_q[k], nFree, optim)
        except Exception:  # a failure inside the local solver: skip this start
            continue
        runs.append(Struct(X=u, Fval=g, Exitflag=flag, X0=starts_q[k]))

    if not runs:
        return None, math.inf, -10, runs

    converged = [r for r in runs if r.Exitflag > 0]
    if converged:
        best = min(converged, key=lambda r: r.Fval)
        flag = 1 if len(converged) == len(runs) else 2
        return best.X, best.Fval, flag, runs

    stalled = [r for r in runs if r.Exitflag == 0]
    if stalled:
        best = min(stalled, key=lambda r: r.Fval)
        return best.X, best.Fval, 0, runs

    return None, math.inf, -8, runs


_parallel_warned = False


# <============================================================================>
# <================================ fit_model =================================>
# <============================================================================>


def fit_model(data1, params0, fixI0, numstartpoints, DT, flagX, forecastingperiod, fitopts=None, method1=0, rng=None):
    """Fit one growth model to one calibration window.

    Port of ``fit_model.m``.

    Parameters
    ----------
    data1 : array, n-by-2
        Time index and incidence of the calibration window.
    params0 : sequence of 7
        Initial guess ``[r, p, a, K, I0, alpha, d]`` (see :func:`initialParams`).
    fixI0 : 0 or 1
        1 fixes the initial state to the first observation, 0 estimates it.
    numstartpoints : int
        Number of start points that cover the parameter box; the seed and up to
        three jittered copies of it are added.
    DT : float
        Time step.
    flagX : int
        Growth model (-1 EXP, 0 GGM, 1 GLM, 2 GRM, 3 LM, 4 RICH, 5 GOM).
    forecastingperiod : int
        Number of steps to forecast beyond the window (0 for a fit only).
    fitopts : dict, optional
        Every key may be omitted:

        ``logscale``      True (default): the optimizer and the random start
                          points work in log10 of r, K and alpha.
        ``useParallel``   accepted for compatibility; the local solves run
                          serially (one warning).
        ``optimOptions``  dict of optimizer overrides using the fmincon names
                          ``FunctionTolerance``, ``FiniteDifferenceStepSize``,
                          ``MaxIterations``, ``MaxFunctionEvaluations``.
        ``LB``, ``UB``    7-vectors that replace the bounds derived from
                          ``data1`` and make ``params0`` the seed as given. The
                          runners pass the bounds of the fit to the observed
                          data to every bootstrap refit.
        ``rBounds``       ``[lower, upper]`` with 0 < lower <= upper: replaces
                          the interval of r derived from ``data1`` (any model)
                          and switches off the retry with a wider interval.
    method1 : int
        Estimation method: 0 least squares, 1 Poisson MLE, 3/4/5 negative
        binomial MLE. (MATLAB reads it from the global ``method1`` that the
        options file sets; here it is an argument.)
    rng : numpy.random.Generator, optional
        Random stream for the start points; default is the toolbox stream
        (see :func:`growthpredict.rng`).

    Returns
    -------
    FitResult
        ``(P, residual, fitcurve, forecastcurve, timevect2, initialguess, fval,
        flagg, fitDiagnostics)``, as in MATLAB. ``P = [r, p, a, K, I0, alpha, d]``.
        ``fitDiagnostics`` holds the bounds, convergence flags, the share of
        converged starts that agree with the best one (``startAgreement``) and
        which estimated coordinates ended on a bound (``atLB``, ``atUB``).
        ``freeIndex`` is zero-based.

    Notes
    -----
    Interval of r. For the models with a free exponent p (GGM, GLM, GRM) r has
    units of cases^(1-p)/time and is searched in ``[1e-4, 10*max(y)]/DT``. For
    EXP, LM and RICH (p fixed at 1) it is within a factor 20 of the slope of
    log incidence over the first 7 points; when the fit ends on that bound it
    is repeated once with r in ``[1e-3, 10]`` and the better fit is kept
    (``fitDiagnostics.rBoundsWidened``). Gompertz: a factor 10 either side of
    the guesses of :func:`initialParams`, per unit of time.
    """
    global _parallel_warned

    rng = get_rng(rng)
    flag1 = flagX

    if method1 is None:
        raise GrowthPredictError(
            "GrowthPredict:fit_model:NoMethod",
            "method1 is empty. It is set inside the options (options_fit / options_forecast).",
        )

    data1 = np.asarray(data1, dtype=float)
    if data1.ndim != 2 or data1.shape[1] < 2:
        raise GrowthPredictError(
            "GrowthPredict:fit_model:InvalidData", "data1 must have two columns: time index and incidence."
        )
    if np.any(~np.isfinite(data1)):
        raise GrowthPredictError(
            "GrowthPredict:fit_model:InvalidData", "fit_model received non-finite data (NaN/Inf)."
        )

    params0 = np.asarray(params0, dtype=float).ravel()
    fitopts = dict(fitopts) if isinstance(fitopts, dict) else None

    timevect = data1[:, 0] * DT

    I0 = data1[0, 1]  # initial condition

    z = np.array(params0[:7], dtype=float, copy=True)
    z[4] = I0

    if method1 == 0 or method1 == 1:
        LBe = [0, 0]
        UBe = [0, 0]
    elif method1 == 3:
        LBe = [1e-8, 1]
        UBe = [1e3, 1]
    elif method1 == 4:
        LBe = [1e-8, 1]
        UBe = [1e5, 1]
    elif method1 == 5:
        LBe = [1e-8, 0.6]  # d>=1
        UBe = [1e5, 1e3]
    else:
        raise GrowthPredictError(
            "GrowthPredict:fit_model:UnsupportedMethod",
            "method1 = %g is not supported; use 0, 1, 3, 4 or 5." % method1,
        )

    w = min(7, data1.shape[0])  # early window

    if w < 2:  # degenerate: too few points
        r0 = params0[0]
        rlb, rub = 1e-3, 10.0
    else:
        y = np.maximum(data1[:w, 1], 0.5)  # guard zeros
        b = np.polyfit(data1[:w, 0] * DT, np.log(y), 1)

        if b[0] <= 0:
            r0 = params0[0]  # no growth signal: keep supplied guess
            rlb, rub = 1e-3, 10.0  # wide default
        else:
            r0 = b[0]
            rlb, rub = r0 / 20, r0 * 20

    # Models with a free growth exponent p (GGM, GLM, GRM): r has units of
    # cases^(1-p)/time, so an interval taken from the slope of log incidence
    # (units 1/time) is only meaningful for p = 1. The interval is taken from the
    # scale of the data instead; r0 above is kept as the seed.
    if flag1 in (0, 1, 2):
        ymax = max(np.max(data1[:, 1]), 1)
        rlb = 1e-4 / DT
        rub = 10 * ymax / DT

    Cobs = np.sum(np.abs(data1[:, 1]))
    Klb = max(20, Cobs)  # can't undershoot observed cumulative
    Kub = min(1e3 * Cobs, 1e9)  # capped
    Kub = max(Kub, Klb)  # keep Klb <= Kub for very small series

    if fixI0 == 1:
        I0lb = I0ub = I0
    else:
        y1 = data1[0, 1]
        I0lb = max(1, 0.2 * y1)
        I0ub = 5 * max(1, y1)
        I0ub = min(I0ub, Klb)  # keep I0 <= K for saturating models
        I0lb = min(I0lb, I0ub)  # preserve I0lb <= I0ub

    if flag1 == -1:  # EXP
        LB = [rlb, 1, 1, 0, I0lb] + LBe
        UB = [rub, 1, 1, 0, I0ub] + UBe
    elif flag1 == 0:  # GGM
        LB = [rlb, 0.01, 1, 0, I0lb] + LBe
        UB = [rub, 1, 1, 0, I0ub] + UBe
    elif flag1 == 1:  # GLM
        LB = [rlb, 0.01, 1, Klb, I0lb] + LBe
        UB = [rub, 1, 1, Kub, I0ub] + UBe
    elif flag1 == 2:  # GRM
        LB = [rlb, 0.01, 0, Klb, I0lb] + LBe
        UB = [rub, 1, 10, Kub, I0ub] + UBe
    elif flag1 == 3:  # Logistic: the exponent a is not part of the model, keep it fixed at 1
        LB = [rlb, 1, 1, Klb, I0lb] + LBe
        UB = [rub, 1, 1, Kub, I0ub] + UBe
    elif flag1 == 4:  # Richards
        LB = [rlb, 1, 0, Klb, I0lb] + LBe
        UB = [rub, 1, 10, Kub, I0ub] + UBe
    elif flag1 == 5:  # Gompertz: r and a are rates, so the bounds are per unit of time
        LB = [max(1e-4, params0[0] / 10) / DT, 1, max(1e-6, params0[2] / 10) / DT, 1, I0lb] + LBe
        UB = [params0[0] * 10 / DT, 1, params0[2] * 10 / DT, 1, I0ub] + UBe
    else:
        raise GrowthPredictError(
            "GrowthPredict:fit_model:UnknownModel",
            "flag1 = %g is not a supported growth model (-1 to 5)." % flag1,
        )

    LB = np.array(LB, dtype=float)
    UB = np.array(UB, dtype=float)

    # Caller-supplied interval for r (fitopts['rBounds']): replaces the interval derived
    # above, for every model, and switches off the retry with a wider interval.
    rBoundsUser = fitopts is not None and fitopts.get("rBounds") is not None and np.size(fitopts["rBounds"]) > 0
    if rBoundsUser:
        try:
            rb = np.asarray(fitopts["rBounds"], dtype=float).ravel()
        except (TypeError, ValueError):
            rb = np.array([])
        if rb.size != 2 or not np.all(np.isfinite(rb)) or rb[0] <= 0 or rb[0] > rb[1]:
            raise GrowthPredictError(
                "GrowthPredict:fit_model:rBounds",
                "fitopts.rBounds must be [lower upper] with 0 < lower <= upper.",
            )
        LB[0] = rb[0]
        UB[0] = rb[1]

    if flag1 != 5:
        z[0] = min(max(r0, LB[0]), UB[0])  # informed r guess
    else:
        z[[0, 2]] = params0[[0, 2]] / DT  # Gompertz guesses of initialParams are per observation step

    # Caller-supplied box (bootstrap refits): use it as is and keep params0 as the seed.
    boundsSupplied = (
        fitopts is not None
        and fitopts.get("LB") is not None
        and fitopts.get("UB") is not None
        and np.size(fitopts["LB"]) > 0
        and np.size(fitopts["UB"]) > 0
    )
    if boundsSupplied:
        if np.size(fitopts["LB"]) != z.size or np.size(fitopts["UB"]) != z.size:
            raise GrowthPredictError(
                "GrowthPredict:fit_model:SuppliedBounds",
                "fitopts.LB and fitopts.UB must have %d entries [r p a K I0 alpha d]." % z.size,
            )
        LB = np.asarray(fitopts["LB"], dtype=float).ravel().copy()
        UB = np.asarray(fitopts["UB"], dtype=float).ravel().copy()
        z = np.array(params0[:7], dtype=float, copy=True)

    nz = z.size

    if LB.size != nz or UB.size != nz:
        raise GrowthPredictError(
            "GrowthPredict:fit_model:Length",
            "Length mismatch: numel(z)=%d, numel(LB)=%d, numel(UB)=%d." % (nz, LB.size, UB.size),
        )
    if np.any(LB > UB) or not np.all(np.isfinite(LB)) or not np.all(np.isfinite(UB)):
        raise GrowthPredictError(
            "GrowthPredict:fit_model:Bounds",
            "The parameter bounds must be finite with LB <= UB (check the data and the initial guesses).",
        )

    z0 = np.minimum(np.maximum(z, LB), UB)  # clamp everything into bounds
    mid = (LB + UB) / 2
    nanmask = np.isnan(z0)
    z0[nanmask] = mid[nanmask]

    ydata = data1[:, 1].copy()

    # <============================================================================>
    # <================= Optimizer coordinates: free parameters, unit box =========>
    # <============================================================================>

    def _flag(name, default):
        if fitopts is not None and fitopts.get(name) is not None:
            return bool(fitopts[name])
        return default

    useLog = _flag("logscale", True)
    useParallel = _flag("useParallel", False)

    if useParallel and not _parallel_warned:
        _parallel_warned = True
        warn(
            "GrowthPredict:fit_model:Parallel",
            "Parallel local solves are not available in the Python port; running the multi-start serially.",
        )

    free = np.flatnonzero(UB > LB)
    nFree = free.size

    tf = np.zeros(nz, dtype=bool)  # coordinates searched in log10
    if useLog:
        logCandidates = np.array([0, 3, 5])  # r, K, alpha
        tf[logCandidates] = (UB[logCandidates] > LB[logCandidates]) & (LB[logCandidates] > 0)

    LBt = _localToT(LB, tf)
    UBt = _localToT(UB, tf)
    spanT = UBt - LBt
    z0T = _localToT(z0, tf)

    def toQ(Z):
        """Rows of Z (original units) -> rows of u (unit box, free coordinates only)."""
        T = _localToT(np.atleast_2d(Z), tf)
        U = (T[:, free] - LBt[free]) / spanT[free]
        return np.minimum(np.maximum(U, 0.0), 1.0)

    def fromQ(U):
        """A point u (unit box, free coordinates) -> z (original units, full vector).

        Fixed coordinates are taken from the seed z0.
        """
        U = np.minimum(np.maximum(np.asarray(U, dtype=float).ravel(), 0.0), 1.0)
        T = z0T.copy()
        T[free] = LBt[free] + U * spanT[free]
        return _localFromT(T, tf)

    # Objective in original units.
    def f(Z):
        return growthObjective(Z, timevect, ydata, flag1, method1)

    def fq(U):
        return f(fromQ(U))

    f0 = f(z0)
    if math.isfinite(f0) and abs(f0) >= 1:
        fscale = abs(f0)
    else:
        fscale = 1.0

    def fqs(U, _scale=None):
        return _localBound(fq(U) / (fscale if _scale is None else _scale))  # what the optimizer minimizes

    optim = dict(_DEFAULT_OPTIM)
    if fitopts is not None and isinstance(fitopts.get("optimOptions"), dict):
        optim.update(fitopts["optimOptions"])

    # <============================================================================>
    # <============================ Start points ==================================>
    # <============================================================================>

    nStarts = max(0, int(math.floor(numstartpoints)))
    nJit = min(3, nStarts)

    def makeStarts():
        return _localFromT(_localMakeStarts(nStarts, nJit, LBt, UBt, z0T, free, rng), tf)

    starts = makeStarts()
    initialguess = starts

    maxAttempts = 3  # Initial attempt plus at most two retries.
    attempt = 0
    flagg = -math.inf
    fval = math.inf
    P = z0.copy()
    bestP = None
    bestFval = math.inf
    startAgreement = math.nan

    if nFree == 0:

        # Nothing to estimate: evaluate the fixed vector once.
        fval = f(z0)
        flagg = 1
        bestP = z0.copy()
        bestFval = fval
        startAgreement = 1.0

    else:

        # Retry with fresh start points when the multi-start reports failure, or when
        # its best point is a failed solve (runs that start where every solve
        # fails sit on a flat plateau and report convergence).
        while (flagg <= 0 or not math.isfinite(fval)) and attempt < maxAttempts:

            attempt += 1

            if attempt > 1:
                starts = makeStarts()
                initialguess = np.vstack([initialguess, starts])

            Pq, fvalq, flagg, allmins = _multistart(fqs, toQ(starts), nFree, optim)

            if Pq is None:
                fval = math.inf
            else:
                P = fromQ(Pq)
                fval = _localUnbound(fvalq) * fscale

            if allmins:
                fvs = np.array([_localUnbound(r.Fval) for r in allmins]) * fscale
                idx = int(np.argmin(fvs))
                fv = fvs[idx]
                if math.isfinite(fv) and fv < bestFval:
                    bestFval = float(fv)
                    bestP = fromQ(allmins[idx].X)
                okRun = np.array([r.Exitflag > 0 for r in allmins]) & np.isfinite(fvs)
                if np.any(okRun):
                    startAgreement = float(
                        np.sum(okRun & (fvs <= bestFval + 1e-4 * max(1, abs(bestFval)))) / np.sum(okRun)
                    )

    if flagg > 0 and not math.isfinite(fval):
        flagg = 0  # a "converged" failed solve is not a solution

    if flagg <= 0:
        if bestP is not None:
            warn(
                "GrowthPredict:fit_model:MultiStart",
                "The multi-start did not converge cleanly after %d attempts; using the best local minimum found (fval = %g)."
                % (attempt, bestFval),
            )
            P = bestP
            fval = bestFval
        else:
            raise GrowthPredictError(
                "GrowthPredict:fit_model:MaxAttemptsReached",
                "Fit failed after %d attempts (last exit flag: %g): no start point gave a valid objective value. "
                "Common causes: a first data point equal to zero with fixI0 = 1, or data the model cannot follow."
                % (attempt, flagg),
            )

    # Rescale-and-confirm pass. The objective was divided by |f(z0)|. When the seed
    # is far from the optimum, the scaled objective at the optimum is tiny and
    # the relative tolerances were much looser than intended (and the other
    # way round when the seed could not be evaluated). If the scale is off by more
    # than a factor of 100, run one more local solve from P with the scale reset.
    rescale = Struct(run=False, scaleBefore=fscale, scaleAfter=fscale, exitflag=math.nan, improved=False)
    fscaleNew = max(abs(fval), 1.0)
    if nFree > 0 and math.isfinite(fval) and (fscaleNew < 1e-2 * fscale or fscaleNew > 1e2 * fscale):
        rescale.run = True
        rescale.scaleAfter = fscaleNew
        try:
            qc, fc2, flc = _local_solve(lambda U: fqs(U, fscaleNew), toQ(P)[0], nFree, optim)
            fc2 = _localUnbound(fc2) * fscaleNew
            rescale.exitflag = flc
            if math.isfinite(fc2) and fc2 < fval:
                P = fromQ(qc)
                fval = fc2
                rescale.improved = True
        except Exception:
            pass  # keep the multi-start point
        fscale = fscaleNew

    # Final evaluation: fval, fitcurve and residual from ONE call at the returned P.
    # This call raises if the curve cannot be evaluated.
    optimizerFval = fval
    fval, fitcurve, _ = growthObjective(P, timevect, ydata, flag1, method1, full_output=True)

    residual = fitcurve - ydata

    tolB = 1e-3
    atLB = np.zeros(nz, dtype=bool)
    atUB = np.zeros(nz, dtype=bool)
    atLB[free] = P[free] <= LB[free] + tolB * np.abs(LB[free])
    atUB[free] = P[free] >= UB[free] - tolB * np.abs(UB[free])

    # Escape from a binding data-derived interval for r. The interval [r0/20, 20*r0]
    # comes from the slope of log incidence over the first 7 points; with small early
    # counts that slope can be far too small and the interval then excludes the
    # growth rate the rest of the window supports. When the fit lands on that bound,
    # refit once from the current solution with r allowed in [1e-3, 10] and keep the
    # result if it is better. The widened bounds are returned in fitDiagnostics, so
    # bootstrap refits use them too.
    if (
        not boundsSupplied
        and not rBoundsUser
        and flag1 != 5
        and (atLB[0] or atUB[0])
        and (LB[0] > 1e-3 or UB[0] < 10)
    ):

        wideopts = dict(fitopts) if fitopts is not None else {}
        wLB = LB.copy()
        wUB = UB.copy()
        wLB[0] = min(LB[0], 1e-3)
        wUB[0] = max(UB[0], 10)
        wideopts["LB"] = wLB
        wideopts["UB"] = wUB

        res2 = fit_model(
            data1, P, fixI0, numstartpoints, DT, flagX, forecastingperiod, wideopts, method1=method1, rng=rng
        )

        if res2.fval < fval:
            warn(
                "GrowthPredict:fit_model:RBoundsWidened",
                "The fit ended on the data-derived bound of r (%g in [%g, %g]); refitted with r in [%g, %g]: "
                "r = %g, objective %g -> %g." % (P[0], LB[0], UB[0], wLB[0], wUB[0], res2.P[0], fval, res2.fval),
            )
            diag2 = res2.fitDiagnostics
            diag2.rBoundsWidened = True
            return FitResult(
                res2.P,
                res2.residual,
                res2.fitcurve,
                res2.forecastcurve,
                res2.timevect2,
                np.vstack([initialguess, res2.initialguess]),
                res2.fval,
                res2.flagg,
                diag2,
            )

    fitDiagnostics = Struct(
        LB=LB,
        UB=UB,
        freeIndex=free,
        nFree=int(nFree),
        logscaleMask=tf,
        seed=z0,
        objectiveAtSeed=f0,
        objectiveScale=fscale,
        optimizerObjective=optimizerFval,
        finalObjective=fval,
        exitflag=flagg,
        attempts=attempt,
        startAgreement=startAgreement,
        numberOfStarts=int(starts.shape[0]),
        rescalePass=rescale,
        atLB=atLB,
        atUB=atUB,
        boundTolerance=tolB,
        boundsSupplied=bool(boundsSupplied),
        rBoundsWidened=False,
    )

    if forecastingperiod < 1:

        forecastcurve = fitcurve.copy()
        timevect2 = timevect.copy()

    else:

        timevect2 = np.arange(data1[0, 0], data1[-1, 0] + forecastingperiod + 0.5) * DT

        forecastcurve = growthIncidenceCurve(P, timevect2, flag1)

    return FitResult(P, residual, fitcurve, forecastcurve, timevect2, initialguess, fval, flagg, fitDiagnostics)
