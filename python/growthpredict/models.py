"""Growth models and the single place where they are integrated.

Python port of ``modifiedLogisticGrowth.m``, ``growthODEOptions.m`` and
``growthIncidenceCurve.m``.

The model evolves a cumulative state C(t):

    flag1  name   dC/dt
    -----  -----  ------------------------
     -1    EXP    r*C
      0    GGM    r*C^p
      1    GLM    r*C^p*(1 - C/K)
      2    GRM    r*C^p*(1 - (C/K)^a)
      3    LM     r*C*(1 - C/K)
      4    RICH   r*C*(1 - (C/K)^a)
      5    GOM    r*C*exp(-a*t)      (t = time since the first point)
"""

from __future__ import annotations

import math
import warnings

import numpy as np
from scipy.integrate import odeint

from .errors import GrowthPredictError

EXP, GGM, GLM, GRM, LM, RICH, GOM = -1, 0, 1, 2, 3, 4, 5

MODEL_NAMES = {EXP: "EXP", GGM: "GGM", GLM: "GLM", GRM: "GRM", LM: "LM", RICH: "RICH", GOM: "GOM"}

_SUPPORTED_FLAGS = (-1, 0, 1, 2, 3, 4, 5)

# Message returned by scipy.integrate.odeint (LSODA) for a complete solve.
_ODE_SUCCESS = "Integration successful."


def modifiedLogisticGrowth(t, x, r, p, a, K, flag1):
    """Right-hand side dC/dt of the growth model selected by ``flag1``.

    Port of ``modifiedLogisticGrowth.m``. ``x`` is the cumulative count C
    (a scalar or a one-element sequence); the return value is a float. An
    unknown ``flag1`` returns 0, as in MATLAB (``dx=zeros(1,1)``).
    """
    try:
        x = x[0]
    except (TypeError, IndexError):
        pass
    x = float(x)

    if flag1 == -1:  # EXP
        return r * x
    if flag1 == 0:  # GGM
        return r * (x ** p)
    if flag1 == 1:  # GLM
        return r * (x ** p) * (1 - (x / K))
    if flag1 == 2:  # GRM
        return r * (x ** p) * (1 - (x / K) ** a)
    if flag1 == 3:  # Logistic
        return r * x * (1 - (x / K))
    if flag1 == 4:  # Richards
        return r * x * (1 - (x / K) ** a)
    if flag1 == 5:  # Gompertz; t is the time since the first point (growthIncidenceCurve)
        return r * x * math.exp(-a * t)
    return 0.0


def growthODEOptions():
    """One numerical policy for every ODE solve in GrowthPredict.

    Port of ``growthODEOptions.m``. The values are used by
    :func:`growthIncidenceCurve` for the fitting objective, the final fitted
    curve and the forecast, so the curve that is scored is the curve that is
    reported. The cumulative curve C(t) is constrained to be nonnegative.

    MATLAB solves with ``ode15s`` at RelTol 1e-6 / AbsTol 1e-8. This port uses
    LSODA (``scipy.integrate.odeint``), which switches between a non-stiff and
    a stiff method automatically, at RelTol 1e-8 / AbsTol 1e-10: at 1e-6 the
    sums of squares of the reference fits were off in the fifth digit, and the
    optimizer differentiates the objective by finite differences, so solver
    error has to stay well below those differences. The tighter tolerances
    cost about 10% in run time.
    """
    return {"RelTol": 1e-8, "AbsTol": 1e-10, "NonNegative": True, "MaxSteps": 5000}


def _rhs(C, t, r, p, a, K, flag1):
    """dC/dt for odeint, with C clipped at zero (the 'NonNegative' option)."""
    c = C[0]
    if not c > 0.0:  # also catches NaN
        c = 0.0
    try:
        if flag1 == 1:
            return (r * (c ** p) * (1.0 - c / K),)
        if flag1 == 0:
            return (r * (c ** p),)
        if flag1 == -1:
            return (r * c,)
        if flag1 == 2:
            return (r * (c ** p) * (1.0 - (c / K) ** a),)
        if flag1 == 3:
            return (r * c * (1.0 - c / K),)
        if flag1 == 4:
            return (r * c * (1.0 - (c / K) ** a),)
        if flag1 == 5:
            return (r * c * math.exp(-a * t),)
    except (OverflowError, ZeroDivisionError, ValueError):
        # The trajectory left the range of double precision: let the solver fail.
        return (math.inf,)
    return (0.0,)


def growthIncidenceCurve(z, tvec, flag1, return_cumulative=False):
    """Solve the growth model and return the incidence curve.

    Port of ``growthIncidenceCurve.m``.

    Parameters
    ----------
    z : sequence
        Parameter vector; ``z[0:5] = [r, p, a, K, I0]`` (longer vectors are
        accepted, the error-model entries alpha and d are ignored).
    tvec : sequence
        Time points (at least 3, increasing).
    flag1 : int
        Growth-model code (-1 EXP, 0 GGM, 1 GLM, 2 GRM, 3 LM, 4 RICH, 5 GOM).
    return_cumulative : bool
        When True, return ``(incidence, C)`` like the two-output MATLAB call.

    Returns
    -------
    incidence : ndarray
        ``abs([C(1); diff(C)])``, the length of ``tvec``.
    C : ndarray
        The cumulative curve (only with ``return_cumulative=True``).

    Raises
    ------
    GrowthPredictError
        ``GrowthPredict:growthIncidenceCurve:TooFewTimes`` for fewer than three
        time points, ``...:InvalidTrajectory`` when the solve is incomplete or
        returns non-finite values. :func:`growthpredict.objective.growthObjective`
        turns that into ``inf`` for the optimizer, so a failed solve can never
        be saved as a fitted curve.

    Notes
    -----
    Gompertz model. Its growth rate ``r*exp(-a*t)`` depends on time, unlike the
    other models, so the time origin matters. Time is measured from
    ``tvec[0]``: ``r`` is the growth rate at the first point of the curve,
    whatever its time label.
    """
    tvec = np.asarray(tvec, dtype=float).ravel()
    if tvec.size < 3:
        raise GrowthPredictError(
            "GrowthPredict:growthIncidenceCurve:TooFewTimes",
            "At least 3 time points are needed (got %d)." % tvec.size,
        )

    if flag1 == 5:
        tvec = tvec - tvec[0]  # elapsed time; the other models do not depend on t

    r, p, a, K, I0 = (float(z[0]), float(z[1]), float(z[2]), float(z[3]), float(z[4]))

    if flag1 not in _SUPPORTED_FLAGS:
        # MATLAB: dx stays zero for an unknown flag, i.e. a constant curve.
        C = np.full(tvec.size, I0, dtype=float)
    else:
        opts = growthODEOptions()
        ok = False
        C = None
        if math.isfinite(I0) and all(math.isfinite(v) for v in (r, p, a, K)):
            with warnings.catch_warnings(), np.errstate(all="ignore"):
                warnings.simplefilter("ignore")
                try:
                    sol, info = odeint(
                        _rhs,
                        [I0],
                        tvec,
                        args=(r, p, a, K, flag1),
                        rtol=opts["RelTol"],
                        atol=opts["AbsTol"],
                        mxstep=opts["MaxSteps"],
                        full_output=True,
                    )
                    ok = info.get("message", "") == _ODE_SUCCESS
                    C = sol[:, 0]
                except Exception:
                    ok = False
        if not ok or C is None or C.shape[0] != tvec.size or not np.all(np.isfinite(C)):
            raise GrowthPredictError(
                "GrowthPredict:growthIncidenceCurve:InvalidTrajectory",
                "The ODE solve must return complete finite real values at the requested times.",
            )
        if opts["NonNegative"]:
            C = np.maximum(C, 0.0)

    incidence = np.abs(np.concatenate(([C[0]], np.diff(C))))

    if return_cumulative:
        return incidence, C
    return incidence
