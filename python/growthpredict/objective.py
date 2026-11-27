"""Objective of the growth-model fit.

Python port of ``growthObjectiveValue.m``, ``growthObjective.m`` and the legacy
wrapper ``plotModifiedLogisticGrowthMethods1.m``.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import gammaln

from .errors import GrowthPredictError
from .models import growthIncidenceCurve


def growthObjectiveValue(fitcurve, observations, method1, alpha=1.0, d=1.0):
    """Score a fitted incidence curve against the data.

    Port of ``growthObjectiveValue.m``.

    ``method1``
        0  sum of squared residuals
        1  Poisson negative log-likelihood
        3  negative binomial NLL, VAR = mean + alpha*mean
        4  negative binomial NLL, VAR = mean + alpha*mean^2
        5  negative binomial NLL, VAR = mean + alpha*mean^d

    The likelihoods leave out the term that depends only on the data,
    ``sum(log(y!))``; :func:`growthpredict.aicc.growthAICcCountAdjustment` adds
    it back for AICc. The negative binomial sums are evaluated with
    ``gammaln(y+m) - gammaln(m)``. Exact zeros of the fitted curve are replaced
    by 0.001. Invalid input or a non-finite score raises
    :class:`GrowthPredictError`; :func:`growthObjective` maps that to ``inf``.
    """
    yfit = np.asarray(fitcurve)
    ydata = np.asarray(observations)

    invalid = (
        yfit.size == 0
        or yfit.size != ydata.size
        or max(yfit.ndim, 1) > 2
        or max(ydata.ndim, 1) > 2
        or (yfit.ndim == 2 and min(yfit.shape) != 1)
        or (ydata.ndim == 2 and min(ydata.shape) != 1)
        or np.iscomplexobj(yfit)
        or np.iscomplexobj(ydata)
    )
    if not invalid:
        yfit = yfit.astype(float).ravel()
        ydata = ydata.astype(float).ravel()
        invalid = not (np.all(np.isfinite(yfit)) and np.all(np.isfinite(ydata)))
    if invalid:
        raise GrowthPredictError(
            "GrowthPredict:growthObjectiveValue:InvalidInput",
            "A score requires matching finite real observation and curve vectors.",
        )

    if np.any(yfit < 0):
        raise GrowthPredictError(
            "GrowthPredict:growthObjectiveValue:NegativeCurve",
            "The fitted incidence curve must be nonnegative.",
        )
    if np.sum(yfit) == 0:
        raise GrowthPredictError(
            "GrowthPredict:growthObjectiveValue:ZeroCurve",
            "The fitted curve is identically zero.",
        )

    yfit = np.where(yfit == 0, 0.001, yfit)  # allows the logarithms below
    yg = np.maximum(ydata, 0)  # argument of the gamma terms

    with np.errstate(all="ignore"):

        if method1 == 0:  # least squares

            objfunction = np.sum((ydata - yfit) ** 2)

        elif method1 == 1:  # Poisson

            objfunction = -np.sum(ydata * np.log(yfit) - yfit)

        elif method1 == 3:  # negative binomial, VAR = mean + alpha*mean

            m = yfit / alpha
            objfunction = -np.sum(
                gammaln(yg + m) - gammaln(m) + ydata * np.log(alpha) - (ydata + m) * np.log(1 + alpha)
            )

        elif method1 == 4:  # negative binomial, VAR = mean + alpha*mean^2

            m = 1 / alpha
            objfunction = -np.sum(
                gammaln(yg + m)
                - gammaln(m)
                + ydata * np.log(alpha * yfit)
                - (ydata + m) * np.log(1 + alpha * yfit)
            )

        elif method1 == 5:  # negative binomial, VAR = mean + alpha*mean^d

            m = (1 / alpha) * yfit ** (2 - d)
            objfunction = -np.sum(
                gammaln(yg + m)
                - gammaln(m)
                + ydata * np.log(alpha * yfit ** (d - 1))
                - (ydata + m) * np.log(1 + alpha * yfit ** (d - 1))
            )

        else:
            raise GrowthPredictError(
                "GrowthPredict:growthObjectiveValue:UnsupportedMethod",
                "method1 = %g is not supported; use 0, 1, 3, 4 or 5." % method1,
            )

    objfunction = float(objfunction)
    if not math.isfinite(objfunction):
        raise GrowthPredictError(
            "GrowthPredict:growthObjectiveValue:InvalidObjective",
            "The parameters do not give a finite real objective.",
        )

    return objfunction


def growthObjective(z, tvec, ydata, flag1, method1, full_output=False):
    """Objective of the growth-model fit, without global variables.

    Port of ``growthObjective.m``; ``z = [r, p, a, K, I0, alpha, d]``.

    With ``full_output=False`` (the optimizer interface, MATLAB's one-output
    call) a trial whose ODE solve fails or whose score is not finite returns
    ``inf``. With ``full_output=True`` (MATLAB's multi-output call, used for
    the final evaluation in ``fit_model``) the same failures raise, so a
    penalty value can never be saved as an accepted fit, and the function
    returns ``(objfunction, fitcurve, C)`` where ``fitcurve`` is exactly the
    curve that was scored.
    """
    try:
        fitcurve, C = growthIncidenceCurve(z, tvec, flag1, return_cumulative=True)
        objfunction = growthObjectiveValue(fitcurve, ydata, method1, z[5], z[6])
    except GrowthPredictError:
        if full_output:
            raise
        return math.inf
    except (FloatingPointError, OverflowError, ZeroDivisionError, ValueError):
        if full_output:
            raise
        return math.inf

    if full_output:
        return objfunction, fitcurve, C
    return objfunction


def plotModifiedLogisticGrowthMethods1(z, tvec, ydata, flag1, method1):
    """Legacy interface to the fitting objective (``plotModifiedLogisticGrowthMethods1.m``).

    The MATLAB function read ``flag1``, ``method1``, ``timevect`` and ``ydata``
    from global variables; the Python port takes them as arguments and forwards
    to :func:`growthObjective`. A trial that cannot be evaluated returns ``inf``.
    """
    return growthObjective(z, tvec, ydata, flag1, method1)
