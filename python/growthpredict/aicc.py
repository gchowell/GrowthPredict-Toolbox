"""Model comparison: parameter counts and the corrected Akaike information criterion.

Python port of ``get_nparams.m``, ``getAICc.m`` and ``growthAICcCountAdjustment.m``.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import gammaln

from .errors import GrowthPredictError, warn
from .objective import growthObjectiveValue
from .utils import Struct


def get_nparams(method1, dist1, flag1, fixI0):
    """Number of estimated quantities of a growth-model fit (``get_nparams.m``).

    Model parameters (by ``flag1``) + the initial value when ``fixI0 = 0`` +
    the parameters of the error model. Least squares counts one parameter for
    the residual variance whatever ``dist1`` is. ``dist1`` is kept in the
    argument list for compatibility; it is not used.
    """
    # flag -1=EXP 0=GGM 1=GLM 2=GRM 3=Logistic 4=RIC 5=Gompertz
    counts = {-1: 1, 0: 2, 1: 3, 2: 4, 3: 2, 4: 3, 5: 2}
    if flag1 not in counts:
        raise GrowthPredictError(
            "GrowthPredict:get_nparams:UnknownModel",
            "flag1 = %g is not a supported growth model (-1 to 5)." % flag1,
        )
    numparams = counts[flag1]

    if fixI0 == 0:  # initial value is estimated
        numparams += 1

    if method1 == 0:  # least squares: residual variance
        numparams += 1
    elif method1 == 1:  # Poisson: no extra parameter
        pass
    elif method1 in (3, 4):  # negative binomial: alpha
        numparams += 1
    elif method1 == 5:  # negative binomial: alpha and d
        numparams += 2
    else:
        raise GrowthPredictError(
            "GrowthPredict:get_nparams:UnsupportedMethod",
            "method1 = %g is not supported; use 0, 1, 3, 4 or 5." % method1,
        )

    return numparams


def getAICc(method1, dist1, flag1, fixI0, fval, n, fvalAdjust=None):
    """Corrected Akaike information criterion of a growth-model fit (``getAICc.m``).

    Parameters
    ----------
    fval : float
        Objective value at the fit (sum of squares for ``method1 = 0``,
        negative log-likelihood without the ``sum(log(y!))`` term otherwise).
    n : int
        Number of observations in the calibration window.
    fvalAdjust : float, optional
        (Count likelihoods) amount that completes ``fval`` to the full negative
        log-likelihood, from :func:`growthAICcCountAdjustment`. Omitted -> 0.

    Returns
    -------
    AICc, part1, part2, numparams, nEff
        ``AICc = part1 + part2``, ``part1 = -2 log L``,
        ``part2 = 2k + 2k(k+1)/(nEff-k-1)``. Least squares uses the complete
        Gaussian log-likelihood,
        ``-2 log L = nEff*log(SSE/nEff) + nEff*log(2*pi) + nEff``.
        ``nEff = n - 1`` when ``fixI0 = 1``: the first fitted value equals the
        first observation by construction. With too few data points AICc is
        ``inf`` with a warning.
    """
    if fvalAdjust is None or np.size(fvalAdjust) == 0:
        fvalAdjust = 0

    numparams = get_nparams(method1, dist1, flag1, fixI0)

    nEff = n
    if fixI0 is not None and fixI0 == 1:
        nEff = n - 1

    if nEff <= numparams + 1:

        warn(
            "GrowthPredict:getAICc:InsufficientSample",
            "AICc is not defined for n=%g informative observations and k=%g; "
            "returning Inf (requires n>k+1)." % (nEff, numparams),
        )
        return math.inf, math.nan, math.inf, numparams, nEff

    part2 = 2 * numparams + (2 * numparams * (numparams + 1)) / (nEff - numparams - 1)

    if method1 == 0:
        with np.errstate(all="ignore"):
            part1 = float(nEff * np.log(np.float64(fval) / nEff) + nEff * np.log(2 * np.pi) + nEff)
    else:
        part1 = 2 * (fval + fvalAdjust)

    AICc = part1 + part2

    return AICc, part1, part2, numparams, nEff


def growthAICcCountAdjustment(method1, ydata, fixI0=None, alpha=None, d=None):
    """Complete the count objective to a full log-likelihood for AICc.

    Port of ``growthAICcCountAdjustment.m``. For the Poisson and negative
    binomial methods (``method1 = 1, 3, 4, 5``) the fitting objective is the
    negative log-likelihood WITHOUT ``sum(log(y!))``. The returned ``adj`` is
    the amount to add so that ``fval + adj`` is the full negative
    log-likelihood of the observations that count toward the sample size:

        adj = sum(gammaln(y + 1)) over the informative observations
              - (objective of the fixed first point)            [fixI0 = 1 only]

    ``adj`` is 0 for ``method1 = 0``. Returns ``(adj, info)``.
    """
    ydata = np.asarray(ydata, dtype=float).ravel()
    adj = 0.0
    info = Struct(logFactorialTerm=0.0, firstPointTerm=0.0, nInformative=int(ydata.size))

    if method1 not in (1, 3, 4, 5):
        return adj, info
    if alpha is None or np.size(alpha) == 0:
        alpha = 1
    if d is None or np.size(d) == 0:
        d = 1

    y = ydata

    pinned = np.zeros(y.size, dtype=bool)
    if fixI0 is not None and fixI0 == 1 and y.size > 0:
        pinned[0] = True

    logFact = float(np.sum(gammaln(np.maximum(y[~pinned], 0) + 1)))

    firstTerm = 0.0
    if np.any(pinned):
        y1 = y[pinned]
        mu1 = y1.copy()
        mu1[mu1 == 0] = 0.001  # the objective's own replacement for a zero mean
        firstTerm = growthObjectiveValue(mu1, y1, method1, alpha, d)

    adj = logFact - firstTerm
    info = Struct(
        logFactorialTerm=logFact, firstPointTerm=firstTerm, nInformative=int(np.count_nonzero(~pinned))
    )
    return adj, info
