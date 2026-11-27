"""Tests of the growth models, the objective, AICc and the parameter counts.

These functions have no MATLAB test file; the targets are closed-form
solutions and likelihoods evaluated independently with scipy.stats.
"""

import math
import warnings

import numpy as np
import pytest
from scipy.special import gammaln
from scipy.stats import nbinom, poisson

import growthpredict as gp
from growthpredict import (
    GrowthPredictError,
    get_nparams,
    getAICc,
    growthAICcCountAdjustment,
    growthIncidenceCurve,
    growthObjective,
    growthObjectiveValue,
    initialParams,
    modifiedLogisticGrowth,
)

T = np.arange(30.0)


def incidence_from(C):
    return np.abs(np.concatenate(([C[0]], np.diff(C))))


# ---------------------------------------------------------------- models


def test_exponential_closed_form():
    r, C0 = 0.21, 3.0
    got = growthIncidenceCurve([r, 1, 1, 0, C0], T, -1)
    np.testing.assert_allclose(got, incidence_from(C0 * np.exp(r * T)), rtol=1e-6)


def test_generalized_growth_closed_form():
    # dC/dt = r*C^p  ->  C = (C0^(1-p) + (1-p)*r*t)^(1/(1-p))
    r, p, C0 = 1.3, 0.6, 5.0
    C = (C0 ** (1 - p) + (1 - p) * r * T) ** (1 / (1 - p))
    got = growthIncidenceCurve([r, p, 1, 0, C0], T, 0)
    np.testing.assert_allclose(got, incidence_from(C), rtol=1e-6)


def test_logistic_closed_form_and_glm_with_p_one():
    r, K, C0 = 0.45, 12000.0, 7.0
    C = K / (1 + (K / C0 - 1) * np.exp(-r * T))
    expected = incidence_from(C)
    np.testing.assert_allclose(growthIncidenceCurve([r, 1, 1, K, C0], T, 3), expected, rtol=1e-6)
    # GLM with p = 1 and GRM / Richards with a = 1 are the logistic model
    np.testing.assert_allclose(growthIncidenceCurve([r, 1, 1, K, C0], T, 1), expected, rtol=1e-6)
    np.testing.assert_allclose(growthIncidenceCurve([r, 1, 1, K, C0], T, 2), expected, rtol=1e-6)
    np.testing.assert_allclose(growthIncidenceCurve([r, 1, 1, K, C0], T, 4), expected, rtol=1e-6)


def test_richards_closed_form():
    # dC/dt = r*C*(1-(C/K)^a)  ->  C = K / (1 + Q*exp(-a*r*t))^(1/a), Q = (K/C0)^a - 1
    r, a, K, C0 = 0.4, 0.6, 9000.0, 11.0
    Q = (K / C0) ** a - 1
    C = K / (1 + Q * np.exp(-a * r * T)) ** (1 / a)
    np.testing.assert_allclose(growthIncidenceCurve([r, 1, a, K, C0], T, 4), incidence_from(C), rtol=1e-6)
    # GRM with p = 1 is the Richards model
    np.testing.assert_allclose(growthIncidenceCurve([r, 1, a, K, C0], T, 2), incidence_from(C), rtol=1e-6)


def test_curve_starts_at_initial_value_and_saturates():
    inc, C = growthIncidenceCurve([0.9, 0.8, 1, 5000, 20], np.arange(200.0), 1, return_cumulative=True)
    assert inc[0] == 20
    assert C[0] == 20
    assert np.all(np.diff(C) >= -1e-6)
    np.testing.assert_allclose(C[-1], 5000, rtol=1e-6)
    np.testing.assert_allclose(np.sum(inc), C[-1], rtol=1e-8)  # abs() of rounding-size decreases at the plateau


def test_right_hand_side():
    assert modifiedLogisticGrowth(0, 10, 0.5, 1, 1, 100, -1) == 5
    assert math.isclose(modifiedLogisticGrowth(0, 16, 0.5, 0.5, 1, 100, 0), 2)
    assert math.isclose(modifiedLogisticGrowth(0, 16, 0.5, 0.5, 1, 32, 1), 1)
    assert math.isclose(modifiedLogisticGrowth(0, 16, 0.5, 0.5, 2, 32, 2), 2 * 0.75)
    assert math.isclose(modifiedLogisticGrowth(0, 16, 0.5, 1, 1, 32, 3), 4)
    assert math.isclose(modifiedLogisticGrowth(0, 16, 0.5, 1, 2, 32, 4), 6)
    assert math.isclose(modifiedLogisticGrowth(2, 10, 0.5, 1, 0.3, 1, 5), 5 * math.exp(-0.6))
    assert modifiedLogisticGrowth(0, [10], 0.5, 1, 1, 100, -1) == 5


def test_incidence_curve_errors():
    with pytest.raises(GrowthPredictError) as e:
        growthIncidenceCurve([1, 1, 1, 1, 1], [0, 1], 1)
    assert e.value.identifier == "GrowthPredict:growthIncidenceCurve:TooFewTimes"

    # a trajectory that leaves the range of double precision is an error, not a curve
    with pytest.raises(GrowthPredictError) as e:
        growthIncidenceCurve([3e4, 1, 1, 0, 14], np.arange(20.0), 0)
    assert e.value.identifier == "GrowthPredict:growthIncidenceCurve:InvalidTrajectory"

    with pytest.raises(GrowthPredictError):
        growthIncidenceCurve([float("nan"), 1, 1, 1, 1], np.arange(5.0), 1)


def test_initial_params():
    y = np.array([4.0, 6, 10, 20])
    s = 40.0
    np.testing.assert_allclose(initialParams(y, -1), [1, 1, 1, s, 4, 1, 1])
    np.testing.assert_allclose(initialParams(y, 0), [1, 0.9, 1, s, 4, 1, 1])
    np.testing.assert_allclose(initialParams(y, 1), [0.1, 0.9, 1, s, 4, 1, 1])
    np.testing.assert_allclose(initialParams(y, 2), [0.1, 0.9, 1, s, 4, 1, 1])
    np.testing.assert_allclose(initialParams(y, 3), [0.1, 1, 1, s, 4, 1, 1])
    np.testing.assert_allclose(initialParams(y, 4), [0.1, 1, 0.9, s, 4, 1, 1])
    r_guess = 1 - 4 / s
    np.testing.assert_allclose(initialParams(y, 5), [r_guess + 0.001, 1, r_guess / math.log(s / 4), s, 4, 1, 1])


# ---------------------------------------------------------------- objective

Y = np.array([3.0, 0, 7, 12, 30, 41, 18, 5])
MU = np.array([2.5, 1.2, 6.0, 14.0, 27.0, 44.0, 20.0, 4.0])
LOGFACT = float(np.sum(gammaln(Y + 1)))


def test_least_squares_objective():
    assert math.isclose(growthObjectiveValue(MU, Y, 0), float(np.sum((Y - MU) ** 2)))


def test_poisson_objective_is_the_likelihood_without_log_factorial():
    full = -float(np.sum(poisson.logpmf(Y, MU)))
    assert math.isclose(growthObjectiveValue(MU, Y, 1) + LOGFACT, full, rel_tol=1e-12)


@pytest.mark.parametrize("method1,alpha,d", [(3, 0.7, 1), (4, 0.05, 1), (5, 0.02, 1.8), (5, 0.3, 0.7)])
def test_negative_binomial_objectives_match_scipy(method1, alpha, d):
    if method1 == 3:
        var = MU + alpha * MU
    elif method1 == 4:
        var = MU + alpha * MU**2
    else:
        var = MU + alpha * MU**d
    p = MU / var
    n = MU * p / (1 - p)
    full = -float(np.sum(nbinom.logpmf(Y, n, p)))
    assert math.isclose(growthObjectiveValue(MU, Y, method1, alpha, d) + LOGFACT, full, rel_tol=1e-10)


def test_zero_mean_is_replaced_and_invalid_input_is_rejected():
    mu = MU.copy()
    mu[1] = 0
    ref = MU.copy()
    ref[1] = 0.001
    assert math.isclose(growthObjectiveValue(mu, Y, 1), growthObjectiveValue(ref, Y, 1))

    for args, ident in (
        ((MU[:-1], Y, 0), "GrowthPredict:growthObjectiveValue:InvalidInput"),
        (([], [], 0), "GrowthPredict:growthObjectiveValue:InvalidInput"),
        ((np.where(MU > 40, np.nan, MU), Y, 0), "GrowthPredict:growthObjectiveValue:InvalidInput"),
        ((-MU, Y, 0), "GrowthPredict:growthObjectiveValue:NegativeCurve"),
        ((np.zeros_like(MU), Y, 0), "GrowthPredict:growthObjectiveValue:ZeroCurve"),
        ((MU, Y, 2), "GrowthPredict:growthObjectiveValue:UnsupportedMethod"),
        ((MU, Y, 3, 0.0), "GrowthPredict:growthObjectiveValue:InvalidObjective"),
    ):
        with pytest.raises(GrowthPredictError) as e:
            growthObjectiveValue(*args)
        assert e.value.identifier == ident


def test_objective_returns_inf_for_failed_solves_and_raises_on_request():
    t = np.arange(20.0)
    y = np.linspace(1, 50, 20)
    z_bad = [3e4, 1, 1, 0, 14, 0, 0]
    assert growthObjective(z_bad, t, y, 0, 0) == math.inf
    with pytest.raises(GrowthPredictError):
        growthObjective(z_bad, t, y, 0, 0, full_output=True)

    z = [0.3, 0.9, 1, 500, 1, 0, 0]
    value, fitcurve, C = growthObjective(z, t, y, 1, 0, full_output=True)
    assert math.isclose(value, float(np.sum((y - fitcurve) ** 2)))
    assert math.isclose(value, growthObjective(z, t, y, 1, 0))
    assert math.isclose(value, gp.plotModifiedLogisticGrowthMethods1(z, t, y, 1, 0))
    np.testing.assert_allclose(np.cumsum(fitcurve), C)


# ---------------------------------------------------------------- parameter counts and AICc


def test_number_of_parameters():
    model = {-1: 1, 0: 2, 1: 3, 2: 4, 3: 2, 4: 3, 5: 2}
    error = {0: 1, 1: 0, 3: 1, 4: 1, 5: 2}
    for flag1, km in model.items():
        for method1, ke in error.items():
            assert get_nparams(method1, method1, flag1, 1) == km + ke
            assert get_nparams(method1, method1, flag1, 0) == km + ke + 1
    # least squares counts the residual variance whatever dist1 is
    assert get_nparams(0, 0, 1, 1) == get_nparams(0, 1, 1, 1) == get_nparams(0, 2, 1, 1) == 4

    with pytest.raises(GrowthPredictError) as e:
        get_nparams(0, 0, 9, 1)
    assert e.value.identifier == "GrowthPredict:get_nparams:UnknownModel"
    with pytest.raises(GrowthPredictError) as e:
        get_nparams(2, 0, 1, 1)
    assert e.value.identifier == "GrowthPredict:get_nparams:UnsupportedMethod"


def test_aicc_least_squares_is_the_gaussian_likelihood():
    n, sse = 20, 3.1779e5
    AICc, part1, part2, k, nEff = getAICc(0, 0, 1, 1, sse, n)
    assert (k, nEff) == (4, 19)
    assert math.isclose(part1, 19 * math.log(sse / 19) + 19 * math.log(2 * math.pi) + 19)
    assert math.isclose(part2, 2 * 4 + (2 * 4 * 5) / (19 - 4 - 1))
    assert math.isclose(AICc, part1 + part2)

    # the initial value is estimated: all n points count and k grows by one
    _, _, _, k, nEff = getAICc(0, 0, 1, 0, sse, n)
    assert (k, nEff) == (5, 20)


def test_aicc_needs_more_points_than_parameters():
    # README: the default GLM fit has k = 4, so it needs at least seven observations
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert getAICc(0, 0, 1, 1, 10.0, 6)[0] == math.inf
    assert math.isfinite(getAICc(0, 0, 1, 1, 10.0, 7)[0])
    with pytest.warns(gp.GrowthPredictWarning, match="InsufficientSample"):
        getAICc(0, 0, 1, 1, 10.0, 6)


@pytest.mark.parametrize("method1,alpha,d", [(1, 1, 1), (3, 0.7, 1), (4, 0.05, 1), (5, 0.02, 1.8)])
def test_count_adjustment_completes_the_likelihood(method1, alpha, d):
    # With fixI0 = 1 the first fitted value equals the first observation.
    mu = MU.copy()
    mu[0] = Y[0]
    fval = growthObjectiveValue(mu, Y, method1, alpha, d)

    adj, info = growthAICcCountAdjustment(method1, Y, 1, alpha, d)

    if method1 == 1:
        full_rest = -float(np.sum(poisson.logpmf(Y[1:], mu[1:])))
    else:
        var = mu + alpha * mu ** (1 if method1 == 3 else (2 if method1 == 4 else d))
        p = mu / var
        nn = mu * p / (1 - p)
        full_rest = -float(np.sum(nbinom.logpmf(Y[1:], nn[1:], p[1:])))

    assert math.isclose(fval + adj, full_rest, rel_tol=1e-10)
    assert info.nInformative == Y.size - 1

    # fixI0 = 0: every observation counts
    adj0, info0 = growthAICcCountAdjustment(method1, Y, 0, alpha, d)
    assert math.isclose(adj0, LOGFACT)
    assert info0.nInformative == Y.size

    AICc, part1, _, _, _ = getAICc(method1, method1, 1, 1, fval, Y.size, adj)
    assert math.isclose(part1, 2 * full_rest, rel_tol=1e-10)


def test_count_adjustment_is_zero_for_least_squares():
    adj, info = growthAICcCountAdjustment(0, Y, 1)
    assert adj == 0 and info.logFactorialTerm == 0
