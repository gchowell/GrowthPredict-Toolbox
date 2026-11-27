"""Tests of the optional settings and input checks of fit_model (no MATLAB test file exists for these)."""

import math
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pytest

import growthpredict as gp
from growthpredict import GrowthPredictError, fit_model, getDoublingTimeCurve, initialParams


@pytest.fixture()
def window(mpox):
    data1 = mpox[:20, :]
    return data1, initialParams(data1[:, 1], 1)


def test_outputs_have_the_matlab_order_and_shapes(window):
    data1, params0 = window
    result = fit_model(data1, params0, 1, 4, 1, 1, 3, None, method1=0)
    P, residual, fitcurve, forecastcurve, timevect2, initialguess, fval, flagg, fitDiagnostics = result

    assert P.shape == (7,) and residual.shape == (20,) and fitcurve.shape == (20,)
    assert forecastcurve.shape == (23,) and timevect2.tolist() == list(range(23))
    np.testing.assert_allclose(residual, fitcurve - data1[:, 1])
    assert math.isclose(fval, float(np.sum(residual**2)))
    np.testing.assert_allclose(forecastcurve[:20], fitcurve, rtol=1e-6)
    assert flagg > 0
    # 4 box-wide starts, 3 jittered copies of the seed and the seed, in original units
    assert initialguess.shape == (8, 7)
    assert np.all(initialguess >= fitDiagnostics.LB - 1e-9) and np.all(initialguess <= fitDiagnostics.UB * (1 + 1e-9))
    assert fitDiagnostics.freeIndex.tolist() == [0, 1, 3] and fitDiagnostics.nFree == 3
    assert fitDiagnostics.logscaleMask.tolist() == [True, False, False, True, False, False, False]
    assert fitDiagnostics.numberOfStarts == 8 and fitDiagnostics.attempts == 1
    assert math.isclose(fitDiagnostics.finalObjective, fval)
    assert result.P is P and result.fitDiagnostics is fitDiagnostics


def test_supplied_box_is_used_as_given(window):
    data1, _ = window
    LB = [1.0, 0.5, 1, 20000, 14, 0, 0]
    UB = [3.0, 1.0, 1, 40000, 14, 0, 0]
    seed = [2.0, 0.9, 1, 30000, 14, 0, 0]
    res = fit_model(data1, seed, 1, 1, 1, 1, 0, {"LB": LB, "UB": UB}, method1=0)
    assert res.fitDiagnostics.boundsSupplied
    np.testing.assert_allclose(res.fitDiagnostics.LB, LB)
    np.testing.assert_allclose(res.fitDiagnostics.seed, seed)
    np.testing.assert_allclose(res.P[0], 1.7302, rtol=0.02)
    assert res.initialguess.shape == (3, 7)  # one box-wide start, one jittered seed, the seed

    # a box that excludes the optimum: the estimate sits on the bound and no retry widens it
    res = fit_model(data1, seed, 1, 1, 1, 1, 0, {"LB": [1.9] + LB[1:], "UB": UB}, method1=0)
    assert res.fitDiagnostics.atLB[0] and not res.fitDiagnostics.rBoundsWidened
    assert math.isclose(res.P[0], 1.9, rel_tol=1e-3)

    # everything fixed: the vector is evaluated once
    res = fit_model(data1, seed, 1, 5, 1, 1, 0, {"LB": seed, "UB": seed}, method1=0)
    assert res.fitDiagnostics.nFree == 0 and res.flagg == 1 and res.fitDiagnostics.startAgreement == 1
    np.testing.assert_allclose(res.P, seed)

    with pytest.raises(GrowthPredictError) as e:
        fit_model(data1, seed, 1, 1, 1, 1, 0, {"LB": LB[:5], "UB": UB[:5]}, method1=0)
    assert e.value.identifier == "GrowthPredict:fit_model:SuppliedBounds"
    with pytest.raises(GrowthPredictError) as e:
        fit_model(data1, seed, 1, 1, 1, 1, 0, {"LB": UB, "UB": LB}, method1=0)
    assert e.value.identifier == "GrowthPredict:fit_model:Bounds"


def test_optimizer_settings(window):
    data1, params0 = window

    # linear coordinates reach the same optimum
    res = fit_model(data1, params0, 1, 10, 1, 1, 0, {"logscale": False}, method1=0)
    assert not np.any(res.fitDiagnostics.logscaleMask)
    np.testing.assert_allclose(res.fval, 3.1779e5, rtol=0.01)

    # a tight budget still returns the best point found
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = fit_model(data1, params0, 1, 2, 1, 1, 0, {"optimOptions": {"MaxFunctionEvaluations": 30}}, method1=0)
    assert math.isfinite(res.fval) and res.fval >= 3.17e5

    # parallel solves are accepted and run serially
    import growthpredict.fitting as fitting

    fitting._parallel_warned = False
    with pytest.warns(gp.GrowthPredictWarning, match="GrowthPredict:fit_model:Parallel"):
        fit_model(data1, params0, 1, 1, 1, 1, 0, {"useParallel": True}, method1=0)


def test_input_checks(window):
    data1, params0 = window
    bad = data1.copy()
    bad[3, 1] = np.nan
    for args, kwargs, ident in (
        ((bad, params0, 1, 2, 1, 1, 0), dict(method1=0), "GrowthPredict:fit_model:InvalidData"),
        ((data1, params0, 1, 2, 1, 1, 0), dict(method1=2), "GrowthPredict:fit_model:UnsupportedMethod"),
        ((data1, params0, 1, 2, 1, 1, 0), dict(method1=None), "GrowthPredict:fit_model:NoMethod"),
        ((data1, params0, 1, 2, 1, 7, 0), dict(method1=0), "GrowthPredict:fit_model:UnknownModel"),
    ):
        with pytest.raises(GrowthPredictError) as e:
            fit_model(*args, **kwargs)
        assert e.value.identifier == ident


def test_negative_binomial_bounds_of_the_error_model(window):
    data1, params0 = window
    for method1, ub_alpha, d_bounds in ((3, 1e3, (1, 1)), (4, 1e5, (1, 1)), (5, 1e5, (0.6, 1e3))):
        res = fit_model(data1, params0, 1, 3, 1, 1, 0, None, method1=method1)
        D = res.fitDiagnostics
        assert D.LB[5] == 1e-8 and D.UB[5] == ub_alpha and (D.LB[6], D.UB[6]) == d_bounds
        assert D.logscaleMask[5]  # alpha is searched on the log scale
        assert D.LB[5] <= res.P[5] <= D.UB[5]


def test_estimated_initial_value_bounds(window):
    data1, params0 = window
    res = fit_model(data1, params0, 0, 3, 1, 1, 0, None, method1=0)
    D = res.fitDiagnostics
    # I0 in [max(1, 0.2*y1), 5*max(1, y1)] with y1 = 14
    assert math.isclose(D.LB[4], 2.8) and math.isclose(D.UB[4], 70)
    assert 4 in D.freeIndex.tolist()


def test_doubling_time_figure():
    plt.figure()
    tds, C0, curve, dt = getDoublingTimeCurve(np.exp(0.2 * np.arange(30.0)), 1, 1)
    axes = plt.gcf().get_axes()
    assert len(axes) == 2
    assert axes[0].get_ylabel() == "Cumulative cases" and axes[1].get_ylabel() == "Doubling time"
    assert curve.shape == (30, 2) and curve[:, 0].tolist() == list(range(30))
    np.testing.assert_allclose(curve[:, 1], np.cumsum(np.exp(0.2 * np.arange(30.0))), rtol=1e-12)
    assert tds[0].tolist() == [0.0, 1.0] and C0 == 1.0
