"""Regression tests for the interval of r in fit_model (port of tests/test_rInterval.m).

The series is the bundled U.S. mpox file, GLM, least squares, fixI0 = 1.
Reference values are those of the MATLAB test file, which reproduce the README
quick-start estimates. Tolerances are loose enough for differences between
optimizers.
"""

import numpy as np
import pytest

from growthpredict import GrowthPredictError, fit_model, initialParams


def local_fit(mpox, firstRow, windowsize1, flag1, fitopts):
    data1 = mpox[firstRow - 1 : firstRow - 1 + windowsize1, :]
    params0 = initialParams(data1[:, 1], flag1)
    res = fit_model(data1, params0, 1, 10, 1, flag1, 0, fitopts, method1=0)
    return res.P, res.fval, res.fitDiagnostics


def test_quick_start_window_is_unchanged(mpox):
    # Rows 1-20: the README quick-start fit.
    P, fval, diagnostics = local_fit(mpox, 1, 20, 1, None)

    np.testing.assert_allclose(P[0], 1.7302, rtol=0.02)
    assert abs(P[1] - 0.8526) < 0.005
    np.testing.assert_allclose(P[3], 28547, rtol=0.01)
    np.testing.assert_allclose(fval, 3.1779e5, rtol=0.01)
    assert not (diagnostics.atLB[0] or diagnostics.atUB[0])


def test_declining_window_is_not_limited_by_r(mpox):
    # Rows 13-32 (README: 'the last 20 observations'). With r capped at 10 the
    # sum of squares was 9.51e5; the data support r of about 800 and 2.32e4.
    P, fval, diagnostics = local_fit(mpox, 13, 20, 1, None)

    assert fval < 3.0e4
    assert P[0] > 10
    assert not (diagnostics.atLB[0] or diagnostics.atUB[0])
    assert not diagnostics.rBoundsWidened

    # MATLAB run of the maintainer: r, p, K, sum of squares = 801.8, 0.178, 18153, 2.323e4
    np.testing.assert_allclose(P[0], 801.8, rtol=0.02)
    assert abs(P[1] - 0.178) < 0.005
    np.testing.assert_allclose(P[3], 18153, rtol=0.005)
    np.testing.assert_allclose(fval, 2.323e4, rtol=0.01)


def test_mid_epidemic_window_is_not_limited_by_r(mpox):
    # Rows 8-27. With r capped at 10 the sum of squares was 5.49e5; reference 9.24e4.
    P, fval, diagnostics = local_fit(mpox, 8, 20, 1, None)

    assert fval < 1.0e5
    np.testing.assert_allclose(P[0], 26.97, rtol=0.10)
    assert not (diagnostics.atLB[0] or diagnostics.atUB[0])


def test_interval_scales_with_the_data(mpox):
    # GGM, GLM and GRM: r in [1e-4, 10*max(y)]/DT.
    ymax = np.max(mpox[12:32, 1])

    for flag1 in (0, 1, 2):
        _, _, diagnostics = local_fit(mpox, 13, 20, flag1, None)
        np.testing.assert_allclose(diagnostics.LB[0], 1e-4, rtol=1e-12)
        np.testing.assert_allclose(diagnostics.UB[0], 10 * ymax, rtol=1e-12)


def test_fixed_exponent_models_keep_their_interval(mpox):
    # Logistic model (p fixed at 1): r is a rate per unit time and its interval
    # still comes from the early slope of log incidence, [slope/20, 20*slope].
    b = np.polyfit(mpox[:7, 0], np.log(np.maximum(mpox[:7, 1], 0.5)), 1)

    _, _, diagnostics = local_fit(mpox, 1, 20, 3, None)

    np.testing.assert_allclose(diagnostics.LB[0], b[0] / 20, rtol=1e-10)
    np.testing.assert_allclose(diagnostics.UB[0], b[0] * 20, rtol=1e-10)


def test_user_interval_is_honoured(mpox):
    # fitopts['rBounds'] replaces the interval and switches off the retry.
    fitopts = {"rBounds": [1e-3, 10]}

    P, fval, diagnostics = local_fit(mpox, 13, 20, 1, fitopts)

    np.testing.assert_allclose(diagnostics.LB[0], 1e-3, rtol=1e-12)
    np.testing.assert_allclose(diagnostics.UB[0], 10, rtol=1e-12)
    assert P[0] <= 10 * (1 + 1e-9)
    assert diagnostics.atUB[0]
    assert not diagnostics.rBoundsWidened
    np.testing.assert_allclose(fval, 9.51e5, rtol=0.03)
    assert abs(P[1] - 0.676) < 0.01  # MATLAB run of the maintainer


def test_invalid_user_interval_is_rejected(mpox):
    data1 = mpox[:20, :]
    params0 = initialParams(data1[:, 1], 1)

    for bad in ([5, 1], [0, 1], [1, 2, 3], "abc"):
        with pytest.raises(GrowthPredictError) as e:
            fit_model(data1, params0, 1, 2, 1, 1, 0, {"rBounds": bad}, method1=0)
        assert e.value.identifier == "GrowthPredict:fit_model:rBounds"
