"""Regression tests for get_Rt (port of tests/test_getRt.m).

Reference values marked 'reference' are those of the MATLAB test file; they
come from an independent re-implementation of get_Rt.m. The other targets are
exact: for incidence growing as exp(r*t) the reproduction number is
(1 + r*b)^a for a gamma generation interval with shape a and scale b, and
exp(r*g) for a fixed interval g.
"""

import numpy as np
import pytest
from scipy.stats import gamma

from growthpredict import GrowthPredictError, get_Rt

# Cholera serial interval (options_Rt_cholera): mean 5 days, SD 8 days, in
# weeks. Its shape is below 1, the hardest case for the discretization.
MEAN_GI = 5 / 7
VAR_GI = (8 / 7) ** 2


def test_result_does_not_depend_on_first_time_label():
    # The same curve labelled 0..19, 1..20 and 12..31 must give the same values.
    t = np.arange(20.0)
    I = 50 * np.exp(0.25 * t) / (1 + np.exp(0.5 * (t - 10)))  # rises, then falls

    for settings in ((1, MEAN_GI, VAR_GI), (2, 1, 0), (3, 1, 0)):
        base, _, _ = get_Rt(t, I, *settings)
        for shift1 in (1, 12, 1000):
            shifted, _, _ = get_Rt(t + shift1, I, *settings)
            assert np.array_equal(np.isnan(shifted), np.isnan(base))
            ok = np.isfinite(base)
            np.testing.assert_allclose(shifted[ok], base[ok], rtol=1e-9)


def test_gamma_agrees_with_exact_value():
    # Exponential growth: exact R = (1 + r*b)^a. Weekly samples are interpolated
    # linearly, which costs a few tenths of a percent.
    b = VAR_GI / MEAN_GI
    a = MEAN_GI / b
    for r in (-0.2, 0.1, 0.3):
        for start1 in (0, 1, 12):
            t = start1 + np.arange(20.0)
            I = 100 * np.exp(r * np.arange(20.0))
            Rs, _, _ = get_Rt(t, I, 1, MEAN_GI, VAR_GI)
            np.testing.assert_allclose(Rs[-1], (1 + r * b) ** a, rtol=0.01)


def test_gamma_and_exponential_match_reference():
    t = np.arange(20.0)
    I = 100 * np.exp(0.3 * t)

    Rs, ts, _ = get_Rt(t, I, 1, MEAN_GI, VAR_GI)

    assert Rs.size == 1900
    assert abs(ts[500] - 5) < 1e-9
    np.testing.assert_allclose(Rs[499], 1.1850088598, rtol=1e-6)  # reference
    np.testing.assert_allclose(Rs[999], 1.1822348870, rtol=1e-6)  # reference
    np.testing.assert_allclose(Rs[-1], 1.1822069459, rtol=1e-6)  # reference

    Rs, _, _ = get_Rt(t, I, 2, 1, 0)

    np.testing.assert_allclose(Rs[499], 1.2945025419, rtol=1e-6)  # reference
    np.testing.assert_allclose(Rs[-1], 1.2925563361, rtol=1e-6)  # reference


def test_fixed_interval_is_exact():
    # R = I(t)/I(t-g) = exp(r*g).
    t = np.arange(20.0)
    I = 100 * np.exp(0.3 * t)

    Rs, _, _ = get_Rt(t, I, 3, 1, 0)

    np.testing.assert_allclose(Rs[-1], np.exp(0.3), rtol=1e-9)
    assert np.flatnonzero(np.isfinite(Rs))[0] == 99  # first estimate at elapsed time 1 (MATLAB index 100)

    Rs, _, _ = get_Rt(t, I, 3, 2, 0)

    np.testing.assert_allclose(Rs[-1], np.exp(0.6), rtol=1e-9)


def test_fixed_intervals_that_used_to_fail():
    # 0.14, 0.28, 0.29 and 0.57 are not exact multiples of the grid step in
    # floating point; 5/7 is not on the grid and is rounded to 0.71.
    t = np.arange(0, 19.005, 0.01)  # dense input: no interpolation error
    I = np.exp(0.3 * t)

    for g in (0.14, 0.28, 0.29, 0.57, 5 / 7):
        Rs, ts, _ = get_Rt(t, I, 3, g, 0)
        assert Rs.size == ts.size - 1
        np.testing.assert_allclose(Rs[-1], np.exp(0.3 * round(g / 0.01) * 0.01), rtol=1e-9)


def test_fixed_interval_limits():
    t = np.arange(20.0)
    I = 100 * np.exp(0.3 * t)

    # Longer than the curve: no lagged incidence, so no estimate.
    Rs, _, _ = get_Rt(t, I, 3, 25, 0)

    assert Rs.size == 1900
    assert np.all(np.isnan(Rs))

    # Shorter than half a grid step.
    with pytest.raises(GrowthPredictError) as e:
        get_Rt(t, I, 3, 0.004, 0)
    assert e.value.identifier == "GrowthPredict:get_Rt:lag"


def test_outputs_are_aligned():
    t = np.arange(12.0, 32.0)
    I = 100 * np.exp(0.1 * np.arange(20.0))

    Rs, ts, giCoverage = get_Rt(t, I, 1, MEAN_GI, VAR_GI)

    assert ts[0] == 12
    assert abs(ts[-1] - 31) < 1e-9
    assert Rs.shape == (ts.size - 1,)
    assert giCoverage.shape == Rs.shape

    # Coverage is the generation-interval distribution function of elapsed time.
    b = VAR_GI / MEAN_GI
    a = MEAN_GI / b

    assert np.min(np.diff(giCoverage)) >= 0
    assert abs(giCoverage[99] - gamma.cdf(1, a, scale=b)) < 1e-10
    assert abs(giCoverage[-1] - gamma.cdf(19, a, scale=b)) < 1e-10


def test_input_checks():
    t = np.arange(20.0)
    I = 100 * np.exp(0.1 * t)

    for args, ident in (
        ((4, 1, 1), "GrowthPredict:get_Rt:type"),
        ((1, -1, 1), "GrowthPredict:get_Rt:mean"),
        ((1, 1, 0), "GrowthPredict:get_Rt:variance"),
        ((1, 1), "GrowthPredict:get_Rt:variance"),
    ):
        with pytest.raises(GrowthPredictError) as e:
            get_Rt(t, I, *args)
        assert e.value.identifier == ident

    # A single observation gives no estimate.
    Rs, ts, cov = get_Rt([5], [10], 1, 1, 1)

    assert Rs.size == 0 and cov.size == 0
    assert ts.tolist() == [5]
