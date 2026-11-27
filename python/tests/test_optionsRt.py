"""Regression tests for the generation-interval options (port of tests/test_optionsRt.m)."""

import math

import numpy as np
from scipy.stats import gamma

import growthpredict as gp
from growthpredict import get_Rt, options_Rt, options_Rt_cholera, options_Rt_mpox
from growthpredict.rt import RT_OPTIONS


def test_default_is_mpox_serial_interval():
    # U.S. case pairs by symptom onset: mean 8.5 days, SD 5.0 days, in weeks.
    type_GId1, mean_GI1, var_GI1 = options_Rt()

    assert type_GId1 == 1
    assert math.isclose(mean_GI1, 8.5 / 7, rel_tol=1e-12)
    assert math.isclose(var_GI1, (5.0 / 7) ** 2, rel_tol=1e-12)


def test_mpox_file_matches_default():
    assert options_Rt_mpox() == options_Rt()


def test_cholera_file_keeps_former_default():
    # Mean 5 days, SD 8 days, in weeks.
    type_GId1, mean_GI1, var_GI1 = options_Rt_cholera()

    assert type_GId1 == 1
    assert math.isclose(mean_GI1, 5 / 7, rel_tol=1e-12)
    assert math.isclose(var_GI1, (8 / 7) ** 2, rel_tol=1e-12)


def test_implied_gamma_distributions():
    # Shape = mean^2/variance. The mpox interval has almost no mass within one
    # day; the cholera interval, applied without its one-day shift, has 41%.
    _, m, v = options_Rt()
    b = v / m
    a = m / b

    assert math.isclose(a, 2.89, rel_tol=1e-10)
    assert abs(gamma.cdf(1 / 7, a, scale=b) - 0.0066) < 0.0005

    _, m, v = options_Rt_cholera()
    b = v / m
    a = m / b

    assert math.isclose(a, 25 / 64, rel_tol=1e-10)
    assert abs(gamma.cdf(1 / 7, a, scale=b) - 0.4072) < 0.0005


def test_Rt_with_default_on_exponential_growth():
    # t = 0..19, incidence 100*exp(0.3*t). Exact R = (1 + r*b)^a = 1.4093.
    type_GId1, mean_GI1, var_GI1 = options_Rt()

    t = np.arange(20.0)
    I = 100 * np.exp(0.3 * t)

    Rs, _, _ = get_Rt(t, I, type_GId1, mean_GI1, var_GI1)

    b = var_GI1 / mean_GI1
    a = mean_GI1 / b

    np.testing.assert_allclose(Rs[-1], (1 + 0.3 * b) ** a, rtol=0.01)
    np.testing.assert_allclose(Rs[499], 1.4006215372, rtol=1e-6)  # reference
    np.testing.assert_allclose(Rs[-1], 1.4004398175, rtol=1e-6)  # reference

    # The cholera interval gives a value much closer to 1 for the same curve.
    Rs, _, _ = get_Rt(t, I, *options_Rt_cholera())

    np.testing.assert_allclose(Rs[-1], 1.1822069459, rtol=1e-6)  # reference


def test_every_Rt_options_function_runs():
    # Each options_Rt* function of the package must return a usable setting.
    found = {name: getattr(gp, name) for name in dir(gp) if name.startswith("options_Rt")}

    assert len(found) >= 3
    assert set(found) == set(RT_OPTIONS)

    for name, optionsFunction in found.items():
        type_GId1, mean_GI1, var_GI1 = optionsFunction()

        assert type_GId1 in (1, 2, 3), name
        assert math.isfinite(mean_GI1) and mean_GI1 > 0, name
        assert math.isfinite(var_GI1) and var_GI1 >= 0, name
        if type_GId1 == 1:
            assert var_GI1 > 0, name
