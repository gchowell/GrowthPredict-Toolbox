"""Regression tests for getDoublingTimeCurve (port of tests/test_doublingTime.m).

Reference values are those of the MATLAB test file. Times are found on a grid
of 0.01, so they are compared with a tolerance of one grid step.
"""

import numpy as np
import pytest

from growthpredict import GrowthPredictError, getDoublingTimeCurve

TOL = 0.011


def test_no_prior_cases_is_the_former_behaviour(mpox):
    # A window that starts at the first row of the series has no earlier cases:
    # three inputs, and a fourth input of zero, give the same result.
    y = mpox[:20, 1]

    tds3, C03, curve3, dt3 = getDoublingTimeCurve(y, 1, 0)
    tds4, C04, curve4, dt4 = getDoublingTimeCurve(y, 1, 0, 0)

    assert np.array_equal(tds4, tds3)
    assert C04 == C03
    assert np.array_equal(curve4, curve3)
    assert np.array_equal(dt4, dt3)

    assert C03 == 14
    np.testing.assert_allclose(dt3, [0.87, 0.89, 0.84, 0.93, 1.06, 1.12, 1.04, 1.17, 1.54, 2.40], atol=TOL)


def test_exponential_growth_inside_an_epidemic():
    # Incidence exp(0.2*k), k = 0..39, window k = 20..39. With the earlier cases
    # the count doubles every ln(2)/0.2 = 3.47 time units; when the count restarts
    # at the window, the first 'doubling' takes 0.81.
    I = np.exp(0.2 * np.arange(40.0))

    window1 = I[20:40]
    prior1 = float(np.sum(I[:20]))

    _, _, _, dtNew = getDoublingTimeCurve(window1, 1, 0, prior1)
    _, _, _, dtOld = getDoublingTimeCurve(window1, 1, 0)

    np.testing.assert_allclose(dtNew, [3.40, 3.46, 3.44, 3.47, 3.46], atol=TOL)
    assert np.max(np.abs(dtNew - np.log(2) / 0.2)) < 0.08

    np.testing.assert_allclose(dtOld[:3], [0.81, 1.34, 1.93], atol=TOL)

    # Same result as following the count of the whole series from the window start.
    _, _, _, dtWhole = getDoublingTimeCurve(np.concatenate(([np.sum(I[:21])], I[21:40])), 1, 0)

    np.testing.assert_allclose(dtNew, dtWhole, atol=1e-9)


def test_falling_incidence_has_no_doubling():
    # Logistic epidemic of 20000 cases, window starting at its midpoint: the count
    # goes from 10000 to 19999 and never doubles. Restarting the count at the
    # window reported doublings after 1.00 and 2.78 time units.
    k = np.arange(41.0)
    C = 20000 / (1 + np.exp(-0.5 * (k - 12)))
    I = np.concatenate(([C[0]], np.diff(C)))

    window1 = I[12:32]
    prior1 = float(np.sum(I[:12]))

    tds, C0data, _, dtNew = getDoublingTimeCurve(window1, 1, 0, prior1)
    _, _, _, dtOld = getDoublingTimeCurve(window1, 1, 0)

    assert dtNew.size == 0
    np.testing.assert_allclose(C0data, 10000, rtol=1e-9)
    assert tds.shape[0] == 1

    np.testing.assert_allclose(dtOld, [1.00, 2.78], atol=TOL)

    # A window in the rising phase still doubles, more slowly than its own running total.
    _, _, _, dtNew = getDoublingTimeCurve(I[6:26], 1, 0, float(np.sum(I[:6])))

    np.testing.assert_allclose(dtNew, [1.43, 1.65, 1.93, 3.31], atol=TOL)


def test_mpox_windows(mpox):
    # Observed U.S. mpox counts. Rows 13-32: 11725 cases before the window; the
    # count goes from 14752 to 29751 and doubles once, after 15.09 weeks.
    y = mpox[:, 1]

    prior1 = float(np.sum(y[:12]))

    assert prior1 == 11725

    tds, C0data, _, dtNew = getDoublingTimeCurve(y[12:32], 1, 0, prior1)
    _, _, _, dtOld = getDoublingTimeCurve(y[12:32], 1, 0)

    assert C0data == 14752
    np.testing.assert_allclose(dtNew, [15.09], atol=TOL)
    assert abs(tds[-1, 0] - 15.09) < TOL
    assert tds[-1, 1] <= 2 * 14752

    np.testing.assert_allclose(dtOld, [1.14, 3.09], atol=TOL)

    # Rows 8-27: 1033 cases before the window.
    _, _, _, dtNew = getDoublingTimeCurve(y[7:27], 1, 0, float(np.sum(y[:7])))

    np.testing.assert_allclose(dtNew, [1.17, 1.67, 2.76], atol=TOL)


def test_input_handling():
    # One point: nothing to follow.
    tds, C0data, _, doublingtimes = getDoublingTimeCurve([5], 1, 0)

    assert tds.size == 0
    assert C0data is None
    assert doublingtimes.size == 0

    # A row (2-D) vector is accepted.
    _, _, _, dtRow = getDoublingTimeCurve(np.exp(0.2 * np.arange(20.0)).reshape(1, -1), 1, 0)
    _, _, _, dtColumn = getDoublingTimeCurve(np.exp(0.2 * np.arange(20.0)).reshape(-1, 1), 1, 0)

    assert np.array_equal(dtRow, dtColumn)

    # The count before the window cannot be negative or missing.
    for bad in (-1, float("nan")):
        with pytest.raises(GrowthPredictError) as e:
            getDoublingTimeCurve(np.ones(10), 1, 0, bad)
        assert e.value.identifier == "GrowthPredict:getDoublingTimeCurve:priorCumulative"
