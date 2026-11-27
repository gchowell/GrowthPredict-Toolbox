"""Regression tests for selectRtRows (port of tests/test_selectRtRows.m).

MATLAB indices are one-based; the Python indices are zero-based, so a MATLAB
reference index k is k - 1 here. ``ts[idxRows + 1]`` gives the same times.
"""

import numpy as np
import pytest
from scipy.stats import gamma

from growthpredict import GrowthPredictError, get_Rt, selectRtRows

# Mpox serial interval of options_Rt and cholera serial interval, in weeks.
MPOX = (8.5 / 7, (5.0 / 7) ** 2)
CHOLERA = (5 / 7, (8 / 7) ** 2)


def first_true(mask):
    return int(np.flatnonzero(mask)[0])


def test_mpox_default_drops_first_two_weeks():
    # Window labelled 12..31. Coverage reaches 95% at elapsed time 2.58, so the
    # rows for t = 13 and t = 14 are left out.
    t = np.arange(12.0, 32.0)
    I = 100 * np.exp(-0.2 * np.arange(20.0))

    Rs, ts, giCoverage = get_Rt(t, I, 1, *MPOX)

    idxRows, reliable = selectRtRows(giCoverage)

    assert reliable.shape == Rs.shape
    assert first_true(reliable) == 258 - 1  # reference
    assert abs(giCoverage[199] - 0.867348) < 1e-5  # reference
    assert abs(giCoverage[299] - 0.976684) < 1e-5  # reference

    assert idxRows.tolist() == list(range(299, 1900, 100))
    np.testing.assert_allclose(ts[idxRows + 1], np.arange(15, 32), atol=1e-9)

    # The first reliable index agrees with the distribution's 95th percentile.
    b = MPOX[1] / MPOX[0]
    a = MPOX[0] / b

    assert first_true(reliable) + 1 == int(np.ceil(gamma.ppf(0.95, a, scale=b) / 0.01))


def test_dropped_rows_are_the_biased_ones():
    # Incidence falling by 20% per week: exact R = (1 + r*b)^a = 0.7759. The rows
    # left out are 150% and 20% too high; every row kept is within 4%.
    b = MPOX[1] / MPOX[0]
    a = MPOX[0] / b
    exact = (1 - 0.2 * b) ** a

    t = np.arange(12.0, 32.0)
    I = 100 * np.exp(-0.2 * np.arange(20.0))

    Rs, _, giCoverage = get_Rt(t, I, 1, *MPOX)

    idxRows, _ = selectRtRows(giCoverage)

    assert abs(Rs[99] / exact - 2.496) < 0.01  # reference: +149.6%
    assert abs(Rs[199] / exact - 1.204) < 0.01  # reference: +20.4%

    assert np.max(np.abs(Rs[idxRows] / exact - 1)) < 0.04
    np.testing.assert_allclose(Rs[299], 0.802078, rtol=1e-5)  # reference


def test_other_generation_intervals():
    t = np.arange(12.0, 32.0)
    I = 100 * np.exp(0.1 * np.arange(20.0))

    # Cholera interval: 95% coverage at elapsed time 2.99.
    _, ts, giCoverage = get_Rt(t, I, 1, *CHOLERA)
    idxRows, reliable = selectRtRows(giCoverage)

    assert first_true(reliable) == 300 - 1  # reference
    np.testing.assert_allclose(ts[idxRows + 1], np.arange(15, 32), atol=1e-9)

    # Fixed interval of 1: estimates exist, and are reliable, from elapsed time 1.
    _, ts, giCoverage = get_Rt(t, I, 3, 1, 0)
    idxRows, reliable = selectRtRows(giCoverage)

    assert first_true(reliable) == 100 - 1
    np.testing.assert_allclose(ts[idxRows + 1], np.arange(13, 32), atol=1e-9)


def test_threshold_and_step_arguments():
    _, _, giCoverage = get_Rt(np.arange(12.0, 32.0), 100 * np.exp(0.1 * np.arange(20.0)), 1, *MPOX)

    # Defaults are 0.95 and 100.
    idxDefault, reliableDefault = selectRtRows(giCoverage)
    idxRows, reliable = selectRtRows(giCoverage, 0.95, 100)

    assert np.array_equal(idxRows, idxDefault)
    assert np.array_equal(reliable, reliableDefault)

    # A stricter threshold drops one more week.
    idxRows, _ = selectRtRows(giCoverage, 0.99)

    assert idxRows.tolist() == list(range(399, 1900, 100))

    # Two rows per unit of time.
    idxRows, _ = selectRtRows(giCoverage, 0.95, 50)

    assert idxRows.tolist() == list(range(299, 1900, 50))


def test_last_time_is_always_the_last_row():
    # A curve whose length is not a whole number of time units still ends with
    # its last time.
    giCoverage = np.linspace(0, 1, 250)

    idxRows, reliable = selectRtRows(giCoverage)

    assert first_true(reliable) == 238 - 1
    assert idxRows.tolist() == [249]

    giCoverage = np.concatenate([np.linspace(0, 0.99, 100), np.ones(150)])

    idxRows, _ = selectRtRows(giCoverage)

    assert idxRows.tolist() == [99, 199, 249]


def test_window_too_short_keeps_last_row_only():
    # Generation interval of mean 10 and SD 5 on a window of 5 time units: coverage
    # ends at 0.14. The table keeps the last time, flagged by its coverage.
    _, _, giCoverage = get_Rt(np.arange(6.0), np.exp(0.2 * np.arange(6.0)), 1, 10, 25)

    idxRows, reliable = selectRtRows(giCoverage)

    assert abs(giCoverage[-1] - 0.1429) < 1e-3  # reference
    assert not np.any(reliable)
    assert idxRows.tolist() == [499]


def test_input_checks():
    giCoverage = np.linspace(0, 1, 300)

    for args, ident in (
        ((0,), "GrowthPredict:selectRtRows:minCoverage"),
        ((1.2,), "GrowthPredict:selectRtRows:minCoverage"),
        ((0.95, 0), "GrowthPredict:selectRtRows:stepsPerUnit"),
        ((0.95, 2.5), "GrowthPredict:selectRtRows:stepsPerUnit"),
    ):
        with pytest.raises(GrowthPredictError) as e:
            selectRtRows(giCoverage, *args)
        assert e.value.identifier == ident

    # No estimates, no rows.
    idxRows, reliable = selectRtRows([])

    assert idxRows.size == 0 and reliable.size == 0

    # A column (2-D) input is accepted.
    idxRows, _ = selectRtRows(np.linspace(0, 1, 300).reshape(-1, 1))

    assert idxRows.tolist() == [299]
