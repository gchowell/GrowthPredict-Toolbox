"""Regression tests for adjustLastWindowStart (port of tests/test_adjustLastWindowStart.m).

adjustLastWindowStart gives the plotting functions the last window start that
the runners use when the requested one does not fit in the series.
"""

import pytest

from growthpredict import GrowthPredictError, adjustLastWindowStart


def runner_rule(tend1, windowsize1, nRows):
    """The lines of Run_Fit_GrowthModels.m and Run_Forecasting_GrowthModels.m."""
    adjusted = False
    if (tend1 + windowsize1 - 1) > nRows:
        tend1 = nRows - windowsize1 + 1
        adjusted = True
    return tend1, adjusted


def test_window_that_fits_is_kept(mpox):
    # U.S. mpox series, 32 rows: the README examples are not adjusted.
    n = mpox.shape[0]

    assert n == 32
    assert adjustLastWindowStart(1, 20, n) == (1, False)  # Run_Fit_GrowthModels(1,1,20)
    assert adjustLastWindowStart(6, 20, n) == (6, False)  # six rolling windows
    assert adjustLastWindowStart(13, 20, n) == (13, False)  # rows 13-32, the last window that fits


def test_window_past_the_end_is_moved_back(mpox):
    # Run_Fit_GrowthModels(1,20,20) saves its files with tend-13; the plotting
    # functions looked for tend-20.
    n = mpox.shape[0]

    assert adjustLastWindowStart(14, 20, n) == (13, True)
    assert adjustLastWindowStart(20, 20, n) == (13, True)

    # options_fit_mpox as shipped asks for row 37 of 32: the last window that
    # fits starts at row 13, below the first window start of 37.
    assert adjustLastWindowStart(37, 20, n) == (13, True)

    # options_fit_cholera: 80 windows of 10 in a series of 107 rows fit.
    assert adjustLastWindowStart(80, 10, 107) == (80, False)


def test_same_as_the_runner_rule():
    # Every combination of a series of 5 to 60 rows, a window of 2 to 40 and a
    # requested start of 1 to 70.
    for n in (5, 16, 32, 50, 60):
        for w in (2, 7, 10, 20, 30, 40):
            for tend1 in range(1, 71):
                actual, adjusted = adjustLastWindowStart(tend1, w, n)
                assert (actual, adjusted) == runner_rule(tend1, w, n), (tend1, w, n)
                if adjusted:
                    assert actual + w - 1 == n  # the adjusted window ends exactly at the last row
                else:
                    assert actual + w - 1 <= n


def test_window_longer_than_the_series():
    # The result is below 1: no window fits, and the caller has nothing to process.
    assert adjustLastWindowStart(1, 40, 32) == (-7, True)


def test_input_checks():
    for args in (([1, 2], 20, 32), (1, None, 32), (1, 20, "n")):
        with pytest.raises(GrowthPredictError) as e:
            adjustLastWindowStart(*args)
        assert e.value.identifier == "GrowthPredict:adjustLastWindowStart:input"
