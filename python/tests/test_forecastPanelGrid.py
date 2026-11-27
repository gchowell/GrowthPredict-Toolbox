"""Regression tests for forecastPanelGrid (port of tests/test_forecastPanelGrid.m).

forecastPanelGrid gives the rows and columns of figure 400 of
plotForecast_GrowthModels, which has one forecast panel per rolling window.
"""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from growthpredict import GrowthPredictError, forecastPanelGrid


def factor(n):
    """Prime factors in ascending order, like MATLAB's factor (factor(1) = 1)."""
    if n < 2:
        return [n]
    out, f = [], 2
    while f * f <= n:
        while n % f == 0:
            out.append(f)
            n //= f
        f += 1
    if n > 1:
        out.append(n)
    return out


def test_every_window_has_a_panel():
    # rows*cols equals the number of windows, so a panel exists for every window
    # and no panel is left empty.
    for n in range(1, 201):
        rows, cols = forecastPanelGrid(n)
        assert rows * cols == n, "%d windows" % n
        assert isinstance(rows, int) and isinstance(cols, int) and min(rows, cols) >= 1, "%d windows" % n


def test_counts_that_used_to_fail():
    # Four or more prime factors: the former grid was f1-by-f2 = 2-by-2 for 16, 24
    # and 32 windows and 2-by-3 for 54, fewer panels than windows.
    counts = [16, 24, 32, 36, 40, 48, 54, 56, 60, 64, 72, 80, 81, 96, 100]
    expected = [(4, 4), (4, 6), (4, 8), (6, 6), (5, 8), (6, 8), (6, 9), (7, 8), (6, 10), (8, 8), (8, 9), (8, 10),
                (9, 9), (8, 12), (10, 10)]  # fmt: skip

    for n, layout in zip(counts, expected):
        assert forecastPanelGrid(n) == layout, "%d windows" % n
        f = factor(n)
        assert f[0] * f[1] < n, "the former grid was too small for this count"


def test_other_counts_are_unchanged():
    # One, two or three prime factors: the layout plotForecast_GrowthModels.m used
    # before. The README example with six windows stays 2-by-3.
    for n in range(1, 201):
        f = factor(n)
        if len(f) == 1:
            expected = (f[0], 1)
        elif len(f) == 2:
            expected = (f[0], f[1])
        elif len(f) == 3:
            expected = (f[0] * f[1], f[2])
        else:
            continue
        assert forecastPanelGrid(n) == expected, "%d windows" % n

    assert forecastPanelGrid(6) == (2, 3)
    assert forecastPanelGrid(1) == (1, 1)
    assert forecastPanelGrid(7) == (7, 1)
    assert forecastPanelGrid(12) == (4, 3)


def test_new_layouts_are_the_most_square():
    # For four or more prime factors rows is the largest divisor not above the
    # square root, so rows <= cols and no other divisor pair is closer to a square.
    for n in range(1, 401):
        if len(factor(n)) < 4:
            continue
        rows, cols = forecastPanelGrid(n)
        assert rows <= cols, "%d windows" % n
        assert rows == max(d for d in range(1, n + 1) if n % d == 0 and d * d <= n), "%d windows" % n


def test_subplot_accepts_every_panel():
    # The calls plotForecast_GrowthModels makes: one panel per window, then the
    # bottom row and the first panel again for the axis labels.
    for n in (6, 12, 16, 24, 32):
        rows, cols = forecastPanelGrid(n)
        fig = plt.figure()
        axes = [fig.add_subplot(rows, cols, k) for k in range(1, n + 1)]
        for c in range(1, cols + 1):
            assert 1 <= (rows - 1) * cols + c <= n
        assert len(fig.get_axes()) == len(axes) == n, "%d windows" % n
        plt.close(fig)


def test_input_checks():
    for bad in (0, -3, 2.5, [4, 6], float("nan"), "6", None, True):
        with pytest.raises(GrowthPredictError) as e:
            forecastPanelGrid(bad)
        assert e.value.identifier == "GrowthPredict:forecastPanelGrid:nWindows"

    # A NumPy integer or a whole float is accepted.
    assert forecastPanelGrid(np.int32(16)) == (4, 4)
    assert forecastPanelGrid(16.0) == (4, 4)
