"""Regression tests for the titles of the parameter histograms (port of tests/test_parameterTitle.m)."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from growthpredict import GrowthPredictError, fitParameterTitles, parameterTitleText, setParameterTitle

# ---------------------------------------------------------------- text


def test_three_significant_digits():
    inputs = [1.7302, 0.8526, 1, 28547.1, 801.8, 0.17743, 0.0012345, 999.6, 9.996, 123456.7, 10, 0.5, -0.52]
    expected = ["1.73", "0.853", "1", "28500", "802", "0.177", "0.00123", "1000", "10", "123000", "10", "0.5", "-0.52"]
    for value, text in zip(inputs, expected):
        lines = parameterTitleText("r", [value, 0, 1])
        assert lines[0] == "r = " + text, "input %g" % value


def test_exponent_outside_range():
    lines = parameterTitleText("K", [1234567, 2.85e8, 3.2e-5])

    assert lines[0] == "K = 1.23e6"
    assert lines[1] == "95% CI (2.85e8, 3.20e-5)"


def test_quick_start_titles():
    assert parameterTitleText("r", [1.7302, 1.49, 2.08]) == ["r = 1.73", "95% CI (1.49, 2.08)"]
    assert parameterTitleText("p", [0.8526, 0.828, 0.874]) == ["p = 0.853", "95% CI (0.828, 0.874)"]
    assert parameterTitleText("a", [1, 1, 1]) == ["a = 1", "95% CI (1, 1)"]
    assert parameterTitleText("K", [28547, 27295, 29915]) == ["K = 28500", "95% CI (27300, 29900)"]


def test_special_values():
    assert parameterTitleText("a", [0, float("nan"), float("inf")]) == ["a = 0", "95% CI (NaN, Inf)"]
    assert parameterTitleText("a", [0, float("-inf"), 0])[1] == "95% CI (-Inf, 0)"


def test_column_input_and_errors():
    assert parameterTitleText("r", np.array([[1.7302], [1.49], [2.08]])) == ["r = 1.73", "95% CI (1.49, 2.08)"]

    for args, ident in (
        (("r", [1, 2]), "GrowthPredict:parameterTitleText:values"),
        (("r", "abc"), "GrowthPredict:parameterTitleText:values"),
        ((1, [1, 2, 3]), "GrowthPredict:parameterTitleText:name"),
    ):
        with pytest.raises(GrowthPredictError) as e:
            parameterTitleText(*args)
        assert e.value.identifier == ident


# ---------------------------------------------------------------- layout

VALUE_SETS = [
    [[1.7302, 1.49, 2.08], [0.8526, 0.828, 0.874], [1, 1, 1], [28547, 27295, 29915]],
    [[801.8, 571, 1210], [0.178, 0.128, 0.223], [1, 1, 1], [18153, 18026, 18470]],
    [[0.00123, 0.00101, 0.00187], [0.999, 0.951, 1], [4.07, 0.0123, 9.87], [2.85e8, 1.02e8, 9.91e8]],
]


def histogram_figure(sizePt):
    """Figure with the top row of the toolbox figures (four histograms above a wide panel)."""
    fig = plt.figure(figsize=(sizePt[0] / 72, sizePt[1] / 72))
    values = [0.9, 1.0, 1.0, 1.1, 1.1, 1.1, 1.2, 1.2, 1.3, 1.4]
    ax = []
    for k, name in enumerate(("r", "p", "a", "K")):
        a = fig.add_subplot(2, 4, k + 1)
        a.hist(values, bins=10)
        a.set_xlabel(name, fontsize=24)
        a._gp_fontsize = 24
        ax.append(a)
    fig.add_subplot(2, 4, (5, 8))
    return fig, ax


def new_titles(ax, values):
    return [setParameterTitle(name, values[k], ax[k]) for k, name in enumerate(("r", "p", "a", "K"))]


def title_edges(fig, titleHandles):
    """Left and right edge of each title, as fractions of the figure width."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width = fig.bbox.width
    left = np.array([t.get_window_extent(renderer).x0 / width for t in titleHandles])
    right = np.array([t.get_window_extent(renderer).x1 / width for t in titleHandles])
    return left, right


@pytest.mark.parametrize("sizePt", [(864, 475), (560, 420), (1200, 700)])
@pytest.mark.parametrize("v", range(len(VALUE_SETS)))
def test_titles_do_not_overlap(sizePt, v):
    # Size of the README figures (12 x 6.6 in), a small figure and a large figure.
    fig, ax = histogram_figure(sizePt)
    titleHandles = new_titles(ax, VALUE_SETS[v])
    fitParameterTitles(fig)

    left, right = title_edges(fig, titleHandles)

    assert np.all(left[1:4] - right[0:3] > 0), "titles overlap"
    assert right[3] < 1, "last title leaves the figure"
    assert left[0] > 0, "first title leaves the figure"


def test_old_titles_overlapped():
    # The check above must be able to see an overlap: one-line bold titles at
    # the font size of the axes overlap at the size of the README figures.
    fig, ax = histogram_figure((864, 475))
    old = ["r=1.7 (95% CI:1.5,2.1)", "p=0.85 (95% CI:0.83,0.87)", "a=1 (95% CI:1,1)", "K=2.9e+04 (95% CI:2.7e+04,3e+04)"]
    titleHandles = [ax[k].set_title(old[k], fontweight="bold", fontsize=24) for k in range(4)]

    left, right = title_edges(fig, titleHandles)

    assert left[1] - right[0] < 0
    assert left[3] - right[2] < 0


def test_one_font_size_within_limits():
    fig, ax = histogram_figure((864, 475))
    titleHandles = new_titles(ax, VALUE_SETS[0])

    fontSize = fitParameterTitles(fig)

    sizes = [t.get_fontsize() for t in titleHandles]
    np.testing.assert_allclose(sizes, [fontSize] * 4, atol=1e-9)  # the four titles share one font size
    assert fontSize <= 0.75 * 24 + 1e-9
    assert fontSize >= 6

    # Rule by character count for this figure: longest line
    # '95% CI (0.828, 0.874)' and '95% CI (27300, 29900)', 21 characters.
    axWidth = ax[0].get_position().width * 864
    expected = min(0.75 * 24, 1.25 * axWidth / (0.6 * 21), 0.031 * 475)
    assert fontSize <= expected + 0.05
    assert fontSize > 0.6 * expected, "the size is close to the rule, not at the floor"

    assert all(t.get_fontweight() == "normal" for t in titleHandles)


def test_titles_follow_figure_size():
    fig, ax = histogram_figure((1440, 800))
    titleHandles = new_titles(ax, VALUE_SETS[1])

    large = fitParameterTitles(fig)

    fig.set_size_inches(480 / 72, 400 / 72)

    small = fitParameterTitles(fig)

    assert small < large

    left, right = title_edges(fig, titleHandles)
    assert np.all(left[1:4] - right[0:3] > 0)


def test_figure_without_parameter_titles():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title("another title")

    assert fitParameterTitles(fig) is None
