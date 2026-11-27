"""Titles and markers of the parameter histograms.

Python port of ``parameterTitleText.m``, ``setParameterTitle.m``,
``fitParameterTitles.m``, ``addEstimateLine.m`` and ``GetAdjustedFontSize.m``.
"""

from __future__ import annotations

import math
import re

import numpy as np

from .errors import GrowthPredictError

_TITLE_TAG = "GrowthPredictParameterTitle"

#: Font size of axes labels and tick labels, in points. The MATLAB toolbox uses 24
#: on its own figure windows; 18 gives the same proportions at the figure sizes
#: used here (12 x 6.6 inches for the histogram-and-fit figures).
BASE_FONT_SIZE = 18


def GetAdjustedFontSize():
    """Font size for axes labels and tick labels, in points (``GetAdjustedFontSize.m``).

    MATLAB scales its base size by the screen resolution, capped at 24 points.
    Matplotlib sizes text in points whatever the screen, so the value is the
    module variable :data:`BASE_FONT_SIZE`, capped at 24; change that variable
    (``growthpredict.titles.BASE_FONT_SIZE = 14``) to restyle all figures.
    """
    return min(BASE_FONT_SIZE, 24)


def _formatNumber(x):
    """Three significant digits; no exponent between 1e-3 and 1e6."""
    x = float(x)

    if math.isnan(x):
        return "NaN"

    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"

    if x == 0:
        return "0"

    y = float("%.3g" % x)  # x rounded to three significant digits
    ay = abs(y)

    if 1e-3 <= ay < 1e6:

        decimals = max(0, 2 - int(math.floor(math.log10(ay))))
        s = "%.*f" % (decimals, y)

        if "." in s:
            s = re.sub(r"0+$", "", s)  # trailing zeros
            s = re.sub(r"\.$", "", s)  # and a decimal point left at the end

    else:

        s = "%.2e" % y  # e.g. 2.85e+08
        s = re.sub(r"e([+-])0*(\d)", r"e\1\2", s)
        s = s.replace("e+", "e")  # 2.85e8

    return s


def parameterTitleText(name, values):
    """Two-line title text for a parameter histogram (``parameterTitleText.m``).

    ``name`` is the parameter name, e.g. ``'r'``; ``values`` is
    ``[estimate, lower, upper]``. Returns a list of two strings, e.g.
    ``['r = 1.73', '95% CI (1.49, 2.08)']``.

    Numbers are written with three significant digits. Values with a magnitude
    from 0.001 up to (not including) 1e6 are written without an exponent and
    without trailing zeros (28547 -> '28500', 0.8526 -> '0.853', 1 -> '1');
    others as '2.85e8' or '3.20e-5'.
    """
    if not isinstance(name, str):
        raise GrowthPredictError(
            "GrowthPredict:parameterTitleText:name", "name must be a string, e.g. 'r'."
        )

    ok = not isinstance(values, (str, bytes))
    if ok:
        try:
            vals = np.asarray(values, dtype=float).ravel()
        except (TypeError, ValueError):
            ok = False
    if not ok or vals.size != 3:
        raise GrowthPredictError(
            "GrowthPredict:parameterTitleText:values", "values must be [estimate lower upper]."
        )

    return [
        "%s = %s" % (name, _formatNumber(vals[0])),
        "95%% CI (%s, %s)" % (_formatNumber(vals[1]), _formatNumber(vals[2])),
    ]


def _parameter_titles(fig):
    return [ax.title for ax in fig.get_axes() if ax.title.get_gid() == _TITLE_TAG]


def fitParameterTitles(fig=None, *_):
    """Size the parameter-histogram titles of a figure to their panels.

    Port of ``fitParameterTitles.m``. Finds the titles created by
    :func:`setParameterTitle` in ``fig`` (default: the current figure) and
    gives them one common font size, in points, at which the longest title is
    no wider than the space of its panel. Returns the font size, or ``None``
    when the figure has no such title.

    Font size, for each title, with W the width of its axes in points, n the
    number of characters of its longest line, F the font size of its axes and
    H the height of the figure in points:

        min( 0.75*F,  1.25*W/(0.6*n),  0.031*H ),   not below 6 points

    The smallest value over the titles is used for all of them. The rendered
    width is then measured and the size reduced if a title is wider than 1.15
    axes widths.
    """
    import matplotlib.pyplot as plt

    if fig is None:
        if not plt.get_fignums():
            return None
        fig = plt.gcf()

    titles = _parameter_titles(fig)

    if not titles:
        return None

    capFraction = 0.75  # of the axes font size
    roomFraction = 1.25  # of the axes width, by character count
    charWidth = 0.6  # of the font size, per character (an upper value for digits)
    heightFraction = 0.031  # of the figure height
    maxFraction = 1.15  # of the axes width, as rendered
    minFontSize = 6

    try:

        figWidth, figHeight = fig.get_size_inches() * 72.0  # points

        sizes = []

        for t in titles:

            ax = t.axes
            axWidth = ax.get_position().width * figWidth

            nchar = max(1, max(len(line) for line in t.get_text().split("\n")))

            axesFont = getattr(ax, "_gp_fontsize", None)
            if axesFont is None:
                axesFont = ax.xaxis.label.get_size()

            sizes.append(
                min(capFraction * axesFont, roomFraction * axWidth / (charWidth * nchar), heightFraction * figHeight)
            )

        fontSize = max(minFontSize, min(sizes))

        for t in titles:
            t.set_fontsize(fontSize)

    except Exception:
        # A figure without a usable size: leave the titles as they are.
        return None

    # Rendered width, as a fraction of the axes width.
    try:

        renderer = fig.canvas.get_renderer()
        widest = 0.0

        for t in titles:
            extent = t.get_window_extent(renderer)
            axBox = t.axes.get_window_extent(renderer)
            if axBox.width > 0:
                widest = max(widest, extent.width / axBox.width)

        if math.isfinite(widest) and widest > maxFraction:

            fontSize = max(minFontSize, fontSize * maxFraction / widest)
            for t in titles:
                t.set_fontsize(fontSize)

    except Exception:
        pass  # keep the size from the character count

    return fontSize


def setParameterTitle(name, values, ax=None):
    """Title of a parameter histogram: estimate and 95% CI on two lines.

    Port of ``setParameterTitle.m``. Writes, for example,

        r = 1.73
        95% CI (1.49, 2.08)

    above the axes, in normal weight, and sizes the text to the panel with
    :func:`fitParameterTitles`. The figure is also connected to
    ``fitParameterTitles`` so that the titles keep fitting when it is resized.
    Call it after the font size of the axes is set. Returns the title object.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        ax = plt.gca()

    h = ax.set_title("\n".join(parameterTitleText(name, values)), fontweight="normal")
    h.set_gid(_TITLE_TAG)

    fig = ax.figure

    fitParameterTitles(fig)

    if not getattr(fig, "_gp_title_callback", None):
        try:
            fig._gp_title_callback = fig.canvas.mpl_connect("resize_event", lambda event: fitParameterTitles(fig))
        except Exception:
            pass  # the titles are sized once, above

    return h


def addEstimateLine(estimate, ax=None):
    """Mark a point estimate on a histogram with a vertical dashed line.

    Port of ``addEstimateLine.m``. Draws a red dashed vertical line at
    ``x = estimate`` spanning the full height of the axes and returns it. The
    y-limits are frozen at their value before the line is drawn. If the
    estimate lies outside the current x-limits, the x-limits are widened to
    include it with a 5% margin. A non-finite estimate draws nothing and
    returns ``None``.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        ax = plt.gca()

    try:
        estimate = float(estimate)
    except (TypeError, ValueError):
        return None

    if not math.isfinite(estimate):
        return None

    yl = ax.get_ylim()
    xl = ax.get_xlim()

    if estimate < xl[0] or estimate > xl[1]:
        lo = min(xl[0], estimate)
        hi = max(xl[1], estimate)
        pad = 0.05 * max(hi - lo, np.spacing(max(abs(lo), abs(hi))))
        ax.set_xlim(lo - pad, hi + pad)
    else:
        ax.set_xlim(xl)

    (h,) = ax.plot([estimate, estimate], yl, "r--", linewidth=2)

    ax.set_ylim(yl)

    return h
