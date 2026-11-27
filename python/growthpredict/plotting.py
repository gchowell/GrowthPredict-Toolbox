"""Figures and per-window tables from the results saved by the runners.

Python port of ``plotFit_GrowthModels.m``, ``plotForecast_GrowthModels.m``,
``plotFit_ReproductionNumber.m``, ``plotForecast_ReproductionNumber.m``,
``plotRt.m``, ``plotGrowthModel.m`` and ``plot_growth_illustration.m``.

Each plotting function reads the MAT files saved by its runner, so run the
runner first and call the plotting function with matching options and window
arguments. The functions draw with Matplotlib and return the same numeric
outputs as in MATLAB; the figures stay open (``matplotlib.pyplot.figure(num)``
gives a figure by its number, ``plt.show()`` displays them in an interactive
session, ``fig.savefig(...)`` writes them to a file).

Figure numbers
    100 + i   parameter histograms, fit and forecast of window i
              (runners, plotForecast_GrowthModels)
    101 + i   the same for a fit-only analysis (plotFit_GrowthModels,
              plotFit_ReproductionNumber)
    200 + i   forecast performance against horizon (plotForecast_GrowthModels)
    300       parameter estimates across rolling windows
    400       one forecast panel per rolling window
    500 + k   incidence and reproduction number of the k-th window
    601-604   plotRt: series and R_t, parameters, calibration performance, AICc
"""

from __future__ import annotations

import math
import os
import warnings

import numpy as np

from . import io as gpio
from .doubling import getDoublingTimeCurve
from .errors import GrowthPredictError, warn
from .io import scalar, vector
from .models import growthIncidenceCurve
from .options import options_fit, options_forecast, resolve_options
from .rt import get_Rt, options_Rt, selectRtRows
from .runners import print_options_summary
from .titles import GetAdjustedFontSize, addEstimateLine, setParameterTitle
from .utils import adjustLastWindowStart, load_timeseries, quantile, smooth

__all__ = [
    "plotFit_GrowthModels",
    "plotForecast_GrowthModels",
    "plotFit_ReproductionNumber",
    "plotForecast_ReproductionNumber",
    "plotRt",
    "plotGrowthModel",
    "plot_growth_illustration",
    "draw_window_figure",
    "forecastPanelGrid",
]

_HIST_COLOR = (0.2422, 0.1504, 0.6603)
_FIT_GRAY = (5 / 7, 5 / 7, 5 / 7)  # gray(8), row 6
_WINDOW_FIGSIZE = (12.0, 6.6)
_RT_FIGSIZE = (9.0, 7.0)
_PERF_FIGSIZE = (8.0, 6.0)
SUMMARY_FIGURE = 300


def _plt():
    import matplotlib.pyplot as plt

    return plt


def _style(ax, fs):
    ax._gp_fontsize = fs
    ax.tick_params(labelsize=fs)
    ax.xaxis.label.set_size(fs)
    ax.yaxis.label.set_size(fs)
    ax.title.set_size(fs)
    ax.xaxis.get_offset_text().set_size(0.7 * fs)
    ax.yaxis.get_offset_text().set_size(0.7 * fs)


def _scaled_axis(ax, which, label):
    """Label an axis, moving a power of ten of large or small values into the label.

    Tick labels such as 28000, 30000 become 2.8, 3.0 with the label 'K (x10^4)',
    so that they fit under a narrow panel.
    """
    from matplotlib.ticker import FuncFormatter

    axis = ax.xaxis if which == "x" else ax.yaxis
    lo, hi = ax.get_xlim() if which == "x" else ax.get_ylim()
    top = max(abs(lo), abs(hi))
    e = int(math.floor(math.log10(top))) if top > 0 and math.isfinite(top) else 0
    if e >= 4 or e <= -3:
        scale = 10.0**e
        axis.set_major_formatter(FuncFormatter(lambda v, pos: "%g" % round(v / scale, 6)))
        label = r"%s ($\times10^{%d}$)" % (label, e)
    else:
        axis.set_major_formatter(FuncFormatter(lambda v, pos: "%g" % round(v, 10)))
    if which == "x":
        ax.set_xlabel(label)
    else:
        ax.set_ylabel(label)


# <============================================================================>
# <=========================== Shared drawing code ============================>
# <============================================================================>


def _prediction_limits(forecast_model12):
    """Median and 95% prediction limits across the noisy trajectories (negative limits -> 0)."""
    lims = np.asarray(quantile(forecast_model12, [0.025, 0.975], axis=1))
    LB1 = np.where(lims[0] >= 0, lims[0], 0.0)
    UB1 = np.where(lims[1] >= 0, lims[1], 0.0)
    median1 = np.median(forecast_model12, axis=1)
    return median1, LB1, UB1, float(np.max(lims[1]))


def _draw_fit_panel(ax, timevect1, timevect2, forecast_model12, fit_model1, timevect_all, data_all, forecastingperiod):
    """Trajectories with noise (cyan), median and 95% limits (red), fits (grey) and data (blue)."""
    from matplotlib.collections import LineCollection

    timevect1 = np.asarray(timevect1, dtype=float).ravel()
    timevect2 = np.asarray(timevect2, dtype=float).ravel()

    segs = np.empty((forecast_model12.shape[1], timevect2.size, 2))
    segs[:, :, 0] = timevect2
    segs[:, :, 1] = forecast_model12.T
    ax.add_collection(LineCollection(segs, colors="c", linewidths=0.5))

    median1, LB1, UB1, top = _prediction_limits(forecast_model12)

    # 95% PI
    ax.plot(timevect2, median1, "r-", linewidth=2)
    ax.plot(timevect2, LB1, "r--", linewidth=2)
    ax.plot(timevect2, UB1, "r--", linewidth=2)

    # model fits
    segs = np.empty((fit_model1.shape[1], timevect1.size, 2))
    segs[:, :, 0] = timevect1
    segs[:, :, 1] = fit_model1.T
    ax.add_collection(LineCollection(segs, colors=[_FIT_GRAY], linewidths=1))

    ax.plot(timevect2, median1, "r-", linewidth=2)

    # the data
    ax.plot(timevect_all, data_all, "bo", markerfacecolor="none", markeredgewidth=2)

    ymax = top * 1.5
    if not np.isfinite(ymax) or ymax <= 0:
        ymax = 1.0

    if forecastingperiod > 0:
        ax.plot([timevect1[-1], timevect1[-1]], [0, ymax], "k--", linewidth=2)

    ax.set_xlim(timevect1[0], timevect2[-1])
    ax.set_ylim(0, ymax)

    from matplotlib.ticker import MaxNLocator

    ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))

    return median1, LB1, UB1


def _histogram(ax, values):
    """Ten-bin histogram of the bootstrap values of one parameter, like MATLAB's ``hist``.

    Values that are equal, or equal up to rounding error (a parameter that every
    refit leaves on the same bound, for example), have no range to divide into
    ten bins; they are drawn as one bar around their common value.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    limits = None
    if v.size:
        lo, hi = float(v.min()), float(v.max())
        if hi - lo <= 1e-10 * max(1.0, abs(lo), abs(hi)):
            centre = 0.5 * (lo + hi)
            half = max(0.5, 1e-6 * abs(centre))
            limits = (centre - half, centre + half)
    ax.hist(v, bins=10, range=limits, color=_HIST_COLOR, edgecolor="k", linewidth=0.5)


def draw_window_figure(
    num, *, Phatss_model1, params, timevect1, timevect2, forecast_model12, fit_model1, timevect_all, data_all,
    forecastingperiod, ylabel, model_name1, all_ylabels=True,
):  # fmt: skip
    """Bootstrap histograms of r, p, a and K above the fitted and forecast incidence.

    ``params`` is ``(param_r, param_p, param_a, param_K)``, each
    ``[estimate, lower, upper]``. The red dashed vertical line marks the fit to
    the observed data and the red dashed horizontal segment its 95% interval.
    Returns ``(figure, median1, LB1, UB1)``.
    """
    plt = _plt()
    fs = GetAdjustedFontSize()

    fig = plt.figure(num=num, figsize=_WINDOW_FIGSIZE, clear=True, facecolor="white")
    gs = fig.add_gridspec(2, 4, left=0.095, right=0.985, top=0.895, bottom=0.115, wspace=0.58, hspace=0.62)

    from matplotlib.ticker import MaxNLocator

    for k, (name, param) in enumerate(zip(("r", "p", "a", "K"), params)):

        ax = fig.add_subplot(gs[0, k])
        values = Phatss_model1[:, k]
        _histogram(ax, values)

        ax.plot([param[1], param[2]], [10, 10], "r--", linewidth=2)
        addEstimateLine(param[0], ax)  # fitted estimate

        if all_ylabels or k == 0:
            ax.set_ylabel("Frequency")

        ax.xaxis.set_major_locator(MaxNLocator(nbins=3, min_n_ticks=2))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
        _scaled_axis(ax, "x", name)
        _style(ax, fs)
        setParameterTitle(name, param[:3], ax)  # two lines, sized to the panel

    # <========================================================================================>
    # <================================ Plot model fit and forecast ===========================>
    # <========================================================================================>

    ax = fig.add_subplot(gs[1, :])

    median1, LB1, UB1 = _draw_fit_panel(
        ax, timevect1, timevect2, forecast_model12, fit_model1, timevect_all, data_all, forecastingperiod
    )

    ax.set_xlabel("Time")
    ax.set_ylabel(ylabel)
    _style(ax, fs)
    ax.set_title(model_name1, fontweight="bold", fontsize=fs)

    return fig, median1, LB1, UB1


def _doubling_summary(fit_model1, DT, priorCumulative1, M):
    """Sequence of doubling times across the bootstrap trajectories.

    Returns ``(seq_doublingtimes, param_doubling)``: rows
    ``[i-th doubling, mean, 95% CI LB, 95% CI UB, share of trajectories reaching it]``
    and ``[mean, 2.5%, 97.5%]`` of the mean doubling time per trajectory
    (trajectories without a doubling are left out).
    """
    doublingtimess = np.full((30, M), np.nan)

    maxd = 1

    meandoublingtime = np.full(M, np.nan)

    for j in range(M):

        _, _, _, doublingtimes = getDoublingTimeCurve(np.maximum(fit_model1[:, j], 0), DT, 0, priorCumulative1)

        doublingtimess[: doublingtimes.size, j] = doublingtimes

        if maxd < doublingtimes.size:
            maxd = doublingtimes.size

        if doublingtimes.size:
            meandoublingtime[j] = np.mean(doublingtimes)

    doublingtimess = doublingtimess[:maxd, :M]

    seq_doublingtimes = np.zeros((maxd, 5))

    for j in range(maxd):

        row = doublingtimess[j, :]
        row = row[~np.isnan(row)]

        if row.size:
            seq_doublingtimes[j, :] = [j + 1, np.mean(row), quantile(row, 0.025), quantile(row, 0.975), row.size / M]
        else:
            seq_doublingtimes[j, :] = [j + 1, np.nan, np.nan, np.nan, 0.0]

    # Mean doubling times
    ok = ~np.isnan(meandoublingtime)
    dmean = float(np.mean(meandoublingtime[ok])) if np.any(ok) else float("nan")
    dLB = quantile(meandoublingtime, 0.025)
    dUB = quantile(meandoublingtime, 0.975)

    return seq_doublingtimes, np.array([dmean, dLB, dUB])


def _write_doubling(names, i, seq_doublingtimes):
    gpio.writetable(
        names.window_csv_tstart_first("doublingtimes", i),
        seq_doublingtimes,
        ["i_th doubling", "db mean", "db 95%CI LB", "db 95% CI UB", "prob. i_th doubling"],
    )


def forecastPanelGrid(nWindows):
    """Rows and columns of the panel with one forecast per window (``forecastPanelGrid.m``).

    Returns ``(rows, cols)`` with ``rows*cols == nWindows``, so that every
    window has a panel in figure 400 and no panel is left empty. The layout
    comes from the prime factors of ``nWindows``, f1 <= f2 <= ...:

        one factor (1 or a prime)   nWindows-by-1
        two factors                 f1-by-f2              (6 windows: 2-by-3)
        three factors               (f1*f2)-by-f3         (12 windows: 4-by-3)
        four or more factors        d-by-(nWindows/d), d the largest divisor of
                                    nWindows not above its square root
                                    (16: 4-by-4, 24: 4-by-6, 32: 4-by-8)
    """
    ok = not isinstance(nWindows, (bool, str, bytes)) and np.ndim(nWindows) == 0
    if ok:
        try:
            ok = math.isfinite(nWindows) and nWindows >= 1 and nWindows == round(nWindows)
        except TypeError:
            ok = False
    if not ok:
        raise GrowthPredictError("GrowthPredict:forecastPanelGrid:nWindows", "nWindows must be a positive integer.")

    nWindows = int(nWindows)

    factors = []
    m, f = nWindows, 2
    while m > 1:
        while m % f == 0:
            factors.append(f)
            m //= f
        f += 1

    if len(factors) <= 1:  # 1 or a prime

        rows, cols = nWindows, 1

    elif len(factors) == 2:

        rows, cols = factors[0], factors[1]

    elif len(factors) == 3:

        rows, cols = factors[0] * factors[1], factors[2]

    else:

        # Four or more prime factors: the pair of divisors closest to a square.
        rows = 1

        for d in range(2, nWindows + 1):

            if d * d > nWindows:
                break

            if nWindows % d == 0:
                rows = d

        cols = nWindows // rows

    return rows, cols


def _summary_figure(data, tstart1, tend1, blocks, smoothed):
    """Observed series above the estimates of r, a, p and K for each window start."""
    plt = _plt()
    fs = GetAdjustedFontSize()

    fig = plt.figure(num=SUMMARY_FIGURE, figsize=_WINDOW_FIGSIZE, clear=True, facecolor="white")
    gs = fig.add_gridspec(2, 4, left=0.095, right=0.985, top=0.96, bottom=0.115, wspace=0.62, hspace=0.45)

    ax = fig.add_subplot(gs[0, :])
    ax.plot(data[:, 0], data[:, 1], "ro-", markerfacecolor="none")
    ax.set_xlabel("Time")
    ax.set_ylabel("Cases")
    _style(ax, fs)

    from matplotlib.ticker import MaxNLocator

    x = np.arange(tstart1, tend1 + 1)

    for k, (label, block) in enumerate(blocks):

        ax = fig.add_subplot(gs[1, k])
        ax.plot(x, block[:, 0], "ro-", markerfacecolor="none")
        ax.plot(x, block[:, 1], "b--")
        ax.plot(x, block[:, 2], "b--")

        if smoothed:
            ax.plot(x, smooth(block[:, 0], 5), "k--", linewidth=3)

        ax.set_xlabel("Time")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        _scaled_axis(ax, "y", label)
        _style(ax, fs)

    return fig


# <============================================================================>
# <============================ Shared set-up code ============================>
# <============================================================================>


class _Setup:
    pass


def _setup(tstart1_pass, tend1_pass, windowsize1_pass, forecastingperiod_pass, options_pass, forecast_mode,
           input_dir, output_dir):  # fmt: skip
    """Options, data, windows and file names of a plotting call."""
    s = _Setup()
    s.opts = opts = resolve_options(options_pass, options_forecast if forecast_mode else options_fit)

    # Results are written to ./output; create it when missing.
    os.makedirs(output_dir, exist_ok=True)

    s.data, s.cadfilename1 = load_timeseries(opts.cadfilename1, input_dir)

    s.tstart1 = int(opts.tstart1 if tstart1_pass is None else tstart1_pass)
    s.tend1 = int(opts.tend1 if tend1_pass is None else tend1_pass)
    s.windowsize1 = int(opts.windowsize1 if windowsize1_pass is None else windowsize1_pass)

    if forecast_mode:
        s.getperformance = bool(opts.getperformance)
        s.forecastingperiod = int(opts.forecastingperiod if forecastingperiod_pass is None else forecastingperiod_pass)
    else:
        s.getperformance = True
        s.forecastingperiod = 0

    # The runners move the last window start back when the series is too short for it and
    # name the files they save with the adjusted value (see adjustLastWindowStart).
    s.tend1, tend1Adjusted = adjustLastWindowStart(s.tend1, s.windowsize1, s.data.shape[0])

    s.nothing_to_plot = s.tstart1 > s.tend1
    if s.nothing_to_plot:
        warn(
            "GrowthPredict:windows",
            "incompatible tstart and tend values: the first window starts at row %d, the last possible "
            "window of %d observations starts at row %d. Nothing to plot." % (s.tstart1, s.windowsize1, s.tend1),
            stacklevel=5,
        )

    s.names = gpio.OutputNames(
        output_dir, s.cadfilename1, opts.flag1, opts.fixI0, opts.method1, opts.dist1, s.windowsize1,
        s.forecastingperiod, opts.caddisease, opts.datatype,
    )  # fmt: skip
    s.windows = list(range(s.tstart1, s.tend1 + 1))
    s.ylabel = "%s %s" % (opts.caddisease, opts.datatype)
    return s


class _Window:
    """Results of one window, read from its MAT file."""

    def __init__(self, W):
        self.raw = W
        self.Phatss_model1 = np.asarray(W.Phatss_model1, dtype=float)
        self.forecast_model12 = np.asarray(W.forecast_model12, dtype=float)
        self.forecast_model1 = np.asarray(W.forecast_model1, dtype=float)
        self.fit_model1 = np.asarray(W.fit_model1, dtype=float)
        self.timevect1 = vector(W.timevect1)
        self.timevect2 = vector(W.timevect2)
        self.timevect_all = vector(W.timevect_all)
        self.data_all = vector(W.data_all)
        self.M = int(scalar(W.M))
        self.DT = scalar(W.DT)
        self.params = tuple(vector(W[k]) for k in ("param_r", "param_p", "param_a", "param_K"))
        self.quantilesc = np.asarray(W.quantilesc, dtype=float)
        self.quantilesf = np.asarray(W.quantilesf, dtype=float)

    def table(self, name, ncols=None):
        a = np.asarray(self.raw[name], dtype=float)
        if ncols is not None:
            a = a.reshape(-1, ncols)
        return a


def _load_window(s, i, required=True):
    path = s.names.window_mat(i, s.tend1)
    if not os.path.isfile(path):
        message = (
            'File "%s" does not exist. Run the matching runner first, with the same options, '
            "window range and horizon." % path
        )
        if required:
            raise FileNotFoundError(message)
        print("Error: " + message)
        return None
    return _Window(gpio.load_mat(path))


def _rolling_summary(W):
    """The rolling summaries stored with a window file, as a dict of 2-D arrays."""
    out = {}
    for key in ("param_rs", "param_ps", "param_as", "param_Ks", "param_I0s", "param_alphas", "param_ds",
                "MCSES", "SCIs"):  # fmt: skip
        a = np.asarray(W.raw[key], dtype=float)
        out[key] = a.reshape(-1, 7 if key in ("MCSES", "SCIs") else 3)
    for key in ("AICcs",):
        out[key] = np.asarray(W.raw[key], dtype=float).reshape(-1, 5)
    for key in ("MAECSS", "MSECSS", "PICSS", "WISCSS", "RMSECSS", "MISCSS",
                "MAEFSS", "MSEFSS", "PIFSS", "WISFSS", "RMSEFSS", "MISFSS"):  # fmt: skip
        out[key] = np.asarray(W.raw[key], dtype=float).reshape(-1, 1)
    return out


def _forecast_table(W, windowsize1, forecastingperiod):
    """[time, data, median, LB, UB] over calibration and forecast rows; NaN beyond the data."""
    median1, LB1, UB1, _ = _prediction_limits(W.forecast_model12)
    n = windowsize1 + forecastingperiod
    time = np.full(n, np.nan)
    obs = np.full(n, np.nan)
    m = min(n, W.data_all.size)
    time[:m] = W.timevect_all[:m]
    obs[:m] = W.data_all[:m]
    return np.column_stack([time, obs, median1, LB1, UB1])


def _performance_figure(num, W, fs):
    """MAE, MSE, coverage of the 95% PI and WIS against forecast horizon."""
    plt = _plt()
    panels = (("MAEFS_model1", "MAE"), ("MSEFS_model1", "MSE"),
              ("PIFS_model1", "Coverage rate of the 95% PI"), ("WISFS_model1", "WIS"))  # fmt: skip
    if W.table("MAEFS_model1", 2).shape[0] == 0:
        return None
    fig = plt.figure(num=num, figsize=_PERF_FIGSIZE, clear=True, facecolor="white")
    from matplotlib.ticker import MaxNLocator

    for k, (key, label) in enumerate(panels):
        ax = fig.add_subplot(2, 2, k + 1)
        t = W.table(key, 2)
        ax.plot(t[:, 0], t[:, 1], "k-o", linewidth=4, markerfacecolor="none")
        ax.set_xlabel("Forecasting horizon")
        ax.set_ylabel(label)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        _scaled_axis(ax, "y", label)
        _style(ax, 0.8 * fs)
    fig.tight_layout()
    return fig


def _rt_block(s, W, curves, cc1, optionsRt, fs):
    """Reproduction number of one window: figure 500+cc1, Rt CSV and the latest estimate."""
    plt = _plt()
    type_GId1, mean_GI1, var_GI1 = optionsRt

    # Estimate effective reproduction number, R_t
    Rss = []
    ts = None
    giCoverage = None
    for realization1 in range(W.M):
        Rs, ts, giCoverage = get_Rt(W.timevect2, curves[:, realization1], type_GId1, mean_GI1, var_GI1)
        Rss.append(Rs)
    Rss = np.column_stack(Rss)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        q = np.asarray(quantile(Rss, [0.5, 0.025, 0.975], axis=1))
    Rtmedian, RtLB, RtUB = q[0], q[1], q[2]

    # Estimates near the start of the window are too high, because the incidence
    # before the window is not used (see get_Rt). Show and export those for which
    # the history covers at least minCoverageRt of the generation-interval distribution.
    minCoverageRt = 0.95

    idxRt, reliableRt = selectRtRows(giCoverage, minCoverageRt, 100)

    tsRt = ts[1:]

    shownRt = reliableRt.copy()

    if not np.any(shownRt):
        warn(
            "GrowthPredict:Rt:coverage",
            "The window is too short for this generation interval: the incidence history never covers %g%% of "
            "its distribution. The R_t estimates, including the one at the last time, are too high."
            % (100 * minCoverageRt),
        )
        shownRt[:] = True

    shownUB = RtUB[shownRt]
    shownUB = shownUB[np.isfinite(shownUB)]
    RtYmax = float(np.max(shownUB)) if shownUB.size else 2.0
    RtYmax = max(1.5, 1.2 * RtYmax)

    fig = plt.figure(num=500 + cc1, figsize=_RT_FIGSIZE, clear=True, facecolor="white")

    ax = fig.add_subplot(2, 1, 1)
    _draw_fit_panel(
        ax, W.timevect1, W.timevect2, W.forecast_model12, W.fit_model1, W.timevect_all, W.data_all,
        s.forecastingperiod,
    )  # fmt: skip
    ax.set_xlabel("Time")
    ax.set_ylabel(s.ylabel)
    _style(ax, fs)

    ax = fig.add_subplot(2, 1, 2)

    ax.plot(tsRt[shownRt], Rtmedian[shownRt], "r-", linewidth=2)
    ax.plot(tsRt[shownRt], RtLB[shownRt], "r--", linewidth=2)
    ax.plot(tsRt[shownRt], RtUB[shownRt], "r--", linewidth=2)

    ax.plot([ts[0], ts[-1]], [1, 1], "k--", linewidth=2)

    ax.set_xlim(W.timevect1[0], W.timevect2[-1])
    ax.set_ylim(0, RtYmax)

    ax.plot([W.timevect1[-1], W.timevect1[-1]], [0, RtYmax], "k--", linewidth=2)

    ax.set_xlabel("Time")
    ax.set_ylabel(r"$R_t$")
    _style(ax, fs)
    fig.tight_layout()

    # <=============================================================================================>
    # <============================== Save file with R_t estimates =================================>
    # <=============================================================================================>

    # Rtmedian[k] is the estimate at ts[k+1] (see get_Rt). idxRt (selectRtRows) holds one row per
    # unit of time after the start of the window, ending at the last time, without the rows whose
    # history covers less than minCoverageRt of the generation-interval distribution. The last
    # column gives that coverage.
    Rtdata = np.column_stack([ts[idxRt + 1], Rtmedian[idxRt], RtLB[idxRt], RtUB[idxRt], giCoverage[idxRt]])

    gpio.writetable(
        s.names.window_csv_tstart_first("Rt", W.i),
        Rtdata,
        ["time", "Rt median", "Rt 95%CI LB", "Rt 95% CI UB", "GI coverage"],
    )

    # most recent R estimate
    last = Rss[-1, :]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        param_R = [
            ts[-1], float(np.median(last)), quantile(last, 0.025), quantile(last, 0.975),
            float(np.std(last, ddof=1)) if last.size > 1 else 0.0,
        ]  # fmt: skip
    return param_R


# <============================================================================>
# <============================== Fit-only plots ==============================>
# <============================================================================>


def _plot_fit(tstart1_pass, tend1_pass, windowsize1_pass, options_pass, optionsRt_pass, with_rt,
              input_dir, output_dir, verbose):  # fmt: skip

    plt = _plt()
    plt.close("all")

    s = _setup(tstart1_pass, tend1_pass, windowsize1_pass, None, options_pass, False, input_dir, output_dir)
    if s.nothing_to_plot:
        return np.zeros((0, 5)), np.zeros(0), np.zeros((0, 5)), np.zeros((0, 0))
    opts = s.opts
    fs = GetAdjustedFontSize()

    optionsRt = None
    if with_rt:
        optionsRt = (options_Rt if optionsRt_pass is None else optionsRt_pass)()

    param_R = []
    param_doubling = None
    seq_doublingtimes = None
    W = None
    cc1 = 1

    for i in s.windows:  # rolling window analysis

        if (s.tend1 - s.tstart1) > 10:
            plt.close("all")

        Wi = _load_window(s, i, required=False)
        if Wi is None:
            continue
        W = Wi
        W.i = i
        if verbose:
            print("File loaded successfully.")

        # <======================================================================================>
        # <======================= Plot parameter distributions and model fit ===================>
        # <======================================================================================>

        _, median1, LB1, UB1 = draw_window_figure(
            101 + i,
            Phatss_model1=W.Phatss_model1,
            params=W.params,
            timevect1=W.timevect1,
            timevect2=W.timevect2,
            forecast_model12=W.forecast_model12,
            fit_model1=W.fit_model1,
            timevect_all=W.timevect_all,
            data_all=W.data_all,
            forecastingperiod=s.forecastingperiod,
            ylabel=s.ylabel,
            model_name1=opts.model_name1,
            all_ylabels=with_rt,
        )

        # compute doubling times

        # Doubling times refer to the cumulative count of the series, so the cases
        # reported before the window are added (see getDoublingTimeCurve). Without
        # them the count restarts at every window.
        priorCumulative1 = float(np.nansum(np.maximum(s.data[: i - 1, 1], 0)))

        seq_doublingtimes, param_doubling = _doubling_summary(W.fit_model1, W.DT, priorCumulative1, W.M)

        _write_doubling(s.names, i, seq_doublingtimes)

        fitdata = np.column_stack([W.timevect1, W.data_all[: W.timevect1.size], median1, LB1, UB1])

        gpio.writetable(s.names.window_csv_tstart_first("Fit", i), fitdata, ["time", "data", "median", "LB", "UB"])

        if with_rt:

            param_R.append(_rt_block(s, W, W.fit_model1, cc1, optionsRt, fs))

        else:

            # Save quantiles of the model fit
            gpio.write_quantiles(s.names.window_csv_tstart_first("quantile-fit", i), W.quantilesc, W.quantilesf)

        cc1 += 1

    if W is None:
        raise FileNotFoundError(
            "No result file was found for windows %d to %d. Run Run_Fit_GrowthModels first, with the same "
            "options and window arguments." % (s.tstart1, s.tend1)
        )

    summary = _rolling_summary(W)

    # <==================================================================================================>
    # <====================== Plot temporal variation of parameters from rolling window analysis =======>
    # <==================================================================================================>

    if s.tend1 > s.tstart1 and summary["param_rs"].shape[0] == len(s.windows):
        _summary_figure(
            s.data, s.tstart1, s.tend1,
            (("r", summary["param_rs"]), ("a", summary["param_as"]), ("p", summary["param_ps"]),
             ("K", summary["param_Ks"])),
            smoothed=False,
        )  # fmt: skip

    # <=====================================================================================================>
    # <============================== Save the CSV files of the rolling window analysis ====================>
    # <=====================================================================================================>

    if summary["param_rs"].shape[0] == len(s.windows):
        gpio.write_performance_calibration(s.names, s.tstart1, s.tend1, s.windowsize1, summary)
        gpio.write_rolling_tables(s.names, s.tstart1, s.tend1, opts.method1, summary)

    if verbose:
        print_options_summary(
            "options_fit", s.cadfilename1, opts.caddisease, opts.datatype, opts.method1, opts.dist1,
            opts.numstartpoints, W.M, opts.flag1, opts.model_name1, opts.fixI0, s.windowsize1, s.tstart1, s.tend1,
        )  # fmt: skip

    return np.array(param_R, dtype=float).reshape(-1, 5), param_doubling, seq_doublingtimes, W.fit_model1


def plotFit_GrowthModels(tstart1_pass=None, tend1_pass=None, windowsize1_pass=None, options_pass=None, *,
                         input_dir="./input", output_dir="./output", verbose=True):  # fmt: skip
    """Plot the fit of each window with quantified uncertainty and export its tables.

    Port of ``plotFit_GrowthModels.m``. Run :func:`Run_Fit_GrowthModels` first
    with the same options and window arguments.

    Draws figure ``101 + i`` for window ``i`` and writes ``Fit-*.csv``,
    ``quantile-fit-*.csv`` and ``doublingtimes-*.csv`` per window, plus the
    rolling CSV summaries. With more than one window, figure 300 shows the
    parameter estimates across windows.

    Returns
    -------
    param_doubling : ndarray
        ``[mean, 2.5%, 97.5%]`` of the mean doubling time (last window).
    seq_doublingtimes : ndarray
        Rows ``[i-th doubling, mean, 95% CI LB, 95% CI UB, share of
        trajectories reaching it]`` (last window).
    fit_model1 : ndarray
        Bootstrap model trajectories of the last window.
    """
    _, param_doubling, seq_doublingtimes, fit_model1 = _plot_fit(
        tstart1_pass, tend1_pass, windowsize1_pass, options_pass, None, False, input_dir, output_dir, verbose
    )
    return param_doubling, seq_doublingtimes, fit_model1


def plotFit_ReproductionNumber(tstart1_pass=None, tend1_pass=None, windowsize1_pass=None, options_pass=None,
                               optionsRt_pass=None, *, input_dir="./input", output_dir="./output",
                               verbose=True):  # fmt: skip
    """Plot the fit and the effective reproduction number of each window.

    Port of ``plotFit_ReproductionNumber.m``. Run :func:`Run_Fit_GrowthModels`
    first. ``optionsRt_pass`` is a function without inputs that returns
    ``(type_GId1, mean_GI1, var_GI1)``, such as ``options_Rt`` (the default).

    Draws figure ``101 + i`` (fit) and figure ``500 + k`` (incidence above
    R_t) for the k-th window, and writes ``Rt-*.csv``, ``Fit-*.csv`` and
    ``doublingtimes-*.csv``. R_t is shown and exported only where the incidence
    history inside the window covers at least 95% of the generation-interval
    distribution.

    Returns
    -------
    param_R : ndarray, one row per window
        ``[last time, median R, 2.5%, 97.5%, SD]`` of the most recent estimate.
    param_doubling, seq_doublingtimes : ndarray
        As in :func:`plotFit_GrowthModels`.
    """
    param_R, param_doubling, seq_doublingtimes, _ = _plot_fit(
        tstart1_pass, tend1_pass, windowsize1_pass, options_pass, optionsRt_pass, True, input_dir, output_dir, verbose
    )
    return param_R, param_doubling, seq_doublingtimes


# <============================================================================>
# <============================== Forecast plots ==============================>
# <============================================================================>


def _plot_forecast(tstart1_pass, tend1_pass, windowsize1_pass, forecastingperiod_pass, options_pass,
                   optionsRt_pass, with_rt, input_dir, output_dir, verbose):  # fmt: skip

    plt = _plt()
    plt.close("all")

    s = _setup(
        tstart1_pass, tend1_pass, windowsize1_pass, forecastingperiod_pass, options_pass, True, input_dir, output_dir
    )
    if s.nothing_to_plot:
        return np.zeros((0, 5)), np.zeros(0), np.zeros((0, 5)), np.zeros((0, 0))
    opts = s.opts
    fs = GetAdjustedFontSize()
    getperformance = s.getperformance
    forecastingperiod = s.forecastingperiod

    optionsRt = None
    if with_rt:
        optionsRt = (options_Rt if optionsRt_pass is None else optionsRt_pass)()

    rows, cols = forecastPanelGrid(len(s.windows))  # layout of figure 400, one panel per window

    acc = {k: [] for k in ("param_rs", "param_as", "param_ps", "param_Ks", "MAECSS", "MSECSS", "PICSS", "WISCSS",
                           "MAEFSS", "MSEFSS", "PIFSS", "WISFSS")}  # fmt: skip

    param_R = []
    param_doubling = None
    seq_doublingtimes = None
    W = None
    panel_axes = []
    fig400 = None
    cc1 = 1

    for i in s.windows:  # rolling window analysis

        if (s.tend1 - s.tstart1) > 10:
            # many windows: keep only the panel of forecasts (figure 400) open
            for num in plt.get_fignums():
                if num != 400:
                    plt.close(num)

        W = _load_window(s, i, required=True)
        W.i = i

        # <======================================================================================>
        # <======================= Plot parameter distributions and model fit and forecast ======>
        # <======================================================================================>

        _, median1, LB1, UB1 = draw_window_figure(
            100 + (cc1 if with_rt else i),
            Phatss_model1=W.Phatss_model1,
            params=W.params,
            timevect1=W.timevect1,
            timevect2=W.timevect2,
            forecast_model12=W.forecast_model12,
            fit_model1=W.fit_model1,
            timevect_all=W.timevect_all,
            data_all=W.data_all,
            forecastingperiod=forecastingperiod,
            ylabel=s.ylabel,
            model_name1=opts.model_name1,
            all_ylabels=True,
        )

        # store parameter estimates
        for key, param in zip(("param_rs", "param_ps", "param_as", "param_Ks"), W.params):
            acc[key].append(param)

        # compute doubling times (of the cumulative count of the series)
        priorCumulative1 = float(np.nansum(np.maximum(s.data[: i - 1, 1], 0)))

        seq_doublingtimes, param_doubling = _doubling_summary(W.fit_model1, W.DT, priorCumulative1, W.M)

        _write_doubling(s.names, i, seq_doublingtimes)

        if with_rt:

            param_R.append(_rt_block(s, W, W.forecast_model1, cc1, optionsRt, fs))

        elif len(s.windows) > 1:

            if fig400 is None:
                fig400 = plt.figure(num=400, figsize=_WINDOW_FIGSIZE, clear=True, facecolor="white")
            ax = fig400.add_subplot(rows, cols, len(panel_axes) + 1)
            panel_axes.append(ax)

            _draw_fit_panel(
                ax, W.timevect1, W.timevect2, W.forecast_model12, W.fit_model1, W.timevect_all, W.data_all,
                forecastingperiod,
            )  # fmt: skip
            _style(ax, min(16, fs))

        if getperformance and forecastingperiod > 0 and W.data_all.size < s.windowsize1 + forecastingperiod:
            warn(
                "GrowthPredict:seriesTooShort",
                "Length of time series data is too short to evaluate the forecasting period indicated in "
                "<forecastingperiod>. Consider setting <getperformance> to 0 in the forecast options or "
                "extending the length of the time series.",
            )

        gpio.writetable(
            s.names.window_csv_tstart_last("Forecast", i),
            _forecast_table(W, s.windowsize1, forecastingperiod),
            ["time", "data", "median", "LB", "UB"],
        )

        if not with_rt:
            # Save quantiles of the model forecast
            gpio.write_quantiles(s.names.window_csv_tstart_last("quantile-forecast", i), W.quantilesc, W.quantilesf)

        # <========================================================================================>
        #                Plot forecasting performance metrics over the predicted horizon
        # <========================================================================================>

        evaluated = W.table("MAEFS_model1", 2).shape[0] > 0

        if forecastingperiod > 0 and getperformance and evaluated:
            _performance_figure(200 + (cc1 if with_rt else i), W, fs)

        # store metrics for calibration
        acc["MAECSS"].append(W.table("MAECS_model1")[-1, -1])
        acc["MSECSS"].append(W.table("MSECS_model1")[-1, -1])
        acc["PICSS"].append(W.table("PICS_model1")[-1, -1])
        acc["WISCSS"].append(scalar(W.raw["WISC_model1"]))

        # store metrics for short-term forecasts
        if forecastingperiod > 0 and getperformance and evaluated:
            acc["MAEFSS"].append(W.table("MAEFS_model1", 2)[-1, -1])
            acc["MSEFSS"].append(W.table("MSEFS_model1", 2)[-1, -1])
            acc["PIFSS"].append(W.table("PIFS_model1", 2)[-1, -1])
            acc["WISFSS"].append(W.table("WISFS_model1", 2)[-1, -1])

        if not with_rt:
            gpio.write_quantiles(s.names.window_csv_tstart_first("quantile", i), W.quantilesc, W.quantilesf)

        cc1 += 1

    if panel_axes:

        for c in range(cols):
            k = (rows - 1) * cols + c
            if k < len(panel_axes):
                panel_axes[k].set_xlabel("Time")

        panel_axes[0].set_ylabel(s.ylabel)
        fig400.tight_layout()

    def col(key):
        return np.asarray(acc[key], dtype=float).reshape(-1, 1)

    blocks = {k: np.vstack(acc[k]) for k in ("param_rs", "param_as", "param_ps", "param_Ks")}

    # <==================================================================================================>
    # <====================== Plot temporal variation of parameters from rolling window analysis =======>
    # <==================================================================================================>

    if s.tend1 > s.tstart1:
        _summary_figure(
            s.data, s.tstart1, s.tend1,
            (("r", blocks["param_rs"]), ("a", blocks["param_as"]), ("p", blocks["param_ps"]),
             ("K", blocks["param_Ks"])),
            smoothed=True,
        )  # fmt: skip

    # <=====================================================================================================>
    # <============================== Save the CSV files of the rolling window analysis ====================>
    # <=====================================================================================================>

    summary = _rolling_summary(W)

    if summary["param_rs"].shape[0] == len(s.windows):
        gpio.write_rolling_tables(s.names, s.tstart1, s.tend1, opts.method1, summary)
        if not with_rt:
            gpio.write_aicc_table(s.names, s.tstart1, s.tend1, summary["AICcs"])

    time = np.arange(s.tstart1, s.tend1 + 1, dtype=float).reshape(-1, 1)

    gpio.writetable(
        s.names.rolling_csv("performance-calibration", s.tstart1, s.tend1),
        np.hstack([time, np.zeros_like(time) + s.windowsize1, col("MAECSS"), col("MSECSS"), col("PICSS"), col("WISCSS")]),
        ["time", "calibration_period", "MAE", "MSE", "Coverage 95%PI", "WIS"],
    )

    if getperformance and forecastingperiod > 0 and len(acc["MAEFSS"]) == len(s.windows):
        gpio.writetable(
            s.names.rolling_csv("performance-forecasting", s.tstart1, s.tend1),
            np.hstack([time, np.zeros_like(time) + forecastingperiod, col("MAEFSS"), col("MSEFSS"), col("PIFSS"), col("WISFSS")]),
            ["time", "Horizon", "MAE", "MSE", "Coverage 95%PI", "WIS"],
        )

    if verbose:
        print_options_summary(
            "options_forecast", s.cadfilename1, opts.caddisease, opts.datatype, opts.method1, opts.dist1,
            opts.numstartpoints, W.M, opts.flag1, opts.model_name1, opts.fixI0, s.windowsize1, s.tstart1, s.tend1,
            forecasting=(getperformance, forecastingperiod),
        )  # fmt: skip

    return np.array(param_R, dtype=float).reshape(-1, 5), param_doubling, seq_doublingtimes, W.forecast_model1


def plotForecast_GrowthModels(tstart1_pass=None, tend1_pass=None, windowsize1_pass=None,
                              forecastingperiod_pass=None, options_pass=None, *, input_dir="./input",
                              output_dir="./output", verbose=True):  # fmt: skip
    """Plot the fit and forecast of each window and export its tables.

    Port of ``plotForecast_GrowthModels.m``. Run
    :func:`Run_Forecasting_GrowthModels` first with the same options, window
    arguments and horizon.

    Draws figure ``100 + i`` for window ``i`` and, when performance metrics
    are enabled, figure ``200 + i`` (MAE, MSE, coverage and WIS against
    horizon). With more than one window, figure 400 has one forecast panel per
    window and figure 300 the parameter estimates across windows with a
    smoothed trend. Writes ``Forecast-*.csv``, ``quantile-forecast-*.csv``,
    ``quantile-*.csv`` and ``doublingtimes-*.csv`` per window, plus the rolling
    CSV summaries.

    Returns
    -------
    param_doubling, seq_doublingtimes : ndarray
        Doubling times of the last window (see :func:`plotFit_GrowthModels`).
    forecast_model1 : ndarray
        Bootstrap model trajectories of the last window.
    """
    _, param_doubling, seq_doublingtimes, forecast_model1 = _plot_forecast(
        tstart1_pass, tend1_pass, windowsize1_pass, forecastingperiod_pass, options_pass, None, False,
        input_dir, output_dir, verbose,
    )  # fmt: skip
    return param_doubling, seq_doublingtimes, forecast_model1


def plotForecast_ReproductionNumber(tstart1_pass=None, tend1_pass=None, windowsize1_pass=None,
                                    forecastingperiod_pass=None, options_pass=None, optionsRt_pass=None, *,
                                    input_dir="./input", output_dir="./output", verbose=True):  # fmt: skip
    """Plot the forecast and the effective reproduction number of each window.

    Port of ``plotForecast_ReproductionNumber.m``. Run
    :func:`Run_Forecasting_GrowthModels` first. R_t is computed from the
    bootstrap model trajectories over the calibration and forecast periods.

    Draws figure ``100 + k`` (fit and forecast), ``500 + k`` (incidence above
    R_t) and ``200 + k`` (forecast performance) for the k-th window, and
    writes ``Rt-*.csv``, ``Forecast-*.csv`` and ``doublingtimes-*.csv``.

    Returns ``(param_R, param_doubling, seq_doublingtimes)`` as
    :func:`plotFit_ReproductionNumber`.
    """
    param_R, param_doubling, seq_doublingtimes, _ = _plot_forecast(
        tstart1_pass, tend1_pass, windowsize1_pass, forecastingperiod_pass, options_pass, optionsRt_pass, True,
        input_dir, output_dir, verbose,
    )  # fmt: skip
    return param_R, param_doubling, seq_doublingtimes


# <============================================================================>
# <================================== plotRt ==================================>
# <============================================================================>


def plotRt(tstart1_pass=None, tend1_pass=None, windowsize1_pass=None, options_pass=None, optionsRt_pass=None, *,
           input_dir="./input", output_dir="./output", verbose=True):  # fmt: skip
    """Reproduction number, parameters and calibration performance across rolling windows.

    Port of ``plotRt.m``. Reads, for a fit-only rolling analysis, the
    ``Rt-*.csv`` file of every window (written by
    :func:`plotFit_ReproductionNumber`), and the ``parameters-``,
    ``performance-calibration-`` and ``AICcs-rollingwindow`` CSV files. Takes
    the last row of each window's Rt file, i.e. the estimate at the end of the
    window.

    Draws figures 601 (series above R_t), 602 (parameters), 603 (calibration
    performance) and 604 (AICc), writes ``Rts-*.csv`` and returns
    ``(Rts, params, performanceC, AICcs)``.
    """
    plt = _plt()
    plt.close("all")

    s = _setup(tstart1_pass, tend1_pass, windowsize1_pass, None, options_pass, False, input_dir, output_dir)
    if s.nothing_to_plot:
        return np.zeros((0, 4)), np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((0, 5))
    opts = s.opts
    fs = 0.75 * GetAdjustedFontSize()
    tstart1, tend1, windowsize1 = s.tstart1, s.tend1, s.windowsize1

    Rts = np.zeros((len(s.windows), 4))

    for cc1, i in enumerate(s.windows):  # rolling window analysis

        # time, Rt median, Rt 95%CI LB, Rt 95% CI UB
        _, Rt = gpio.readtable(s.names.window_csv_tstart_first("Rt", i))
        Rts[cc1, :] = Rt[-1, 0:4]

    _, params = gpio.readtable(s.names.rolling_csv("parameters-rollingwindow", tstart1, tend1))

    M = opts.B
    if verbose:
        print_options_summary(
            "options_fit", s.cadfilename1, opts.caddisease, opts.datatype, opts.method1, opts.dist1,
            opts.numstartpoints, M, opts.flag1, opts.model_name1, opts.fixI0, windowsize1, tstart1, tend1,
        )  # fmt: skip

    # plot epidemic curve and reproduction number

    fig = plt.figure(num=601, figsize=_RT_FIGSIZE, clear=True, facecolor="white")

    ax = fig.add_subplot(2, 1, 1)
    curve = s.data
    rows = slice(tstart1 - 1, tend1 + windowsize1 - 1)
    ax.plot(curve[rows, 0], curve[rows, 1], "bo-", markerfacecolor="none")
    ax.set_xlabel("Time")
    ax.set_ylabel(s.ylabel)
    _style(ax, fs)

    ax = fig.add_subplot(2, 1, 2)
    ax.plot(Rts[:, 0], Rts[:, 1], "r.-")
    ax.plot(Rts[:, 0], Rts[:, 2:4], "r--")
    ax.plot([tstart1, Rts[-1, 0]], [1, 1], "k--")
    ax.set_xlabel("Time")
    if Rts[-1, 0] > tstart1:
        ax.set_xlim(tstart1, Rts[-1, 0])
    lo, hi = np.nanmin(Rts[:, 2]) - 0.1, np.nanmax(Rts[:, 3]) + 0.1
    if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
        ax.set_ylim(lo, hi)
    ax.set_ylabel(r"Reproduction number, $R_t$")
    _style(ax, fs)
    fig.tight_layout()

    # save the Rts values in a csv file
    gpio.writetable(s.names.rts_csv(tstart1, tend1), Rts, ["time", "mean Rt", "Rt LB", "Rt UB"])

    # plot parameter estimates

    x = np.arange(tstart1, tend1 + 1)

    fig = plt.figure(num=602, figsize=_PERF_FIGSIZE, clear=True, facecolor="white")
    for k, label in enumerate(("r", "p", "a", "K")):
        ax = fig.add_subplot(2, 2, k + 1)
        c = 1 + 3 * k
        ax.plot(x, params[:, c], "ro-", markerfacecolor="none")
        ax.plot(x, params[:, c + 1], "b--")
        ax.plot(x, params[:, c + 2], "b--")
        ax.set_ylabel(label)
        ax.set_xlabel("Time")
        _style(ax, fs)
    fig.tight_layout()

    # plot calibration performance metrics

    _, performanceC = gpio.readtable(s.names.rolling_csv("performance-calibration", tstart1, tend1))

    xe = np.arange(tstart1 + windowsize1 - 1, tend1 + windowsize1)

    fig = plt.figure(num=603, figsize=_PERF_FIGSIZE, clear=True, facecolor="white")
    for k, (label, c) in enumerate((("MAE", 2), ("MSE", 3), ("Coverage 95% PI", 4), ("WIS", 5))):
        ax = fig.add_subplot(2, 2, k + 1)
        ax.plot(xe, performanceC[:, c], "r-")
        ax.set_ylabel(label)
        ax.set_xlabel("Time")
        _style(ax, fs)
        if verbose:
            print("mean %s=%s" % (label[0].lower() + label[1:] if label.startswith("Cov") else label, "%.5g" % np.mean(performanceC[:, c])))
    fig.tight_layout()

    # plot AICc

    _, AICcs = gpio.readtable(s.names.rolling_csv("AICcs-rollingwindow", tstart1, tend1))

    fig = plt.figure(num=604, figsize=(6.4, 4.8), clear=True, facecolor="white")
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(xe, AICcs[:, 1], "r-")
    ax.set_ylabel("AICc")
    ax.set_xlabel("Time")
    _style(ax, fs)
    fig.tight_layout()
    if verbose:
        print("mean AICc=%.5g" % np.mean(AICcs[:, 1]))

    return Rts, params, performanceC, AICcs


# <============================================================================>
# <=============================== Illustrations ==============================>
# <============================================================================>


def plotGrowthModel(flag1_pass=None, r_pass=None, p_pass=None, a_pass=None, K_pass=None, I0_pass=None,
                    windowsize1_pass=None, ax=None):  # fmt: skip
    """Plot the solution of a growth model for the parameter values supplied.

    Port of ``plotGrowthModel.m``. The curve is solved at times
    ``0, 1, ..., windowsize1`` from ``C(0) = I0`` and drawn in the current
    axes (or ``ax``).

    Required values: GGM (0): r, p; GLM (1): r, p, K; GRM (2): r, p, a, K;
    LM (3): r, K; RICH (4): r, a, K; GOM (5): r, a, K. A missing value raises
    :class:`GrowthPredictError` (MATLAB prints a message and returns).

    Returns ``(curve, params)``: the incidence curve ``abs([C(1); diff(C)])``
    and the vector of the model's own parameters.
    """

    def need(value, text):
        if value is None:
            raise GrowthPredictError("GrowthPredict:plotGrowthModel:missing", text)
        return value

    flag1 = need(flag1_pass, "flag1 value is required")

    a1 = 1.0
    p1 = 1.0
    K1 = 1.0

    if flag1 == 0:
        r1 = need(r_pass, "r value is required for the GGM model")
        p1 = need(p_pass, "p value is required for the GGM model")
        model_name1 = "Generalized Growth Model"
        params = [r1, p1]
    elif flag1 == 1:
        r1 = need(r_pass, "r value is required for the GLM model")
        p1 = need(p_pass, "p value is required for the GLM model")
        K1 = need(K_pass, "K value is required for the GLM model")
        model_name1 = "Generalized Logistic Growth Model"
        params = [r1, p1, K1]
    elif flag1 == 2:
        r1 = need(r_pass, "r value is required for the Generalized Richards model")
        p1 = need(p_pass, "p value is required for the Generalized Richards model")
        a1 = need(a_pass, "a value is required for the Generalized Richards model")
        K1 = need(K_pass, "K value is required for the Generalized Richards model")
        model_name1 = "Generalized Richards Model"
        params = [r1, p1, a1, K1]
    elif flag1 == 3:
        r1 = need(r_pass, "r value is required for the LM model")
        K1 = need(K_pass, "K value is required for the LM model")
        model_name1 = "Logistic Growth Model"
        params = [r1, K1]
    elif flag1 == 4:
        r1 = need(r_pass, "r value is required for the Richards model")
        a1 = need(a_pass, "a value is required for the Richards model")
        K1 = need(K_pass, "K value is required for the Richards model")
        model_name1 = "Richards Model"
        params = [r1, a1, K1]
    elif flag1 == 5:  # Gompertz
        r1 = need(r_pass, "r value is required for the Gompertz model")
        a1 = need(a_pass, "a value is required for the Gompertz model")
        K1 = need(K_pass, "K value is required for the Gompertz model")
        model_name1 = "Gompertz Model"
        params = [r1, a1, K1]
    else:
        raise GrowthPredictError(
            "GrowthPredict:plotGrowthModel:model", "flag1 = %g is not a growth model of plotGrowthModel (0 to 5)." % flag1
        )

    I0 = need(I0_pass, "Initial value C(0) is required")

    windowsize1 = need(windowsize1_pass, "Simulation duration is required")

    timevect = np.arange(0, int(windowsize1) + 1, dtype=float)

    curve = growthIncidenceCurve([r1, p1, a1, K1, I0], timevect, flag1)

    plt = _plt()
    if ax is None:
        ax = plt.gca()

    ax.plot(timevect, curve, "b-", linewidth=2)

    ax.set_xlabel("Time")

    ax.set_ylabel("dC(t)/dt")

    ax.set_title(model_name1)

    _style(ax, GetAdjustedFontSize())
    ax.figure.set_facecolor("white")

    return curve, np.array(params, dtype=float)


def plot_growth_illustration():
    """Four growth models with the same parameters, in a 2-by-2 panel (``plot_growth_illustration.m``).

    GLM, GRM, LM and Richards with r = 0.18, p = 0.9, a = 0.6, K = 10000,
    C(0) = 1 over 200 time units. Returns the figure.
    """
    plt = _plt()
    plt.close("all")

    fig = plt.figure(figsize=(14, 10), facecolor="white")

    legends = ["A", "B", "C", "D"]

    flag1s = [1, 2, 3, 4]

    for cc1 in range(4):

        ax = fig.add_subplot(2, 2, cc1 + 1)

        curve, _ = plotGrowthModel(flag1s[cc1], 0.18, 0.9, 0.6, 10000, 1, 200, ax=ax)

        ax.text(5, np.max(curve) - 20, legends[cc1], fontsize=24)

    fig.tight_layout()
    return fig
