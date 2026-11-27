"""Fit and forecast runners: rolling-window calibration with a parametric bootstrap.

Python port of ``Run_Fit_GrowthModels.m`` and ``Run_Forecasting_GrowthModels.m``.

The two MATLAB runners share almost all of their code; here both call one
implementation, :func:`_run`, and differ only where the MATLAB files differ.

Outputs are written to ``./output`` with the MATLAB file names: one MAT file
per window (``Forecast-growthModel-*.mat``, written with ``scipy.io.savemat``),
MAT files with the rolling summaries, and the CSV files
``parameters-``, ``MCSES-``, ``SCIS-``, ``diagnostics-``, ``AICcs-``,
``performance-calibration-`` and ``performance-forecasting-rollingwindow``.
"""

from __future__ import annotations

import os

import numpy as np

from . import io as gpio
from .aicc import get_nparams, getAICc, growthAICcCountAdjustment
from .errors import GrowthPredictError, warn
from .fitting import fit_model, initialParams
from .metrics import boundHitFractions, computeforecastperformance, computeQuantiles, computeWIS
from .noise import AddErrorStructure, estimatePearsonDispersion
from .options import options_fit, options_forecast, resolve_options
from .utils import adjustLastWindowStart, get_rng, load_timeseries, quantile

__all__ = ["Run_Fit_GrowthModels", "Run_Forecasting_GrowthModels"]

_PARAMLABELS = ["r", "p", "a", "K", "I0", "alpha", "d"]


def Run_Fit_GrowthModels(
    tstart1_pass=None,
    tend1_pass=None,
    windowsize1_pass=None,
    options_pass=None,
    fitopts_pass=None,
    *,
    input_dir="./input",
    output_dir="./output",
    plot=True,
    verbose=True,
    rng=None,
):
    """Fit a growth model to one or more calibration windows, with quantified uncertainty.

    Port of ``Run_Fit_GrowthModels.m``.

    Parameters
    ----------
    tstart1_pass, tend1_pass : int, optional
        One-based row where the first and the last window start. ``None`` uses
        the value of the options.
    windowsize1_pass : int, optional
        Calibration length.
    options_pass : function or Options, optional
        Options function (e.g. ``options_fit_mpox``) or :class:`Options`
        object; default ``options_fit``.
    fitopts_pass : dict, optional
        Optimizer settings added to the defaults, e.g.
        ``{'rBounds': [1e-3, 50]}`` to set the interval of r (see ``fit_model``).
    input_dir, output_dir : str
        Folders for the data file and the results (MATLAB: ``./input``, ``./output``).
    plot : bool
        Draw the figure of each window (figure number 100 + window start).
    verbose : bool
        Print the summary of the run.
    rng : numpy.random.Generator, optional
        Random stream; default is the toolbox stream (``growthpredict.rng(seed)``).

    Returns
    -------
    AICcs : ndarray, one row per window
        ``[window start, AICc, part1, part2, numparams]``.
    fit_model1 : ndarray
        Bootstrap model trajectories of the last window (one column per refit).
    """
    return _run(
        "fit", tstart1_pass, tend1_pass, windowsize1_pass, None, options_pass, fitopts_pass,
        input_dir, output_dir, plot, verbose, rng,
    )  # fmt: skip


def Run_Forecasting_GrowthModels(
    tstart1_pass=None,
    tend1_pass=None,
    windowsize1_pass=None,
    forecastingperiod_pass=None,
    options_pass=None,
    fitopts_pass=None,
    *,
    input_dir="./input",
    output_dir="./output",
    plot=True,
    verbose=True,
    rng=None,
):
    """Fit a growth model and forecast, for one or more windows, with quantified uncertainty.

    Port of ``Run_Forecasting_GrowthModels.m``. The arguments are those of
    :func:`Run_Fit_GrowthModels` plus ``forecastingperiod_pass``, the forecast
    horizon; the default options are ``options_forecast``.

    Returns
    -------
    AICcs : ndarray, one row per window
        ``[window start, AICc, part1, part2, numparams]``.
    forecast_model1 : ndarray
        Bootstrap model trajectories (calibration and forecast rows) of the
        last window, one column per refit.
    """
    return _run(
        "forecast", tstart1_pass, tend1_pass, windowsize1_pass, forecastingperiod_pass, options_pass,
        fitopts_pass, input_dir, output_dir, plot, verbose, rng,
    )  # fmt: skip


def _ci(P, Phatss, k):
    """[estimate, 2.5% and 97.5% quantile of the bootstrap refits] of parameter k."""
    return np.array([P[k], quantile(Phatss[:, k], 0.025), quantile(Phatss[:, k], 0.975)])


def _run(
    mode, tstart1_pass, tend1_pass, windowsize1_pass, forecastingperiod_pass, options_pass, fitopts_pass,
    input_dir, output_dir, plot, verbose, rng,
):  # fmt: skip

    forecast_mode = mode == "forecast"

    rng = get_rng(rng)

    if plot:
        import matplotlib.pyplot as plt

        from .plotting import draw_window_figure

        plt.close("all")

    # <============================================================================>
    # <=================== Load parameter values supplied by user =================>
    # <============================================================================>

    opts = resolve_options(options_pass, options_forecast if forecast_mode else options_fit)

    # Results are written to ./output; create it when missing.
    os.makedirs(output_dir, exist_ok=True)

    # <============================================================================>
    # <================================ Datasets properties =======================>
    # <============================================================================>

    DT = 1

    caddisease = opts.caddisease
    datatype = opts.datatype

    # <=============================================================================>
    # <=========================== Parameter estimation ============================>
    # <=============================================================================>

    method1 = opts.method1  # Type of estimation method: 0 = LSQ

    d = 1

    dist1 = opts.dist1  # type of error structure:

    # LSQ=0,
    # MLE Poisson=1,
    # Pearson chi-squard=2,
    # MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
    # MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
    # MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

    numstartpoints = opts.numstartpoints  # Number of initial guesses for the multi-start optimization

    numstartpoints_boot = 1  # Box-wide start points for each bootstrap refit (a jittered copy of the best fit and the best fit itself are always added)

    fitopts = {"logscale": True, "useParallel": False}  # optimizer settings, see fit_model

    # Optional last input: a dict whose entries are added to (or replace) the optimizer
    # settings above, e.g. {'rBounds': [1e-3, 50]} to set the interval of r.
    if fitopts_pass is not None and (not isinstance(fitopts_pass, dict) or len(fitopts_pass) > 0):
        if not isinstance(fitopts_pass, dict):
            raise GrowthPredictError(
                "GrowthPredict:fitopts", "fitopts_pass must be a dict, e.g. {'rBounds': [1e-3, 50]}."
            )
        fitopts.update(fitopts_pass)

    M = int(opts.B)  # number of bootstrap realizations to characterize parameter uncertainty

    # <==============================================================================>
    # <============================== Growth model ==================================>
    # <==============================================================================>

    flag1 = opts.flag1  # growth model
    model_name1 = opts.model_name1

    fixI0 = opts.fixI0  # 0=Estimate the initial number of cases; 1 = Fix it to the first data point

    # <==============================================================================>
    # <======================== Load epidemic data ==================================>
    # <==============================================================================>

    data, cadfilename1 = load_timeseries(opts.cadfilename1, input_dir)

    # <==============================================================================>
    # <========================== Forecasting parameters ============================>
    # <==============================================================================>

    if forecast_mode:
        getperformance = bool(opts.getperformance)
        forecastingperiod = int(opts.forecastingperiod if forecastingperiod_pass is None else forecastingperiod_pass)
    else:
        getperformance = True
        forecastingperiod = 0

    # <==================================================================================>
    # <========================== Parameters of the rolling window analysis =============>
    # <==================================================================================>

    tstart1 = int(opts.tstart1 if tstart1_pass is None else tstart1_pass)
    tend1 = int(opts.tend1 if tend1_pass is None else tend1_pass)
    windowsize1 = int(opts.windowsize1 if windowsize1_pass is None else windowsize1_pass)

    # <===========================================================================================================>
    # <====== Check that the number of estimated parameters is smaller than the number of data points ============>
    # <===========================================================================================================>

    numparams = get_nparams(method1, dist1, flag1, fixI0)

    if forecast_mode:
        if numparams >= windowsize1:
            raise GrowthPredictError(
                "GrowthPredict:calibrationPeriod",
                "Number of estimated parameters should be smaller than the calibration period. "
                "Consider increasing the length of the calibration period.",
            )
    elif numparams > windowsize1:
        raise GrowthPredictError(
            "GrowthPredict:calibrationPeriod",
            "Number of estimated parameters should not exceed the length of the calibration period. "
            "Consider increasing the length of the calibration period.",
        )

    # <==================================================================================>
    # ============================ Rolling window analysis ==============================>
    # <==================================================================================>

    acc = {
        name: []
        for name in (
            "param_rs", "param_as", "param_ps", "param_Ks", "param_I0s", "param_alphas", "param_ds",
            "MCSES", "SCIs", "param_bootmedians", "estimOutsideCIs", "fracLBs", "fracUBs",
            "startAgreements", "bestFitOnBounds", "RMSECSS", "MSECSS", "MAECSS", "PICSS", "MISCSS",
            "RMSEFSS", "MSEFSS", "MAEFSS", "PIFSS", "MISFSS", "WISCSS", "WISFSS", "quantilescs",
            "quantilesfs", "AICcs",
        )
    }  # fmt: skip

    nrows = data.shape[0]

    tend1, tend1Adjusted = adjustLastWindowStart(tend1, windowsize1, nrows)

    if tend1Adjusted and verbose:
        print("adjusting tend1")

    if tstart1 > tend1:
        warn(
            "GrowthPredict:windows",
            "incompatible tstart and tend values: the first window starts at row %d, the last possible "
            "window of %d observations starts at row %d. Nothing was run." % (tstart1, windowsize1, tend1),
        )
        return np.zeros((0, 5)), np.zeros((0, 0))

    names = gpio.OutputNames(
        output_dir, cadfilename1, flag1, fixI0, method1, dist1, windowsize1, forecastingperiod, caddisease, datatype
    )

    windows = list(range(tstart1, tend1 + 1))

    # Bootstrap counts across all rolling windows (reporting only).
    requestedBootstrap = M * len(windows)
    acceptedBootstrap = 0

    fit_model1 = np.zeros((0, 0))
    forecast_model1 = np.zeros((0, 0))

    for i in windows:  # rolling window analysis

        if plot and (tend1 - tstart1) > 10:
            plt.close("all")

        t_window = np.arange(i - 1, i - 1 + windowsize1)  # zero-based rows of the window

        # calibration data
        datac = data[t_window, 1].copy()
        timevect1 = DT * data[t_window, 0]

        # calibration and future data
        data_all = data[i - 1 :, 1].copy()
        timevect_all = DT * data[i - 1 :, 0]

        # first data point cannot be zero
        if datac[0] == 0:
            warn(
                "GrowthPredict:FirstPointZero",
                "Window %d: the first data point is zero; it is replaced by 0.01 for the fit. "
                "Consider fixI0 = 0 or a later tstart1." % i,
            )
            datac[0] = 0.01

        if forecast_mode and getperformance and data_all.size < windowsize1 + forecastingperiod:
            raise GrowthPredictError(
                "GrowthPredict:seriesTooShort",
                "Length of time series data is too short to evaluate the forecasting period indicated in "
                "<forecastingperiod> (%d observations from the window start, %d needed). "
                "Set getperformance = 0 to forecast beyond the data."
                % (data_all.size, windowsize1 + forecastingperiod),
            )

        data1 = np.column_stack([data[t_window, 0], datac])

        params0 = initialParams(datac, flag1)

        fit_d = fit_model(data1, params0, fixI0, numstartpoints, DT, flag1, 0, fitopts, method1=method1, rng=rng)

        P_model1d = fit_d.P
        fitcurve_model1d = fit_d.fitcurve
        fval = fit_d.fval
        exitflag_model1d = fit_d.flagg
        fitDiagnostics_model1d = fit_d.fitDiagnostics
        initialguess = fit_d.initialguess

        # Count likelihoods: complete fval to the full log-likelihood of the informative
        # observations so that AICc is comparable across estimation methods (zero for least squares).
        fvalAdjustAICc, _ = growthAICcCountAdjustment(method1, data1[:, 1], fixI0, P_model1d[5], P_model1d[6])

        AICc, part1, part2, numparams, nEff = getAICc(
            method1, dist1, flag1, fixI0, fval, data1.shape[0], fvalAdjustAICc
        )

        acc["AICcs"].append([i, AICc, part1, part2, numparams])

        poissonFallback = False  # set when dist1=2 is requested but the data are not over-dispersed

        if method1 == 0 and dist1 == 1:

            factor1 = 1.0

        elif method1 == 0 and dist1 == 2:  # estimate <factor1>, the variance to mean ratio that scales the NB error structure

            # Pearson dispersion of the data around the FITTED curve (quasi-Poisson estimator),
            # so that the epidemic trend is not counted as noise.
            # numparams also counts the Gaussian variance of the least-squares fit, hence the -1.
            factor1, dispInfo = estimatePearsonDispersion(data1[:, 1], fitcurve_model1d, numparams - 1, 1, fixI0)

            if not np.isfinite(factor1):
                warn(
                    "GrowthPredict:Dispersion",
                    "Window %d: the dispersion cannot be estimated (%d usable points for %d curve parameters). "
                    "Using Poisson bootstrap noise for this window." % (i, dispInfo.nUsed, numparams - 1),
                )

            # VAR = factor1*mean needs factor1 > 1. Not over-dispersed -> Poisson noise for THIS
            # window's bootstrap only (dist1_boot); dist1 and the output file names are unchanged.
            if not np.isfinite(factor1) or factor1 <= 1:
                poissonFallback = True
                factor1 = 1.0

        elif method1 == 0 and dist1 == 0:  # normal distribution of the error structure

            with np.errstate(all="ignore"):
                var1 = np.sum((fitcurve_model1d - datac) ** 2) / np.float64(fitcurve_model1d.size - numparams)
                factor1 = float(np.sqrt(var1))

        elif method1 == 1:

            dist1 = 1
            factor1 = 1.0

        elif method1 == 3:

            dist1 = 3
            factor1 = P_model1d[5]  # alpha; VAR=mean+alpha*mean

        elif method1 == 4:

            dist1 = 4
            factor1 = P_model1d[5]  # alpha; VAR=mean+alpha*mean^2

        elif method1 == 5:

            dist1 = 5
            factor1 = P_model1d[5]  # alpha; VAR=mean+alpha*mean^d
            d = P_model1d[6]

        else:
            raise GrowthPredictError(
                "GrowthPredict:errorStructure",
                "dist1 = %g cannot be combined with method1 = %g. With least squares (method1 = 0) "
                "use dist1 = 0, 1 or 2." % (dist1, method1),
            )

        names.dist1 = dist1

        # <===========================================================================================>
        # <======= Derive parameter uncertainty of the best fitting models and save results ==========>
        # <===========================================================================================>

        ntot = windowsize1 + max(forecastingperiod, 0)

        fit_model1 = np.zeros((windowsize1, M))
        forecast_model1 = np.zeros((ntot, M))
        forecast_model12 = np.zeros((ntot, 20 * M))

        Phatss_model1 = np.zeros((M, 7))

        dist1_boot = dist1  # noise model of this window's bootstrap
        if poissonFallback:
            dist1_boot = 1

        # Every bootstrap refit uses the parameter bounds of the fit to the observed data and
        # starts from that fit (see fit_model, 'Bounds of bootstrap refits').
        fitopts_boot = dict(fitopts)
        fitopts_boot["LB"] = fitDiagnostics_model1d.LB
        fitopts_boot["UB"] = fitDiagnostics_model1d.UB

        bootLB = np.full((M, 7), np.nan)  # optimizer box of every bootstrap refit
        bootUB = np.full((M, 7), np.nan)

        f_model1_sims = np.zeros((windowsize1, M))

        cumfit = np.cumsum(fitcurve_model1d)

        timevect2 = timevect1

        for j in range(M):

            f_model1_sim = AddErrorStructure(cumfit, 1, dist1_boot, factor1, d, rng=rng)[:, 0]

            f_model1_sims[:, j] = f_model1_sim

            # Fit the model to the bootstrap data
            data1b = np.column_stack([timevect1, f_model1_sim])

            fit_b = fit_model(
                data1b, P_model1d, fixI0, numstartpoints_boot, DT, flag1, forecastingperiod, fitopts_boot,
                method1=method1, rng=rng,
            )  # fmt: skip

            P_model1 = fit_b.P
            forecastcurve_model1 = fit_b.forecastcurve
            timevect2 = fit_b.timevect2

            bootLB[j, :] = fit_b.fitDiagnostics.LB
            bootUB[j, :] = fit_b.fitDiagnostics.UB
            acceptedBootstrap += int(fit_b.flagg > 0 and np.isfinite(fit_b.fval))

            fit_model1[:, j] = fit_b.fitcurve

            forecast_model1[:, j] = forecastcurve_model1

            if method1 == 0:

                noisy = AddErrorStructure(np.cumsum(forecastcurve_model1), 20, dist1_boot, factor1, 0, rng=rng)

            else:

                noisy = AddErrorStructure(
                    np.cumsum(forecastcurve_model1), 20, dist1_boot, P_model1[-2], P_model1[-1], rng=rng
                )

            forecast_model12[:, 20 * j : 20 * (j + 1)] = noisy

            Phatss_model1[j, :] = P_model1

        # end bootstrapping loop

        datalatest = data[i - 1 :, 0:2]

        # <==================================================================================================>
        # <========== Get forecast performance metrics for the model =======================================>
        # <==================================================================================================>

        (RMSECS_model1, MSECS_model1, MAECS_model1, PICS_model1, MISCS_model1,
         RMSEFS_model1, MSEFS_model1, MAEFS_model1, PIFS_model1, MISFS_model1) = computeforecastperformance(
            data1, datalatest, forecast_model1, forecast_model12, forecastingperiod
        )  # fmt: skip

        WISC_model1, WISFS_model1, _, _ = computeWIS(data1, datalatest, forecast_model12, forecastingperiod)

        # store metrics for calibration
        acc["RMSECSS"].append(RMSECS_model1[-1, -1])
        acc["MSECSS"].append(MSECS_model1[-1, -1])
        acc["MAECSS"].append(MAECS_model1[-1, -1])
        acc["PICSS"].append(PICS_model1[-1, -1])
        acc["MISCSS"].append(MISCS_model1[-1, -1])
        acc["WISCSS"].append(WISC_model1)

        # store metrics for short-term forecasts
        if forecastingperiod > 0 and RMSEFS_model1.shape[0] > 0:

            acc["RMSEFSS"].append(RMSEFS_model1[-1, -1])
            acc["MSEFSS"].append(MSEFS_model1[-1, -1])
            acc["MAEFSS"].append(MAEFS_model1[-1, -1])
            acc["PIFSS"].append(PIFS_model1[-1, -1])
            acc["MISFSS"].append(MISFS_model1[-1, -1])
            acc["WISFSS"].append(WISFS_model1[-1, -1])

        # <==================================================================================================>
        # <========== Compute quantiles of the calibration and forecasting periods and store ================>
        # <==================================================================================================>

        quantilesc, quantilesf, _ = computeQuantiles(data1, forecast_model12, forecastingperiod)

        acc["quantilescs"].append(quantilesc)
        acc["quantilesfs"].append(quantilesf)

        # <========================================================================================>
        # <================================ Parameter estimates ===================================>
        # <========================================================================================>

        # Point estimate = fit to the observed data (P_model1d); the 95% CI is the 2.5-97.5%
        # percentile range of the bootstrap refits. The bootstrap medians are kept in
        # param_bootmedians.
        param_r = _ci(P_model1d, Phatss_model1, 0)
        param_p = _ci(P_model1d, Phatss_model1, 1)
        param_a = _ci(P_model1d, Phatss_model1, 2)
        param_K = _ci(P_model1d, Phatss_model1, 3)
        param_I0 = _ci(P_model1d, Phatss_model1, 4)
        param_alpha = _ci(P_model1d, Phatss_model1, 5)
        param_d = _ci(P_model1d, Phatss_model1, 6)

        # <========================================================================================>
        # <=========================== Diagnostics of the estimates ===============================>
        # <========================================================================================>

        # A best fit on a bound means the bound, not the data, determined that estimate.
        bestFitOnBound = fitDiagnostics_model1d.atLB | fitDiagnostics_model1d.atUB

        if np.any(bestFitOnBound):
            warn(
                "GrowthPredict:BestFitOnBound",
                "Window %d: the fitted value of %s is on its parameter bound "
                "(bounds are set in fit_model from the data)." % (i, _join(bestFitOnBound)),
            )

        param_bootmedian = np.median(Phatss_model1, axis=0)

        # A fitted value outside its own bootstrap interval signals bootstrap bias (often an
        # estimate on a bound or along a ridge), not a coding error: inspect the fit.
        estimOutsideCI = (P_model1d < quantile(Phatss_model1, 0.025, axis=0)) | (
            P_model1d > quantile(Phatss_model1, 0.975, axis=0)
        )

        if np.any(estimOutsideCI):
            warn(
                "GrowthPredict:EstimateOutsideCI",
                "Window %d: the fitted value of %s lies outside its bootstrap 95%% interval. "
                "Inspect the fit before reporting." % (i, _join(estimOutsideCI)),
            )

        # Share of bootstrap estimates sitting on the optimizer bounds (NaN for fixed parameters).
        # A large share means the interval is set by the box rather than by the data.
        fracLB, fracUB, _, _ = boundHitFractions(Phatss_model1, bootLB, bootUB)

        with np.errstate(invalid="ignore"):
            boundLimited = (fracLB > 0.05) | (fracUB > 0.05)

        if np.any(boundLimited):
            warn(
                "GrowthPredict:BoundLimited",
                "Window %d: more than 5%% of the bootstrap estimates of %s are on a parameter bound; "
                "the confidence interval is limited by the bounds set in fit_model." % (i, _join(boundLimited)),
            )

        if fitDiagnostics_model1d.startAgreement < 0.5:
            warn(
                "GrowthPredict:StartAgreement",
                "Window %d: only %.0f%% of the converged optimizer starts reached the best objective value; "
                "the likelihood surface may be multimodal. Consider increasing numstartpoints."
                % (i, 100 * fitDiagnostics_model1d.startAgreement),
            )

        acc["param_bootmedians"].append(param_bootmedian)
        acc["estimOutsideCIs"].append(estimOutsideCI.astype(float))
        acc["fracLBs"].append(fracLB)
        acc["fracUBs"].append(fracUB)
        acc["startAgreements"].append(fitDiagnostics_model1d.startAgreement)
        acc["bestFitOnBounds"].append(bestFitOnBound.astype(float))

        MCSE = (np.std(Phatss_model1, axis=0, ddof=1) if M > 1 else np.zeros(7)) / np.sqrt(M)

        params_ci = np.vstack([param_r, param_p, param_a, param_K, param_I0, param_alpha, param_d])
        with np.errstate(all="ignore"):
            SCI = np.log10(params_ci[:, 2] / params_ci[:, 1])

        acc["param_rs"].append(param_r)
        acc["param_as"].append(param_a)
        acc["param_ps"].append(param_p)
        acc["param_Ks"].append(param_K)
        acc["param_I0s"].append(param_I0)
        acc["param_alphas"].append(param_alpha)
        acc["param_ds"].append(param_d)

        acc["MCSES"].append(MCSE)
        acc["SCIs"].append(SCI)

        # <========================================================================================>
        # <======================= Plot parameter distributions and the fit / forecast ============>
        # <========================================================================================>

        if plot:
            draw_window_figure(
                100 + i,
                Phatss_model1=Phatss_model1,
                params=(param_r, param_p, param_a, param_K),
                timevect1=timevect1,
                timevect2=timevect2,
                forecast_model12=forecast_model12,
                fit_model1=fit_model1,
                timevect_all=timevect_all,
                data_all=data_all,
                forecastingperiod=forecastingperiod,
                ylabel="%s %s" % (caddisease, datatype),
                model_name1=model_name1,
                all_ylabels=True,
            )

        # <=========================================================================================>
        # <================================ Save the results of this window ========================>
        # <=========================================================================================>

        summary = {key: _stack(val) for key, val in acc.items()}

        window = dict(
            # settings
            cadfilename1=cadfilename1, caddisease=caddisease, datatype=datatype, method1=method1, dist1=dist1,
            dist1_boot=dist1_boot, poissonFallback=poissonFallback, numstartpoints=numstartpoints,
            numstartpoints_boot=numstartpoints_boot, M=M, flag1=flag1, model_name1=model_name1, fixI0=fixI0,
            getperformance=int(getperformance), forecastingperiod=forecastingperiod, windowsize1=windowsize1,
            tstart1=tstart1, tend1=tend1, i=i, DT=DT, d=d, factor1=factor1,
            # data
            data=data, datac=_col(datac), data1=data1, data_all=_col(data_all), datalatest=datalatest,
            timevect1=_col(timevect1), timevect2=_col(timevect2), timevect_all=_col(timevect_all),
            # fit to the observed data
            params0=params0, P_model1d=P_model1d, fitcurve_model1d=_col(fitcurve_model1d), fval=fval,
            exitflag_model1d=exitflag_model1d, fitDiagnostics_model1d=fitDiagnostics_model1d,
            initialguess=initialguess, AICc=AICc, part1=part1, part2=part2, numparams=numparams, nEff=nEff,
            fvalAdjustAICc=fvalAdjustAICc,
            # bootstrap
            Phatss_model1=Phatss_model1, bootLB=bootLB, bootUB=bootUB, f_model1_sims=f_model1_sims,
            fit_model1=fit_model1, forecast_model1=forecast_model1, forecast_model12=forecast_model12,
            # estimates and diagnostics of this window
            param_r=param_r, param_p=param_p, param_a=param_a, param_K=param_K, param_I0=param_I0,
            param_alpha=param_alpha, param_d=param_d, param_bootmedian=param_bootmedian, MCSE=MCSE, SCI=SCI,
            estimOutsideCI=estimOutsideCI, fracLB=fracLB, fracUB=fracUB, bestFitOnBound=bestFitOnBound,
            # performance of this window
            RMSECS_model1=RMSECS_model1, MSECS_model1=MSECS_model1, MAECS_model1=MAECS_model1,
            PICS_model1=PICS_model1, MISCS_model1=MISCS_model1, WISC_model1=WISC_model1,
            RMSEFS_model1=RMSEFS_model1, MSEFS_model1=MSEFS_model1, MAEFS_model1=MAEFS_model1,
            PIFS_model1=PIFS_model1, MISFS_model1=MISFS_model1, WISFS_model1=WISFS_model1,
            quantilesc=quantilesc, quantilesf=quantilesf,
            # rolling summaries up to this window
            **summary,
        )  # fmt: skip

        gpio.save_mat(names.window_mat(i, tend1), window)

    # end rolling window analysis

    summary = {key: _stack(val) for key, val in acc.items()}

    AICcs = summary["AICcs"]

    # save model parameters from tstart1 to tend1
    gpio.save_mat(
        names.summary_mat("parameters", tstart1, tend1),
        {
            key: summary[key]
            for key in (
                "param_rs", "param_as", "param_ps", "param_Ks", "param_I0s", "param_alphas", "param_ds",
                "param_bootmedians", "estimOutsideCIs", "fracLBs", "fracUBs", "startAgreements", "bestFitOnBounds",
            )
        },
    )  # fmt: skip

    # save calibration performance metrics from tstart1 to tend1 (AICc, MSE, MAE, coverage 95% PI, MIS, WIS)
    gpio.save_mat(
        names.summary_mat("performanceCalibration", tstart1, tend1),
        {key: summary[key] for key in ("AICcs", "RMSECSS", "MSECSS", "MAECSS", "PICSS", "MISCSS", "WISCSS")},
    )

    # save quantiles of the model fits from tstart1 to tend1
    gpio.save_mat(names.summary_mat("QuantilesCalibration", tstart1, tend1), {"quantilescs": summary["quantilescs"]})

    if forecastingperiod > 0:

        # save forecasting performance metrics from tstart1 to tend1
        gpio.save_mat(
            names.summary_mat("performanceForecasting", tstart1, tend1),
            {key: summary[key] for key in ("RMSEFSS", "MSEFSS", "MAEFSS", "PIFSS", "MISFSS", "WISFSS")},
        )

        # save quantiles of the model forecasts from tstart1 to tend1
        gpio.save_mat(
            names.summary_mat("QuantilesForecastingPerformance", tstart1, tend1),
            {"quantilesfs": summary["quantilesfs"]},
        )

    # <=====================================================================================================>
    # <============================== Save the CSV files of the rolling window analysis ====================>
    # <=====================================================================================================>

    gpio.write_rolling_tables(names, tstart1, tend1, method1, summary)

    gpio.write_diagnostics_table(names, tstart1, tend1, summary)

    gpio.write_performance_calibration(names, tstart1, tend1, windowsize1, summary)

    if forecast_mode and getperformance and forecastingperiod > 0 and summary["MAEFSS"].size > 0:
        gpio.write_performance_forecasting(names, tstart1, tend1, forecastingperiod, summary)

    if verbose:
        print_options_summary(
            "options_forecast" if forecast_mode else "options_fit", cadfilename1, caddisease, datatype, method1,
            dist1, numstartpoints, M, flag1, model_name1, fixI0, windowsize1, tstart1, tend1,
            forecasting=(getperformance, forecastingperiod) if forecast_mode else None,
        )  # fmt: skip

        # Print counts only; they do not change fitting or retained results.
        print("\nRequested bootstrap datasets: %d" % requestedBootstrap)
        print("Accepted bootstrap datasets: %d" % acceptedBootstrap)
        print("Unresolved bootstrap datasets: %d\n" % (requestedBootstrap - acceptedBootstrap))

    return AICcs, (forecast_model1 if forecast_mode else fit_model1)


# <============================================================================>
# <================================ Helpers ===================================>
# <============================================================================>


def _join(mask):
    return ", ".join(name for name, hit in zip(_PARAMLABELS, mask) if hit)


def _col(x):
    return np.asarray(x, dtype=float).reshape(-1, 1)


def _stack(rows):
    """Rows collected over the windows -> one array (windows down the rows)."""
    if len(rows) == 0:
        return np.zeros((0, 0))
    first = np.asarray(rows[0])
    if first.ndim == 0:
        return np.asarray(rows, dtype=float).reshape(-1, 1)
    if first.ndim == 1:
        return np.vstack([np.asarray(r, dtype=float) for r in rows])
    return np.vstack(rows)  # blocks of quantiles


def print_options_summary(
    title, cadfilename1, caddisease, datatype, method1, dist1, numstartpoints, M, flag1, model_name1, fixI0,
    windowsize1, tstart1, tend1, forecasting=None,
):  # fmt: skip
    """Print the settings of a run, as the MATLAB runners and plotting functions do."""
    separator = "<=============================================================================>"

    print(separator)
    print("**************************** %s ***********************************" % title)
    print(separator)

    print("<================================ Datasets properties ========================>")
    print(separator)

    print("cadfilename1: %-30s\t(String) Name of the data file containing time-series data." % cadfilename1)
    print("caddisease:    %-30s\t(String) Name of the disease or subject related to the time-series data." % caddisease)
    print("datatype:      %-30s\t(String) Nature of the data (cases, deaths, hospitalizations, etc)." % datatype)

    print(separator)

    print("<=========================== Parameter estimation ============================>")
    print(separator)

    print("method1:        %-20s\t(Method for optimization)." % method1)
    print("dist1:          %-20s\t(Distribution used in the estimation)." % dist1)
    print("numstartpoints: %-20s\t(Number of initial guesses for optimization)." % numstartpoints)
    print("B:              %-20s\t(Number of bootstrap realizations for uncertainty)." % M)

    print(separator)

    print("<=========================== Growth model ====================================>")
    print(separator)

    print("flag1:       %-20s\t(Integer) Growth model to fit the time-series data." % flag1)
    print("model_name1: %-20s\t(String) Name of the model." % model_name1)
    print("fixI0:       %-20s\t(Boolean) Fix initial value in time-series (true/false)." % fixI0)

    print(separator)

    if forecasting is not None:
        getperformance, forecastingperiod = forecasting
        print("<============================================================================>")
        print("                          Forecasting Parameters                             ")
        print("<============================================================================>")
        print("  - Forecasting Performance Metrics: %s" % ("Enabled" if getperformance else "Disabled"))
        print("  - Forecast Horizon: %s time units ahead" % forecastingperiod)
        print("<============================================================================>")

    print("<======================= Sliding window parameters ===========================>")
    print(separator)

    print("windowsize1: %-20s\t(Integer) Moving window size." % windowsize1)
    print("tstart11:    %-20s\t(Integer) Start time point for rolling window analysis." % tstart1)
    print("tend1:       %-20s\t(Integer) End time point for rolling window analysis." % tend1)

    print(separator)
