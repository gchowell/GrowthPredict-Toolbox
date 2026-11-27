"""End-to-end tests of the runners and the plotting functions.

Each test runs the bootstrap with a small number of replicates, so the whole
file takes a few minutes. Deselect with ``-m "not slow"``.

Targets are the README of the toolbox (file names, table layouts, the
quick-start numbers) and the MATLAB values recorded by the maintainer in
docs/validation/README.md for the U.S. mpox series.
"""

import math
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pytest

import growthpredict as gp
from growthpredict import io as gpio

pytestmark = pytest.mark.slow

QUICK = "flag1-1-fixI0-1-method-0-dist-0"


def run_quiet(function, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return function(*args, verbose=False, **kwargs)


def table(folder, name):
    path = os.path.join(folder, name)
    assert os.path.isfile(path), "missing %s; folder holds %s" % (name, sorted(os.listdir(folder)))
    names, array = gpio.readtable(path)
    return names, array


# ---------------------------------------------------------------- quick start


@pytest.fixture(scope="module")
def quick(tmp_path_factory):
    """README quick start: rows 1-20 of the U.S. mpox series, 4-step forecast (60 replicates)."""
    out = str(tmp_path_factory.mktemp("quick"))
    gp.rng(1)
    opts = gp.options_forecast().replace(B=60)
    AICcs, forecast_model1 = run_quiet(gp.Run_Forecasting_GrowthModels, 1, 1, 20, 4, opts, output_dir=out)
    figs_after_run = plt.get_fignums()
    plt.close("all")
    return dict(out=out, opts=opts, AICcs=AICcs, forecast_model1=forecast_model1, figs=figs_after_run)


def test_quick_start_returns_and_files(quick):
    out = quick["out"]
    assert quick["AICcs"].shape == (1, 5)
    assert quick["forecast_model1"].shape == (24, 60)
    assert quick["figs"] == [101]  # the runner draws figure 100 + window start

    tail = "tstart-1-tend-1-calibrationperiod-20-horizon-4-Mpox-cases.csv"
    mat = "Most_Recent_Timeseries_US-CDC.txt-%s-tstart-1-tend-1-calibrationperiod-20-forecastingperiod-4.mat" % QUICK
    expected = [
        "Forecast-growthModel-" + mat,
        "parameters-growthModel-" + mat,
        "performanceCalibration-growthModel-" + mat,
        "performanceForecasting-growthModel-" + mat,
        "QuantilesCalibration-growthModel-" + mat,
        "QuantilesForecastingPerformance-growthModel-" + mat,
        "parameters-rollingwindow-%s-%s" % (QUICK, tail),
        "MCSES-rollingwindow-%s-%s" % (QUICK, tail),
        "SCIS-rollingwindow-%s-%s" % (QUICK, tail),
        "diagnostics-rollingwindow-%s-%s" % (QUICK, tail),
        "AICcs-rollingwindow-%s-%s" % (QUICK, tail),
        "performance-calibration-%s-%s" % (QUICK, tail),
        "performance-forecasting-%s-%s" % (QUICK, tail),
    ]
    assert sorted(os.listdir(out)) == sorted(expected)


def test_quick_start_estimates(quick):
    out = quick["out"]
    names, params = table(out, "parameters-rollingwindow-%s-tstart-1-tend-1-calibrationperiod-20-horizon-4-Mpox-cases.csv" % QUICK)
    assert names == ["time", "r estimate", "r LB", "r UB", "p estimate", "p LB", "p UB", "a estimate", "a LB", "a UB",
                     "K0 estimate", "K0 LB", "K0 UB", "I0 estimate", "I0 LB", "I0 UB"]  # fmt: skip
    row = dict(zip(names, params[0]))

    # the fit to the observed data (README figure: r = 1.73, p = 0.853, K = 28600; MATLAB test: 1.7302, 0.8526, 28547)
    np.testing.assert_allclose(row["r estimate"], 1.7302, rtol=0.02)
    assert abs(row["p estimate"] - 0.8526) < 0.005
    np.testing.assert_allclose(row["K0 estimate"], 28547, rtol=0.01)
    assert row["time"] == 1 and row["a estimate"] == 1 and row["I0 estimate"] == 14

    # bootstrap percentile limits (MATLAB: r 1.47-2.15 with B = 100, 1.49-2.08 with B = 300)
    assert 1.30 < row["r LB"] < 1.60 and 1.95 < row["r UB"] < 2.35
    assert 0.80 < row["p LB"] < 0.845 and 0.862 < row["p UB"] < 0.90
    assert 26300 < row["K0 LB"] < 28000 and 29200 < row["K0 UB"] < 31200
    assert row["a LB"] == row["a UB"] == 1 and row["I0 LB"] == row["I0 UB"] == 14

    # AICc of the least-squares fit: Gaussian likelihood on n - 1 = 19 points, k = 4
    i, AICc, part1, part2, k = quick["AICcs"][0]
    sse = 3.1779e5
    assert (i, k) == (1, 4)
    np.testing.assert_allclose(part1, 19 * math.log(sse / 19) + 19 * math.log(2 * math.pi) + 19, rtol=1e-3)
    assert math.isclose(part2, 8 + 40 / 14) and math.isclose(AICc, part1 + part2)
    _, aicc_csv = table(out, "AICcs-rollingwindow-%s-tstart-1-tend-1-calibrationperiod-20-horizon-4-Mpox-cases.csv" % QUICK)
    np.testing.assert_allclose(aicc_csv, quick["AICcs"], rtol=1e-12)

    names, perf = table(out, "performance-forecasting-%s-tstart-1-tend-1-calibrationperiod-20-horizon-4-Mpox-cases.csv" % QUICK)
    assert names == ["time", "Horizon", "MAE", "MSE", "Coverage 95%PI", "WIS"]
    assert perf.shape == (1, 6) and perf[0, 1] == 4 and 0 <= perf[0, 4] <= 100
    assert perf[0, 2] > 0 and perf[0, 3] >= perf[0, 2] ** 2 * 0.999 and perf[0, 5] > 0

    names, diag = table(out, "diagnostics-rollingwindow-%s-tstart-1-tend-1-calibrationperiod-20-horizon-4-Mpox-cases.csv" % QUICK)
    assert len(names) == 2 + 5 * 7 and names[1] == "startAgreement" and names[-1] == "d best fit on bound"
    assert 0.5 <= diag[0, 1] <= 1


def test_quick_start_mat_file(quick, mpox):
    W = gpio.load_mat(
        os.path.join(
            quick["out"],
            "Forecast-growthModel-Most_Recent_Timeseries_US-CDC.txt-%s-tstart-1-tend-1-calibrationperiod-20-forecastingperiod-4.mat" % QUICK,
        )
    )
    # the variables named in the README
    assert W.P_model1d.shape == (1, 7)
    assert W.Phatss_model1.shape == (60, 7)
    assert W.forecast_model12.shape == (24, 20 * 60)  # 20 noise draws per refit
    assert W.forecast_model1.shape == (24, 60) and W.fit_model1.shape == (20, 60)
    np.testing.assert_allclose(W.timevect2.ravel(), np.arange(24))
    assert gpio.scalar(W.dist1_boot) == 0 and gpio.scalar(W.poissonFallback) == 0
    diag = W.fitDiagnostics_model1d[0, 0]
    assert diag.LB.shape == (1, 7) and diag.UB.shape == (1, 7)
    np.testing.assert_allclose(diag.UB.ravel()[0], 10 * np.max(mpox[:20, 1]))  # r in [1e-4, 10*max(y)]
    np.testing.assert_allclose(W.forecast_model1, quick["forecast_model1"])
    # every refit starts at the fixed first observation
    assert np.all(W.forecast_model1[0, :] == 14) and np.all(W.forecast_model12[0, :] == 14)


def test_plot_forecast_and_trajectory_table(quick):
    out = quick["out"]
    param_doubling, seq_doublingtimes, forecast_model1 = run_quiet(
        gp.plotForecast_GrowthModels, 1, 1, 20, 4, quick["opts"], output_dir=out
    )
    assert plt.get_fignums() == [101, 201]
    np.testing.assert_allclose(forecast_model1, quick["forecast_model1"])
    assert param_doubling.shape == (3,) and seq_doublingtimes.shape[1] == 5
    assert seq_doublingtimes[0, 4] == 1.0 and np.all(np.diff(seq_doublingtimes[:, 0]) == 1)

    # README: 'Open the trajectory file to see the forecasts beside the four held-out counts (669, 546, 404 and 329)'
    names, traj = table(out, "Forecast-%s-tstart-1-calibrationperiod-20-horizon-4-Mpox-cases.csv" % QUICK)
    assert names == ["time", "data", "median", "LB", "UB"]
    assert traj.shape == (24, 5)
    assert traj[-4:, 0].tolist() == [20, 21, 22, 23]
    assert traj[-4:, 1].tolist() == [669, 546, 404, 329]
    assert np.all(traj[:, 3] <= traj[:, 2]) and np.all(traj[:, 2] <= traj[:, 4]) and np.all(traj[:, 3] >= 0)
    # the forecast follows the decline (MATLAB figure: about 470, 330, 235, 160 with limits a few hundred wide)
    assert np.all(np.diff(traj[-4:, 2]) < 0) and 380 < traj[-4, 2] < 560 and 100 < traj[-1, 2] < 230
    assert np.mean((traj[:20, 1] >= traj[:20, 3]) & (traj[:20, 1] <= traj[:20, 4])) >= 0.85

    names, q = table(out, "quantile-forecast-%s-tstart-1-calibrationperiod-20-horizon-4-Mpox-cases.csv" % QUICK)
    assert names[0] == "Q_0.010" and names[11] == "Q_0.500" and names[-1] == "Q_0.990" and q.shape == (24, 23)
    np.testing.assert_allclose(q[:, 1], traj[:, 3], rtol=1e-10)  # Q_0.025 is LB
    np.testing.assert_allclose(q[:, 21], traj[:, 4], rtol=1e-10)  # Q_0.975 is UB
    assert os.path.isfile(os.path.join(out, "quantile-flag1-1-tstart-1-fixI0-1-method-0-dist-0-calibrationperiod-20-horizon-4-Mpox-cases.csv"))
    assert os.path.isfile(os.path.join(out, "doublingtimes-flag1-1-tstart-1-fixI0-1-method-0-dist-0-calibrationperiod-20-horizon-4-Mpox-cases.csv"))

    # forecast with the reproduction number: R_t over calibration and forecast periods
    param_R, _, _ = run_quiet(gp.plotForecast_ReproductionNumber, 1, 1, 20, 4, quick["opts"], gp.options_Rt, output_dir=out)
    assert plt.get_fignums() == [101, 201, 501]
    assert param_R.shape == (1, 5) and param_R[0, 0] == 23 and 0.3 < param_R[0, 1] < 0.9
    names, rt = table(out, "Rt-flag1-1-tstart-1-fixI0-1-method-0-dist-0-calibrationperiod-20-horizon-4-Mpox-cases.csv")
    assert rt[:, 0].tolist() == list(range(3, 24))


def test_same_seed_same_results(tmp_path):
    opts = gp.options_forecast().replace(B=8)
    results = []
    for k in range(2):
        gp.rng(123)
        out = str(tmp_path / ("run%d" % k))
        _, fm = run_quiet(gp.Run_Forecasting_GrowthModels, 1, 1, 20, 4, opts, output_dir=out, plot=False)
        results.append(fm)
    assert np.array_equal(results[0], results[1])
    gp.rng(124)
    _, other = run_quiet(gp.Run_Forecasting_GrowthModels, 1, 1, 20, 4, opts, output_dir=str(tmp_path / "run2"), plot=False)
    assert not np.array_equal(results[0], other)


# ---------------------------------------------------------------- fit only, reproduction number, doubling times


@pytest.fixture(scope="module")
def last20(tmp_path_factory):
    """README: 'Run_Fit_GrowthModels(13, 13, 20, @options_fit_mpoxCDC)', the last 20 observations (100 replicates)."""
    out = str(tmp_path_factory.mktemp("last20"))
    gp.rng(1)
    opts = gp.options_fit_mpoxCDC().replace(B=100)
    AICcs, fit_model1 = run_quiet(gp.Run_Fit_GrowthModels, 13, 13, 20, opts, output_dir=out, plot=False)
    return dict(out=out, opts=opts, AICcs=AICcs, fit_model1=fit_model1)


def test_fit_of_the_last_20_observations(last20):
    out = last20["out"]
    assert last20["fit_model1"].shape == (20, 100)
    names, params = table(out, "parameters-rollingwindow-%s-tstart-13-tend-13-calibrationperiod-20-horizon-0-Mpox-cases.csv" % QUICK)
    row = dict(zip(names, params[0]))
    # MATLAB (maintainer, 25 Nov 2026): r, p, K = 801.8, 0.178, 18153; bootstrap r 8e2 (5.7e2, 1.2e3)
    np.testing.assert_allclose(row["r estimate"], 801.8, rtol=0.02)
    assert abs(row["p estimate"] - 0.178) < 0.005
    np.testing.assert_allclose(row["K0 estimate"], 18153, rtol=0.005)
    assert 450 < row["r LB"] < 680 and 1000 < row["r UB"] < 1500
    assert row["time"] == 13 and row["I0 estimate"] == 3027
    assert not os.path.exists(os.path.join(out, "performance-forecasting-%s-tstart-13-tend-13-calibrationperiod-20-horizon-0-Mpox-cases.csv" % QUICK))


def test_plot_fit_tables(last20, mpox):
    out = last20["out"]
    param_doubling, seq_doublingtimes, fit_model1 = run_quiet(gp.plotFit_GrowthModels, 13, 13, 20, last20["opts"], output_dir=out)
    assert plt.get_fignums() == [114]  # 101 + window start
    np.testing.assert_allclose(fit_model1, last20["fit_model1"])

    # MATLAB: doubling time of rows 13-32, count of the series: 15.61 (13.44, 17.32), reached by 100%
    assert seq_doublingtimes.shape == (1, 5)
    assert seq_doublingtimes[0, 0] == 1 and seq_doublingtimes[0, 4] == 1.0
    np.testing.assert_allclose(seq_doublingtimes[0, 1], 15.61, rtol=0.04)
    assert 12.3 < seq_doublingtimes[0, 2] < 14.8 and 16.3 < seq_doublingtimes[0, 3] < 18.6
    np.testing.assert_allclose(param_doubling[0], seq_doublingtimes[0, 1])

    names, fit = table(out, "Fit-flag1-1-tstart-13-fixI0-1-method-0-dist-0-calibrationperiod-20-horizon-0-Mpox-cases.csv")
    assert names == ["time", "data", "median", "LB", "UB"] and fit.shape == (20, 5)
    assert fit[:, 0].tolist() == list(range(12, 32)) and fit[:, 1].tolist() == mpox[12:32, 1].tolist()
    names, q = table(out, "quantile-fit-flag1-1-tstart-13-fixI0-1-method-0-dist-0-calibrationperiod-20-horizon-0-Mpox-cases.csv")
    assert q.shape == (20, 23)


def test_reproduction_number_of_the_last_20_observations(last20):
    out = last20["out"]
    param_R, param_doubling, seq_doublingtimes = run_quiet(
        gp.plotFit_ReproductionNumber, 13, 13, 20, last20["opts"], gp.options_Rt, output_dir=out
    )
    assert plt.get_fignums() == [114, 501]

    names, rt = table(out, "Rt-flag1-1-tstart-13-fixI0-1-method-0-dist-0-calibrationperiod-20-horizon-0-Mpox-cases.csv")
    assert names == ["time", "Rt median", "Rt 95%CI LB", "Rt 95% CI UB", "GI coverage"]

    # MATLAB (maintainer): 17 rows exported; first row (t = 15) 0.850 (0.837, 0.863), coverage 0.9767;
    # last row (t = 31) 0.721 (0.706, 0.735)
    assert rt.shape == (17, 5)
    assert rt[:, 0].tolist() == list(range(15, 32))
    np.testing.assert_allclose(rt[0, 1], 0.850, rtol=0.015)
    np.testing.assert_allclose(rt[0, 4], 0.9767, atol=2e-4)
    np.testing.assert_allclose(rt[-1, 1], 0.721, rtol=0.02)
    assert 0.82 < rt[0, 2] < rt[0, 1] < rt[0, 3] < 0.88
    assert 0.68 < rt[-1, 2] < rt[-1, 1] < rt[-1, 3] < 0.76
    assert np.all(rt[:, 4] >= 0.95)

    assert param_R.shape == (1, 5)
    np.testing.assert_allclose(param_R[0, :4], rt[-1, :4], rtol=1e-12)
    assert seq_doublingtimes.shape == (1, 5)

    # another generation interval: the cholera serial interval pulls R_t toward 1
    _, rt_default = rt, rt
    run_quiet(gp.plotFit_ReproductionNumber, 13, 13, 20, last20["opts"], gp.options_Rt_cholera, output_dir=out)
    _, rt_cholera = table(out, "Rt-flag1-1-tstart-13-fixI0-1-method-0-dist-0-calibrationperiod-20-horizon-0-Mpox-cases.csv")
    assert abs(rt_cholera[-1, 1] - 1) < abs(rt_default[-1, 1] - 1)

    # a fixed interval too long for the window: every estimate is kept out, one flagged row remains
    with pytest.warns(gp.GrowthPredictWarning, match="GrowthPredict:Rt:coverage"):
        gp.plotFit_ReproductionNumber(13, 13, 20, last20["opts"], lambda: (1, 30.0, 100.0), output_dir=out, verbose=False)
    _, rt_long = table(out, "Rt-flag1-1-tstart-13-fixI0-1-method-0-dist-0-calibrationperiod-20-horizon-0-Mpox-cases.csv")
    assert rt_long.shape == (1, 5) and rt_long[0, 0] == 31 and rt_long[0, 4] < 0.95


# ---------------------------------------------------------------- rolling windows


@pytest.fixture(scope="module")
def rolling(tmp_path_factory):
    """Three rolling windows (start rows 1-3), 4-step forecasts, 20 replicates."""
    out = str(tmp_path_factory.mktemp("rolling"))
    gp.rng(1)
    opts = gp.options_forecast().replace(B=20)
    AICcs, _ = run_quiet(gp.Run_Forecasting_GrowthModels, 1, 3, 20, 4, opts, output_dir=out, plot=False)
    return dict(out=out, opts=opts, AICcs=AICcs)


def test_rolling_forecasts(rolling):
    out = rolling["out"]
    assert rolling["AICcs"][:, 0].tolist() == [1, 2, 3]
    tail = "tstart-1-tend-3-calibrationperiod-20-horizon-4-Mpox-cases.csv"

    names, params = table(out, "parameters-rollingwindow-%s-%s" % (QUICK, tail))
    assert params.shape == (3, 16) and params[:, 0].tolist() == [1, 2, 3]
    assert params[:, 13].tolist() == [14, 16, 34]  # I0 = first observation of each window
    np.testing.assert_allclose(params[0, 1], 1.7302, rtol=0.02)

    for prefix, ncol in (("MCSES-rollingwindow", 6), ("SCIS-rollingwindow", 6), ("diagnostics-rollingwindow", 37),
                         ("performance-calibration", 6), ("performance-forecasting", 6), ("AICcs-rollingwindow", 5)):  # fmt: skip
        _, t = table(out, "%s-%s-%s" % (prefix, QUICK, tail))
        assert t.shape == (3, ncol), prefix
        assert t[:, 0].tolist() == [1, 2, 3], prefix

    # one MAT file per window, named with its own start row and the last start row
    for i in (1, 2, 3):
        assert os.path.isfile(os.path.join(
            out, "Forecast-growthModel-Most_Recent_Timeseries_US-CDC.txt-%s-tstart-%d-tend-3-calibrationperiod-20-forecastingperiod-4.mat" % (QUICK, i)))  # fmt: skip

    run_quiet(gp.plotForecast_GrowthModels, 1, 3, 20, 4, rolling["opts"], output_dir=out)
    assert plt.get_fignums() == [101, 102, 103, 201, 202, 203, 300, 400]
    assert len(plt.figure(400).get_axes()) == 3  # one panel per window
    assert len(plt.figure(300).get_axes()) == 5  # the series above r, a, p and K

    for i, first in ((1, 0), (2, 1), (3, 2)):
        _, traj = table(out, "Forecast-%s-tstart-%d-calibrationperiod-20-horizon-4-Mpox-cases.csv" % (QUICK, i))
        assert traj.shape == (24, 5) and traj[0, 0] == first and traj[-1, 0] == first + 23

    # the plotting function rewrites the rolling tables from the saved results
    _, params2 = table(out, "parameters-rollingwindow-%s-%s" % (QUICK, tail))
    np.testing.assert_allclose(params2, params, rtol=1e-12)


def test_rolling_fits_and_plotRt(tmp_path):
    out = str(tmp_path)
    opts = gp.options_fit().replace(B=12, numstartpoints=6)
    AICcs, _ = run_quiet(gp.Run_Fit_GrowthModels, 5, 7, 20, opts, output_dir=out, plot=False)
    assert AICcs.shape == (3, 5)

    param_R, _, _ = run_quiet(gp.plotFit_ReproductionNumber, 5, 7, 20, opts, gp.options_Rt, output_dir=out)
    assert param_R.shape == (3, 5) and param_R[:, 0].tolist() == [23, 24, 25]
    assert plt.get_fignums() == [106, 107, 108, 300, 501, 502, 503]

    Rts, params, performanceC, AICcs2 = run_quiet(gp.plotRt, 5, 7, 20, opts, gp.options_Rt, output_dir=out)
    assert plt.get_fignums() == [601, 602, 603, 604]
    assert Rts.shape == (3, 4) and Rts[:, 0].tolist() == [23, 24, 25]
    np.testing.assert_allclose(Rts, param_R[:, :4], rtol=1e-10)
    assert params.shape == (3, 16) and performanceC.shape == (3, 6)
    np.testing.assert_allclose(AICcs2, AICcs, rtol=1e-12)
    names, rts = table(out, "Rts-flag1-1-tstart-5-tend-7-fixI0-1-method-0-dist-0-calibrationperiod-20-horizon-0-Mpox-cases.csv")
    assert names == ["time", "mean Rt", "Rt LB", "Rt UB"]
    np.testing.assert_allclose(rts, Rts, rtol=1e-12)


# ---------------------------------------------------------------- other estimation methods and options


@pytest.mark.parametrize(
    "options,method1,ncols",
    [(gp.options_forecast_pois, 1, 16), (gp.options_forecast_nb, 3, 19),
     (lambda: gp.options_forecast().replace(method1=4), 4, 19),
     (lambda: gp.options_forecast().replace(method1=5), 5, 22)],
)  # fmt: skip
def test_likelihood_methods(tmp_path, options, method1, ncols):
    out = str(tmp_path)
    opts = options().replace(B=15)
    AICcs, forecast_model1 = run_quiet(gp.Run_Forecasting_GrowthModels, 1, 1, 20, 4, opts, output_dir=out, plot=False)
    model = "flag1-1-fixI0-1-method-%d-dist-%d" % (method1, method1)

    names, params = table(out, "parameters-rollingwindow-%s-tstart-1-tend-1-calibrationperiod-20-horizon-4-Mpox-cases.csv" % model)
    assert params.shape == (1, ncols)
    row = dict(zip(names, params[0]))
    # count likelihoods weight the early small counts more than least squares does
    assert 1.0 < row["r estimate"] < 2.0 and 0.83 < row["p estimate"] < 0.92
    np.testing.assert_allclose(row["K0 estimate"], 28400, rtol=0.03)
    if method1 >= 3:
        assert row["alpha estimate"] > 0 and row["alpha LB"] <= row["alpha UB"]
    if method1 == 5:
        assert 0.6 <= row["d estimate"] <= 1e3

    W = gpio.load_mat(os.path.join(
        out, "Forecast-growthModel-Most_Recent_Timeseries_US-CDC.txt-%s-tstart-1-tend-1-calibrationperiod-20-forecastingperiod-4.mat" % model))  # fmt: skip
    draws = W.forecast_model12
    assert np.all(draws[1:, :] == np.round(draws[1:, :])) and np.all(draws >= 0)  # counts

    # AICc from the full count likelihood: k = 3 growth parameters + those of the error model
    k = {1: 3, 3: 4, 4: 4, 5: 5}[method1]
    assert AICcs[0, 4] == k and math.isfinite(AICcs[0, 1])
    assert math.isclose(AICcs[0, 3], 2 * k + 2 * k * (k + 1) / (19 - k - 1))

    run_quiet(gp.plotForecast_GrowthModels, 1, 1, 20, 4, opts, output_dir=out)
    assert plt.get_fignums() == [101, 201]


def test_likelihood_ranking_is_consistent(tmp_path):
    # Poisson is nested in the negative binomial models: their likelihood part cannot be worse.
    part1 = {}
    for method1 in (1, 3, 4):
        opts = gp.options_fit().replace(method1=method1, B=2)
        AICcs, _ = run_quiet(gp.Run_Fit_GrowthModels, 1, 1, 20, opts, output_dir=str(tmp_path), plot=False)
        part1[method1] = AICcs[0, 2]
    assert part1[3] <= part1[1] + 1e-6 and part1[4] <= part1[1] + 1e-6
    # this series is strongly over-dispersed relative to Poisson
    assert part1[1] - part1[4] > 100


def test_least_squares_with_poisson_and_negative_binomial_noise(tmp_path):
    out = str(tmp_path)
    for dist1 in (1, 2):
        opts = gp.options_fit().replace(dist1=dist1, B=10)
        run_quiet(gp.Run_Fit_GrowthModels, 1, 1, 20, opts, output_dir=out, plot=False)
        W = gpio.load_mat(os.path.join(
            out, "Forecast-growthModel-Most_Recent_Timeseries_US-CDC.txt-flag1-1-fixI0-1-method-0-dist-%d-tstart-1-tend-1-calibrationperiod-20-forecastingperiod-0.mat" % dist1))  # fmt: skip
        np.testing.assert_allclose(W.P_model1d.ravel()[0], 1.7302, rtol=0.02)  # the fit itself does not depend on dist1
        if dist1 == 2:
            # Pearson dispersion of the residuals around the fitted curve: far above 1 for this series
            assert gpio.scalar(W.factor1) > 5 and gpio.scalar(W.dist1_boot) == 2 and gpio.scalar(W.poissonFallback) == 0
        else:
            assert gpio.scalar(W.factor1) == 1 and gpio.scalar(W.dist1_boot) == 1


def test_negative_binomial_noise_falls_back_to_poisson_without_overdispersion(tmp_path):
    # Counts that follow the model curve up to rounding are not over-dispersed (Pearson
    # dispersion below 1), so dist1 = 2 uses Poisson draws for the bootstrap of that window.
    folder = tmp_path / "input"
    folder.mkdir()
    y = np.round(gp.growthIncidenceCurve([0.5, 0.9, 1, 4000, 5], np.arange(30.0), 1))
    np.savetxt(folder / "smooth.txt", np.column_stack([np.arange(30), y]), fmt="%d")
    out = str(tmp_path / "output")
    opts = gp.Options(cadfilename1="smooth", caddisease="Sim", datatype="cases", dist1=2, B=4, windowsize1=30)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gp.Run_Fit_GrowthModels(1, 1, 30, opts, input_dir=str(folder), output_dir=out, plot=False, verbose=False)
    name = [f for f in os.listdir(out) if f.startswith("Forecast-growthModel-")][0]
    assert "-dist-2-" in name  # the file name keeps dist-2
    W = gpio.load_mat(os.path.join(out, name))
    assert gpio.scalar(W.poissonFallback) == 1 and gpio.scalar(W.dist1_boot) == 1 and gpio.scalar(W.factor1) == 1
    np.testing.assert_allclose(W.P_model1d.ravel()[[0, 1, 3]], [0.5, 0.9, 4000], rtol=0.02)  # the curve is recovered


def test_estimated_initial_value(tmp_path):
    opts = gp.options_forecast().replace(fixI0=0, B=15)
    AICcs, _ = run_quiet(gp.Run_Forecasting_GrowthModels, 1, 1, 20, 4, opts, output_dir=str(tmp_path), plot=False)
    names, params = table(str(tmp_path), "parameters-rollingwindow-flag1-1-fixI0-0-method-0-dist-0-tstart-1-tend-1-calibrationperiod-20-horizon-4-Mpox-cases.csv")
    row = dict(zip(names, params[0]))
    assert AICcs[0, 4] == 5  # r, p, K, I0 and the residual variance
    assert 2.8 <= row["I0 estimate"] <= 70 and row["I0 LB"] < row["I0 UB"]  # searched in [0.2, 5] times the first count


def test_parameter_recovery_on_the_simulated_series(tmp_path):
    # options_fit_GLM: Poisson data simulated from the GLM with r = 0.18, p = 0.94, K = 10000
    opts = gp.options_fit_GLM().replace(B=25)
    run_quiet(gp.Run_Fit_GrowthModels, None, None, None, opts, output_dir=str(tmp_path), plot=False)
    name = [f for f in os.listdir(tmp_path) if f.startswith("parameters-rollingwindow")][0]
    assert name == "parameters-rollingwindow-flag1-1-fixI0-1-method-1-dist-1-tstart-1-tend-1-calibrationperiod-120-horizon-0-GLM-cases.csv"
    names, params = table(str(tmp_path), name)
    row = dict(zip(names, params[0]))
    np.testing.assert_allclose(row["r estimate"], 0.18, rtol=0.12)
    np.testing.assert_allclose(row["p estimate"], 0.94, rtol=0.03)
    np.testing.assert_allclose(row["K0 estimate"], 10000, rtol=0.03)
    assert row["K0 LB"] < 10000 * 1.02 and row["K0 UB"] > 10000 * 0.98


def test_forecast_beyond_the_data(tmp_path):
    # options_forecast_HIV: Gompertz, 12 yearly observations, 11 steps ahead of a 16-row series (getperformance = 0)
    out = str(tmp_path)
    opts = gp.options_forecast_HIV().replace(B=15)
    AICcs, forecast_model1 = run_quiet(gp.Run_Forecasting_GrowthModels, None, None, None, None, opts, output_dir=out, plot=False)
    assert forecast_model1.shape == (23, 15)
    assert not any(f.startswith("performance-forecasting") for f in os.listdir(out))

    run_quiet(gp.plotForecast_GrowthModels, None, None, None, None, opts, output_dir=out)
    assert plt.get_fignums() == [101]  # no performance figure
    names, traj = table(out, "Forecast-flag1-5-fixI0-0-method-0-dist-0-tstart-1-calibrationperiod-12-horizon-11-HIV-cases.csv")
    assert traj.shape == (23, 5)
    assert np.all(np.isfinite(traj[:16, :2])) and np.all(np.isnan(traj[16:, :2]))  # no time label or count beyond the data
    assert np.all(np.isfinite(traj[:, 2:]))


def test_other_models_run_end_to_end(tmp_path):
    for flag1 in (gp.EXP, gp.GGM, gp.GRM, gp.LM, gp.RICH, gp.GOM):
        opts = gp.options_forecast().replace(flag1=flag1, B=4, numstartpoints=4)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            AICcs, fm = gp.Run_Forecasting_GrowthModels(1, 1, 20, 2, opts, output_dir=str(tmp_path), plot=False, verbose=False)
        assert fm.shape == (22, 4) and np.all(np.isfinite(fm)), flag1
        assert AICcs[0, 4] == gp.get_nparams(0, 0, flag1, 1)


# ---------------------------------------------------------------- messages


def test_messages_for_unusable_requests(tmp_path):
    out = str(tmp_path)

    # more estimated quantities than observations
    with pytest.raises(gp.GrowthPredictError, match="calibration period"):
        gp.Run_Forecasting_GrowthModels(1, 1, 4, 2, output_dir=out, plot=False, verbose=False)

    # the horizon cannot be scored: rows 13-32 leave nothing to hold out
    with pytest.raises(gp.GrowthPredictError, match="too short"):
        gp.Run_Forecasting_GrowthModels(13, 13, 20, 4, gp.options_forecast().replace(B=2), output_dir=out, plot=False, verbose=False)

    # first window beyond the last possible one (options_fit_mpox as shipped starts at row 37 of 32)
    with pytest.warns(gp.GrowthPredictWarning, match="incompatible tstart and tend"):
        AICcs, curves = gp.Run_Fit_GrowthModels(None, None, None, gp.options_fit_mpox, output_dir=out, plot=False, verbose=False)
    assert AICcs.shape == (0, 5)

    with pytest.raises(FileNotFoundError):
        gp.Run_Fit_GrowthModels(1, 1, 20, gp.options_fit().replace(cadfilename1="no-such-file"), output_dir=out, plot=False)

    # plotting before running
    with pytest.raises(FileNotFoundError, match="Run"):
        gp.plotForecast_GrowthModels(1, 1, 20, 4, output_dir=out, verbose=False)
    with pytest.raises(FileNotFoundError, match="Run"):
        gp.plotFit_GrowthModels(1, 1, 20, output_dir=out, verbose=False)

    with pytest.raises(gp.GrowthPredictError) as e:
        gp.Run_Fit_GrowthModels(1, 1, 20, None, "rBounds", output_dir=out, plot=False)
    assert e.value.identifier == "GrowthPredict:fitopts"

    with pytest.raises(gp.GrowthPredictError) as e:
        gp.Run_Fit_GrowthModels(1, 1, 20, gp.options_fit().replace(dist1=4, B=2), output_dir=out, plot=False, verbose=False)
    assert e.value.identifier == "GrowthPredict:errorStructure"


def test_user_interval_of_r_reaches_the_bootstrap(tmp_path):
    # README: Run_Fit_GrowthModels(1, 1, 20, [], struct('rBounds', [1e-3 50]))
    opts = gp.options_fit().replace(B=5)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        gp.Run_Fit_GrowthModels(13, 13, 20, opts, {"rBounds": [1e-3, 10]}, output_dir=str(tmp_path), plot=False, verbose=False)
    idents = {getattr(w.message, "identifier", "") for w in caught}
    assert "GrowthPredict:BestFitOnBound" in idents  # r = 10 sits on the user's upper bound
    name = [f for f in os.listdir(tmp_path) if f.startswith("Forecast-growthModel-")][0]
    W = gpio.load_mat(os.path.join(str(tmp_path), name))
    np.testing.assert_allclose(W.bootUB[:, 0], 10)
    np.testing.assert_allclose(W.bootLB[:, 0], 1e-3)
    assert np.all(W.Phatss_model1[:, 0] <= 10 * (1 + 1e-9))
    _, diag = table(str(tmp_path), [f for f in os.listdir(tmp_path) if f.startswith("diagnostics")][0])
    assert diag[0, 2 + 4 * 7] == 1  # 'r best fit on bound'


def test_first_observation_zero_is_replaced(tmp_path):
    folder = tmp_path / "input"
    folder.mkdir()
    y = np.round(gp.growthIncidenceCurve([0.6, 0.9, 1, 3000, 2], np.arange(25.0), 1))
    y[0] = 0
    np.savetxt(folder / "zero-start.txt", np.column_stack([np.arange(25), y]), fmt="%d")
    opts = gp.Options(cadfilename1="zero-start", caddisease="Sim", datatype="cases", B=3, windowsize1=25)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        AICcs, fits = gp.Run_Fit_GrowthModels(1, 1, 25, opts, input_dir=str(folder), output_dir=str(tmp_path / "output"), plot=False, verbose=False)
    assert "GrowthPredict:FirstPointZero" in {getattr(w.message, "identifier", "") for w in caught}
    assert np.all(fits[0, :] == 0.01)
    assert os.path.isfile(str(tmp_path / "output" / "parameters-rollingwindow-flag1-1-fixI0-1-method-0-dist-0-tstart-1-tend-1-calibrationperiod-25-horizon-0-Sim-cases.csv"))


# ---------------------------------------------------------------- illustrations


def test_growth_model_illustrations():
    curve, params = gp.plotGrowthModel(1, 0.18, 0.9, None, 10000, 1, 200)
    assert curve.size == 201 and params.tolist() == [0.18, 0.9, 10000]
    np.testing.assert_allclose(np.sum(curve), 10000, rtol=1e-3)
    assert gp.plotGrowthModel(4, 0.18, None, 0.6, 10000, 1, 50)[1].tolist() == [0.18, 0.6, 10000]
    with pytest.raises(gp.GrowthPredictError, match="K value is required"):
        gp.plotGrowthModel(1, 0.18, 0.9, None, None, 1, 200)
    fig = gp.plot_growth_illustration()
    assert len(fig.get_axes()) == 4
    assert [ax.get_title() for ax in fig.get_axes()] == [
        "Generalized Logistic Growth Model", "Generalized Richards Model", "Logistic Growth Model", "Richards Model"]  # fmt: skip


# ---------------------------------------------------------------- three corrections to the plotting functions
# (1) figure 400 with 16 or 24 windows, (2) rolling-parameter figure, (3) last window start beyond the series.
# The MATLAB plotting functions were corrected in the same way; see docs/validation/plot-windows.md.


def _lines(ax):
    """x and y data of the plotted lines of an axes, in drawing order."""
    return [(np.asarray(line.get_xdata(), dtype=float), np.asarray(line.get_ydata(), dtype=float)) for line in ax.get_lines()]


@pytest.mark.parametrize(
    "nwin,grid,options,windowsize1",
    [(16, (4, 4), gp.options_forecast, 10), (24, (4, 6), gp.options_forecast_cholera, 10)],
)
def test_panel_of_forecasts_with_many_windows(tmp_path, nwin, grid, options, windowsize1):
    # MATLAB v1.1.1 stopped here: with 16 or 24 windows the grid of figure 400 was 2-by-2.
    out = str(tmp_path)
    opts = options().replace(B=2, numstartpoints=3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        AICcs, _ = gp.Run_Forecasting_GrowthModels(1, nwin, windowsize1, 4, opts, output_dir=out, plot=False, verbose=False)
        gp.plotForecast_GrowthModels(1, nwin, windowsize1, 4, opts, output_dir=out, verbose=False)
    assert AICcs.shape[0] == nwin

    # more than 11 windows: the figures of earlier windows are closed, figure 400 keeps every panel
    assert plt.get_fignums() == [100 + nwin, 200 + nwin, 300, 400]
    panels = plt.figure(400).get_axes()
    assert len(panels) == nwin
    assert panels[0].get_subplotspec().get_gridspec().get_geometry() == grid
    data, _ = gp.load_timeseries(opts.cadfilename1)
    for k, ax in enumerate(panels):
        # panel k shows the window that starts at row k + 1: its x-axis runs over calibration and forecast times
        np.testing.assert_allclose(ax.get_xlim(), [data[k, 0], data[k + windowsize1 + 3, 0]])
        assert len(ax.collections) == 2 and len(ax.get_lines()) >= 5  # trajectories, fits, limits, data, window end
    assert [ax.get_xlabel() for ax in panels] == [""] * (nwin - grid[1]) + ["Time"] * grid[1]  # bottom row only
    assert panels[0].get_ylabel() == "%s %s" % (opts.caddisease, opts.datatype)

    # one trajectory table per window
    for i in (1, nwin):
        assert os.path.isfile(os.path.join(
            out, "Forecast-flag1-1-fixI0-1-method-0-dist-0-tstart-%d-calibrationperiod-%d-horizon-4-%s-%s.csv"
            % (i, windowsize1, opts.caddisease, opts.datatype)))  # fmt: skip


def test_rolling_parameter_figure_reads_one_set_of_arrays(rolling):
    # Each of the four lower panels (r, a, p, K) shows the estimate, its two limits and the smoothed
    # estimate of the SAME parameter, taken from the windows that were plotted.
    out = rolling["out"]
    run_quiet(gp.plotForecast_GrowthModels, 1, 3, 20, 4, rolling["opts"], output_dir=out)
    _, params = table(out, "parameters-rollingwindow-%s-tstart-1-tend-3-calibrationperiod-20-horizon-4-Mpox-cases.csv" % QUICK)
    columns = {"r": 1, "p": 4, "a": 7, "K": 10}

    axes = plt.figure(300).get_axes()
    assert [ax.get_ylabel().split(" ")[0] for ax in axes[1:]] == ["r", "a", "p", "K"]
    for ax in axes[1:]:
        c = columns[ax.get_ylabel().split(" ")[0]]
        lines = _lines(ax)
        assert len(lines) == 4
        for x, _ in lines:
            assert x.tolist() == [1, 2, 3]
        np.testing.assert_allclose(lines[0][1], params[:, c], rtol=1e-12)  # estimate
        np.testing.assert_allclose(lines[1][1], params[:, c + 1], rtol=1e-12)  # lower limit
        np.testing.assert_allclose(lines[2][1], params[:, c + 2], rtol=1e-12)  # upper limit, same parameter
        np.testing.assert_allclose(lines[3][1], gp.smooth(params[:, c], 5), rtol=1e-12)  # smoothed estimate

    # the same holds for the forecast plot with the reproduction number
    run_quiet(gp.plotForecast_ReproductionNumber, 1, 3, 20, 4, rolling["opts"], gp.options_Rt, output_dir=out)
    axes = plt.figure(300).get_axes()
    a_panel = [ax for ax in axes[1:] if ax.get_ylabel() == "a"][0]
    np.testing.assert_allclose(_lines(a_panel)[2][1], params[:, 9], rtol=1e-12)


def test_last_window_start_beyond_the_series_fit(tmp_path, mpox):
    # 32 rows, windows of 20: the last window that fits starts at row 13. The runner moves tend1 = 20
    # back to 13 and names its files with tend-13; the plotting functions must look for the same names.
    out = str(tmp_path)
    opts = gp.options_fit().replace(B=3, numstartpoints=4)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        AICcs, _ = gp.Run_Fit_GrowthModels(12, 20, 20, opts, output_dir=out, plot=False, verbose=False)
    assert AICcs[:, 0].tolist() == [12, 13]
    mats = sorted(f for f in os.listdir(out) if f.startswith("Forecast-growthModel-"))
    assert [("-tstart-%d-tend-13-" % i) in f for i, f in zip((12, 13), mats)] == [True, True]
    assert not any("-tend-20-" in f for f in os.listdir(out))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        param_doubling, seq, fits = gp.plotFit_GrowthModels(12, 20, 20, opts, output_dir=out, verbose=False)
        assert plt.get_fignums() == [113, 114, 300]
        assert fits.shape == (20, 3) and param_doubling.shape == (3,)
        param_R, _, _ = gp.plotFit_ReproductionNumber(12, 20, 20, opts, gp.options_Rt, output_dir=out, verbose=False)
        assert plt.get_fignums() == [113, 114, 300, 501, 502]
        assert param_R[:, 0].tolist() == [30, 31]  # last time of the windows starting at rows 12 and 13
        Rts, params, performanceC, AICcs2 = gp.plotRt(12, 20, 20, opts, gp.options_Rt, output_dir=out, verbose=False)
    assert Rts[:, 0].tolist() == [30, 31] and params[:, 0].tolist() == [12, 13] and performanceC.shape == (2, 6)
    np.testing.assert_allclose(AICcs2, AICcs, rtol=1e-12)

    # every table is named with the adjusted value, and the same result comes from asking for it directly
    assert all("-tend-13-" in f for f in os.listdir(out) if "rollingwindow" in f or f.startswith(("Rts-", "performance-")))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, _, fits13 = gp.plotFit_GrowthModels(12, 13, 20, opts, output_dir=out, verbose=False)
    np.testing.assert_allclose(fits13, fits)


def test_last_window_start_beyond_the_series_forecast(tmp_path):
    # Forecasts beyond the data (getperformance = 0): tend1 = 20 is moved back to 13.
    out = str(tmp_path)
    opts = gp.options_forecast().replace(B=3, numstartpoints=4, getperformance=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        AICcs, _ = gp.Run_Forecasting_GrowthModels(12, 20, 20, 3, opts, output_dir=out, plot=False, verbose=False)
        _, _, fc = gp.plotForecast_GrowthModels(12, 20, 20, 3, opts, output_dir=out, verbose=False)
        assert plt.get_fignums() == [112, 113, 300, 400]
        param_R, _, _ = gp.plotForecast_ReproductionNumber(12, 20, 20, 3, opts, gp.options_Rt, output_dir=out, verbose=False)
    assert AICcs[:, 0].tolist() == [12, 13] and fc.shape == (23, 3) and param_R.shape == (2, 5)
    assert len(plt.figure(300).get_axes()) == 5
    _, traj = table(out, "Forecast-flag1-1-fixI0-1-method-0-dist-0-tstart-13-calibrationperiod-20-horizon-3-Mpox-cases.csv")
    assert traj.shape == (23, 5) and np.all(np.isnan(traj[20:, :2])) and np.all(np.isfinite(traj[:20, :2]))


def test_nothing_to_plot_when_the_first_window_does_not_fit(tmp_path):
    # options_fit_mpox as shipped starts at row 37 of a 32-row series: runner and plotting functions
    # say so and return empty results instead of looking for files.
    out = str(tmp_path)
    for function, args, shapes in (
        (gp.plotFit_GrowthModels, (gp.options_fit_mpox,), [(0,), (0, 5), (0, 0)]),
        (gp.plotFit_ReproductionNumber, (gp.options_fit_mpox, gp.options_Rt), [(0, 5), (0,), (0, 5)]),
        (gp.plotRt, (gp.options_fit_mpox, gp.options_Rt), [(0, 4), (0, 0), (0, 0), (0, 5)]),
    ):
        with pytest.warns(gp.GrowthPredictWarning, match="incompatible tstart and tend"):
            result = function(None, None, None, *args, output_dir=out, verbose=False)
        assert [r.shape for r in result] == shapes
    far = gp.options_forecast().replace(tstart1=40, tend1=45, windowsize1=20, getperformance=0)
    for function, extra, shapes in (
        (gp.plotForecast_GrowthModels, (), [(0,), (0, 5), (0, 0)]),
        (gp.plotForecast_ReproductionNumber, (gp.options_Rt,), [(0, 5), (0,), (0, 5)]),
    ):
        with pytest.warns(gp.GrowthPredictWarning, match="incompatible tstart and tend"):
            result = function(None, None, None, None, far, *extra, output_dir=out, verbose=False)
        assert [r.shape for r in result] == shapes
    assert plt.get_fignums() == [] and os.listdir(out) == []


def test_histogram_of_values_equal_up_to_rounding():
    # Bootstrap values of a parameter that every refit leaves on the same bound can differ in the last
    # digits only; such a column has no range to divide into ten bins.
    from growthpredict.plotting import draw_window_figure

    rng = np.random.default_rng(0)
    n, B = 12, 6
    Phatss = np.column_stack([
        rng.normal(2, 0.1, B),                       # r
        np.full(B, 0.01),                            # p on its lower bound
        np.ones(B),                                  # a fixed
        1e9 * (1 + 1e-16 * np.arange(B)),            # K on its upper bound, equal up to rounding
    ])  # fmt: skip
    Phatss[1, 3] = np.nextafter(1e9, 2e9)
    curves = np.abs(rng.normal(50, 5, size=(n, 20 * B)))
    fig, median1, LB1, UB1 = draw_window_figure(
        1, Phatss_model1=Phatss, params=tuple(np.array([c.mean(), c.min(), c.max()]) for c in Phatss.T),
        timevect1=np.arange(n), timevect2=np.arange(n), forecast_model12=curves, fit_model1=curves[:, :B],
        timevect_all=np.arange(n), data_all=np.full(n, 50.0), forecastingperiod=0, ylabel="cases", model_name1="GLM",
    )  # fmt: skip
    hist_axes = fig.get_axes()[:4]
    assert [sum(p.get_height() for p in ax.patches) for ax in hist_axes] == [B, B, B, B]  # every value is counted
    lo, hi = hist_axes[3].get_xlim()
    assert lo < 1e9 < hi and median1.shape == (n,)
