"""Tests of the helpers, the data loader and the options."""

import numpy as np
import pytest

import growthpredict as gp
from growthpredict import GrowthPredictError, Options, days2weeks, load_timeseries, smooth
from growthpredict.options import FIT_OPTIONS, FORECAST_OPTIONS, resolve_options


def test_smooth_is_the_moving_average_of_the_curve_fitting_toolbox():
    y = np.array([1.0, 4, 2, 8, 5, 7, 3])
    # ends: 1 point, then 3 points, then the full span of 5
    expected = [1, (1 + 4 + 2) / 3, (1 + 4 + 2 + 8 + 5) / 5, (4 + 2 + 8 + 5 + 7) / 5, (2 + 8 + 5 + 7 + 3) / 5,
                (5 + 7 + 3) / 3, 3]  # fmt: skip
    np.testing.assert_allclose(smooth(y, 5), expected)
    np.testing.assert_allclose(smooth(y, 4), smooth(y, 3))  # an even span is reduced by one
    np.testing.assert_allclose(smooth(y[:2], 5), y[:2])
    np.testing.assert_allclose(smooth(y, 1), y)


def test_days2weeks():
    assert days2weeks(np.ones(14)).tolist() == [7, 7]
    assert days2weeks(np.arange(1, 11)).tolist() == [28, 27]  # the last week may be short


def test_rng_repeats_a_stream():
    a = gp.rng(3).random(4)
    b = gp.rng(3, "twister").random(4)
    assert np.array_equal(a, b)
    with pytest.raises(ValueError):
        gp.rng(1, "no-such-generator")


def test_example_datasets_are_found(mpox):
    assert gp.example_input_dir() is not None
    assert mpox.shape == (32, 2)
    assert mpox[:5, 1].tolist() == [14, 16, 34, 80, 149]
    # README: the four counts held out by the quick start
    assert mpox[20:24, 1].tolist() == [669, 546, 404, 329]
    data, name = load_timeseries("Most_Recent_Timeseries_US-CDC.txt")
    assert name == "Most_Recent_Timeseries_US-CDC.txt" and np.array_equal(data, mpox)


def test_local_input_folder_and_cumulative_files(tmp_path):
    folder = tmp_path / "input"
    folder.mkdir()
    C = np.array([3.0, 5, 10, 18, 18, 25])
    np.savetxt(folder / "Cumulative-cases.txt", np.column_stack([np.arange(6), C]))
    np.savetxt(folder / "cases.txt", np.column_stack([np.arange(6), C]))

    data, name = load_timeseries("Cumulative-cases", str(folder))
    assert data[:, 1].tolist() == [3, 2, 5, 8, 0, 7]  # [C(1); diff(C)]

    data, _ = load_timeseries("cases", str(folder))
    assert data[:, 1].tolist() == C.tolist()  # any other name is incidence as is

    with pytest.raises(FileNotFoundError):
        load_timeseries("no-such-series", str(folder))


def test_every_options_function_names_an_existing_dataset():
    assert len(FIT_OPTIONS) == 11 and len(FORECAST_OPTIONS) == 8
    for name, function in {**FIT_OPTIONS, **FORECAST_OPTIONS}.items():
        opts = function()
        assert isinstance(opts, Options), name
        assert getattr(gp, name) is function
        data, _ = load_timeseries(opts.cadfilename1)
        assert data.shape[0] >= opts.windowsize1, name
        assert opts.flag1 in gp.MODEL_NAMES and opts.model_name1 == gp.MODEL_NAMES[opts.flag1], name
        assert opts.method1 in (0, 1, 3, 4, 5), name
        if opts.method1 > 0:
            assert opts.dist1 == opts.method1, name


def test_defaults_of_the_readme_table():
    fit, forecast = gp.options_fit(), gp.options_forecast()
    assert (fit.flag1, fit.fixI0, fit.method1, fit.dist1, fit.numstartpoints) == (1, 1, 0, 0, 10)
    assert (fit.tstart1, fit.tend1, fit.windowsize1, fit.B) == (1, 1, 20, 300)
    assert (forecast.tstart1, forecast.tend1, forecast.windowsize1, forecast.forecastingperiod, forecast.B) == (1, 1, 10, 4, 100)
    assert forecast.getperformance == 1
    assert (gp.options_forecast_pois().method1, gp.options_forecast_pois().dist1) == (1, 1)
    assert (gp.options_forecast_nb().method1, gp.options_forecast_nb().dist1) == (3, 3)


def test_error_structure_follows_the_likelihood_and_replace():
    assert Options(method1=4, dist1=0).dist1 == 4
    assert Options(method1=0, dist1=2).dist1 == 2
    opts = gp.options_fit().replace(flag1=gp.GRM, method1=5)
    assert (opts.model_name1, opts.dist1) == ("GRM", 5)
    assert gp.options_fit().flag1 == gp.GLM  # the original is not changed

    assert resolve_options(None, gp.options_fit) == gp.options_fit()
    assert resolve_options(opts, gp.options_fit) is opts
    assert resolve_options(gp.options_fit_GLM, gp.options_fit).method1 == 1
    for bad in ("options_fit", 3, lambda: 3):
        with pytest.raises(GrowthPredictError) as e:
            resolve_options(bad, gp.options_fit)
        assert e.value.identifier == "GrowthPredict:options"
