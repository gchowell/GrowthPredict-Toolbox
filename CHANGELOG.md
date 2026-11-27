# Changelog

Versions follow the tags of this repository. Entries list what changed for a user of the toolbox; "Rerun" says which earlier results are affected.

## Unreleased

Checked in Python. In MATLAB `runtests('tests')` passes 60 of 60, the three fit plotting functions were run with an adjusted last window start and `plotForecast_GrowthModels.m` with 16 windows (figure 400: 16 panels, 4-by-4). See `docs/validation/plot-windows.md`.

### Rerun

- **Nothing.** No fitted value, forecast or exported number changes.

### Fixed

- **`plotForecast_GrowthModels.m`: panel of forecasts (figure 400) with 16, 24, 32, 36, 40, 48, … windows.** The grid had fewer panels than windows whenever the number of windows has four or more prime factors, and the function stopped at the fifth window. The grid now comes from `forecastPanelGrid.m`. Layouts for other counts are unchanged (six windows: 2-by-3).
- **`plotForecast_GrowthModels.m`: figure 400 with more than 11 windows** kept the panel of the last window only, because `close all` at each window closed it too. The other figures are still closed at each window; figure 400 is kept.
- **Plotting functions after a run in which the runner moved the last window start back** (`plotFit_GrowthModels.m`, `plotForecast_GrowthModels.m`, `plotFit_ReproductionNumber.m`, `plotForecast_ReproductionNumber.m`, `plotRt.m`). When the series is too short for the requested `tend1`, the runners use the last window that fits and name their files with it; the plotting functions used the requested value and looked for files that had not been written. They now apply the same rule (`adjustLastWindowStart.m`). With no window left they print `incompatible tstart and tend values` and return empty outputs, like the runners.
- **`plotForecast_GrowthModels.m` and `plotForecast_ReproductionNumber.m`: rolling-parameter figure.** The upper limit of `a` was read from `param_as` instead of `param_as2`. The two arrays hold the same rows when the figure is drawn, so no figure changes.

### Added

- **`forecastPanelGrid.m`** gives the rows and columns of figure 400 for any number of windows.
- **`adjustLastWindowStart.m`** holds the rule of the runners for the last window start.
- **`tests/test_forecastPanelGrid.m` and `tests/test_adjustLastWindowStart.m`** (11 tests; no toolbox needed).

- **Python port** in the `python` folder: the package `growthpredict`, with the functions of the code folder under their MATLAB names, the same output file names, 155 tests and its own README. The MATLAB toolbox does not depend on it.
- 
## v1.1.1 (2026-11-26)

### Rerun

- **Nothing.** No fitting, forecasting or plotting function changed, so results are the same as with v1.1.0.

### Changed

- **Option files declare their own name.** The 17 bundled variants of `options_fit.m` and `options_forecast.m` (for example `options_fit_GLM.m` and `options_forecast_pois.m`) declared a function named `options_fit` or `options_forecast`. MATLAB goes by the filename, so a handle such as `@options_fit_GLM` returns the same values as before.
- **`options_fit_EXP_choleraSally.m` is now `options_fit_EXP_cholera.m`.** Pass `@options_fit_EXP_cholera`.
- **Twelve scripts from earlier analyses moved to `forecasting_growthmodels code/legacy`,** with a README describing each: the five `temp_Rt_slidingwindow*.m` scripts, `plotRt_cholera.m`, `Read_Rs_choleradistricts.m`, `estimateGrowthRatesCholera.m`, `plotParamsByErrorStructure_simulated1.m`, `timingStudy.m`, `get_Rt_old.m` and `options_Rt_xinyi.m`. No runner, plotting function or test calls them. To use one, call `addpath('legacy')` from the code folder.
- **`CITATION.cff`** names the tutorial paper (Chowell et al., *Scientific Reports*, 2024) as the preferred citation, so "Cite this repository" on GitHub gives the paper.

### Removed

- **The `output` folder** and its results from earlier versions. The runners and plotting functions create the folder on first use. `make_readme_figures(false)` needs a full run of `make_readme_figures` first.
- **`options_fit_EXP.m`.** For an exponential fit, see `options_fit_EXP_cholera.m`.

## v1.1.0 (2026-11-25)

### Rerun

- **Fits with GGM, GLM or GRM whose `r` ended on its upper bound.** The runners warned about these (`BestFitOnBound`), and `diagnostics-rollingwindow-*.csv` flags them. Typical cases are windows that start inside an epidemic and series with large counts and sub-exponential growth.
- **Reproduction numbers from any window whose first time label is not zero**, including every rolling window after the first and series whose time index starts at 1. These values were far too large.
- **Reproduction numbers computed with the default `options_Rt.m`,** whose generation interval changed.
- **Doubling times from any window that does not start at the first row of the series.**
- **Gompertz fits from any window whose first time label is not zero.** The fitted curve changes where `r` had reached its upper bound, and the reported `r` changes in every such window.

### Changed

- **Interval of `r` for GGM, GLM and GRM** (`fit_model.m`). `r` is searched in `[1e-4, 10*max(y)]/DT`. Before, the interval came from the early slope of log incidence, with `[1e-3, 10]` as fallback, which is only meaningful when `p = 1`. EXP, LM, RICH and GOM are unchanged.
- **Default generation interval** (`options_Rt.m`). Now the serial interval of mpox from U.S. case pairs by symptom onset, mean 8.5 days and SD 5.0 days (Madewell et al., 2023), to match the default data file. The former default, mean 5 days and SD 8 days, is the cholera serial interval and remains in `options_Rt_cholera.m`.
- **Reproduction-number plots and `Rt-*.csv`** (`plotFit_ReproductionNumber.m`, `plotForecast_ReproductionNumber.m`). Estimates are shown and exported only where the incidence history inside the window covers at least 95% of the generation-interval distribution. The CSV has one row per unit of time, labelled with the time each estimate refers to, and a new last column `GI coverage`. The y-axis of the plot follows the curve shown.
- **Doubling times** (`getDoublingTimeCurve.m` and the four plotting functions). The count that doubles is the cumulative count of the series: cases before the window are included. A window with falling incidence now reports few doublings or none.
- **Mean doubling time** in the plotting functions leaves out bootstrap trajectories without a doubling instead of returning `NaN`.
- **Time origin of the Gompertz model** (`growthIncidenceCurve.m`). In `dC/dt = r*C*exp(-a*t)`, `t` is now the time since the first observation of the window, so `r` is the growth rate at the start of the window. Before, the time labels were used as they are: the same observations with a later first label needed a larger `r`, which reached its upper bound for later windows. The Gompertz bounds are also divided by `DT`.
- **Titles of the parameter histograms** (`r`, `p`, `a`, `K`; both runners and the four plotting functions). The estimate and its 95% CI are written on two lines, in normal weight, at a font size that follows the width of the panel, for example `r = 1.73` over `95% CI (1.49, 2.08)`. Before, each title was one bold line at about 26 points, wider than its panel, and the four titles ran into each other. Numbers now have three significant digits instead of two, and no exponent below one million (`K = 28500`, before `K=2.9e+04`).

### Fixed

- **`get_Rt.m`: generation interval evaluated at the time labels.** The gamma and exponential distributions are now evaluated at elapsed time since the start of the curve.
- **`get_Rt.m`: fixed interval one grid step late.** The estimate is `I(t)/I(t-g)`; it was `I(t)/I(t-g-0.01)`.
- **`get_Rt.m`: indexing error for fixed intervals** such as 0.29 or 0.57, where `g/0.01` is not an integer in floating point. The interval is rounded to the nearest 0.01.
- **`Rt-*.csv` time labels** were 0.01 units earlier than the estimates.
- **`options_Rt_mpox.m`** contained a syntax error (`1//3.8`) and a mean and variance that did not describe one distribution. It now holds the values of `options_Rt.m`.
- **`options_Rt_cholera.m`** declared a function named `options_Rt`.
- **`getDoublingTimeCurve.m`** left its outputs unassigned for a curve of one point.
- **`plotFit_ReproductionNumber.m`** printed the axis limits to the console.
- **`Run_Fit_GrowthModels.m` and `Run_Forecasting_GrowthModels.m`** printed the title strings (`cad1` to `cad6`) and the handle of each plotted line to the console.

### Added

- **`fitopts.rBounds = [lower upper]`** in `fit_model.m` sets the interval of `r`. Both runners accept a struct of optimizer settings as an optional last argument, for example `Run_Fit_GrowthModels(1, 1, 20, [], struct('rBounds', [1e-3 50]))`.
- **Third output of `get_Rt.m`,** `giCoverage`: the fraction of the generation-interval distribution covered by the incidence history at each estimate.
- **`selectRtRows.m`** selects the estimates with enough history.
- **Fourth input of `getDoublingTimeCurve.m`,** `priorCumulative`: cases counted before the first point of the curve.
- **`setParameterTitle.m`, `parameterTitleText.m` and `fitParameterTitles.m`** write and size the titles of the parameter histograms. `fitParameterTitles(fig)` can be called after resizing a figure from a script; figures resized by hand follow on their own.
- **Regression tests** in `forecasting_growthmodels code/tests` (49 tests in seven files); run `runtests('tests')` from the code folder.
- **Validation notes** in `docs/validation`, with the Python re-implementations used to check each change.
- **This changelog.**

### Known limits

- With GRM, `r` can still reach its upper bound when `a` is close to zero, because the two trade off.
- With `r` free, `p` reaches its lower bound of 0.01 more often, mostly for GGM on windows that level off or decline, which that model cannot follow.
- Reproduction numbers just above the 95% coverage threshold can be 10% to 20% too high for heavy-tailed generation intervals when incidence is falling. The exported coverage allows a stricter cut.
- The first doublings of a series are short when the series begins after the outbreak did, because the earlier cases are not counted.
- With the Gompertz model, `a` stays on its lower bound for windows of early growth with no visible slowing. A wider interval was tested and not adopted, because it made short-term forecasts worse.

## v1.0.2 and earlier

See the commit history.
