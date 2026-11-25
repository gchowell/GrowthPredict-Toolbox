# Changelog

Versions follow the tags of this repository. Entries list what changed for a user of the toolbox; "Rerun" says which earlier results are affected.

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

### Fixed

- **`get_Rt.m`: generation interval evaluated at the time labels.** The gamma and exponential distributions are now evaluated at elapsed time since the start of the curve.
- **`get_Rt.m`: fixed interval one grid step late.** The estimate is `I(t)/I(t-g)`; it was `I(t)/I(t-g-0.01)`.
- **`get_Rt.m`: indexing error for fixed intervals** such as 0.29 or 0.57, where `g/0.01` is not an integer in floating point. The interval is rounded to the nearest 0.01.
- **`Rt-*.csv` time labels** were 0.01 units earlier than the estimates.
- **`options_Rt_mpox.m`** contained a syntax error (`1//3.8`) and a mean and variance that did not describe one distribution. It now holds the values of `options_Rt.m`.
- **`options_Rt_cholera.m`** declared a function named `options_Rt`.
- **`getDoublingTimeCurve.m`** left its outputs unassigned for a curve of one point.
- **`plotFit_ReproductionNumber.m`** printed the axis limits to the console.

### Added

- **`fitopts.rBounds = [lower upper]`** in `fit_model.m` sets the interval of `r`. Both runners accept a struct of optimizer settings as an optional last argument, for example `Run_Fit_GrowthModels(1, 1, 20, [], struct('rBounds', [1e-3 50]))`.
- **Third output of `get_Rt.m`,** `giCoverage`: the fraction of the generation-interval distribution covered by the incidence history at each estimate.
- **`selectRtRows.m`** selects the estimates with enough history.
- **Fourth input of `getDoublingTimeCurve.m`,** `priorCumulative`: cases counted before the first point of the curve.
- **Regression tests** in `forecasting_growthmodels code/tests` (39 tests in six files); run `runtests('tests')` from the code folder.
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
