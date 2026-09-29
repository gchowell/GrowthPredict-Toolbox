# GrowthPredict

**A MATLAB toolbox for fitting and forecasting epidemic growth trajectories with quantified uncertainty.**

GrowthPredict brings phenomenological growth models, parameter estimation, parametric bootstrapping, and forecast evaluation into one workflow. It is designed for researchers, students, and instructors working with epidemic incidence curves and other time series that can be represented by growth models based on ordinary differential equations (ODEs).

Fit a selected model to a calibration window, examine parameter uncertainty, and generate short-term forecasts. Repeat the analysis across rolling windows to study how model estimates and forecast performance change as observations accumulate. The toolbox also includes reproduction-number and doubling-time routines; review the [implementation notes](#implementation-notes-and-limitations) before interpreting these diagnostics.

[Quick start](#quick-start) · [Tutorial paper](https://doi.org/10.1038/s41598-024-51852-8) · [Video tutorials](https://www.youtube.com/watch?v=op93_wUeXXA&list=PLiMOXVNNZfvYLdwNKrIdBmH5NTvGk6IG2) · [Citation](#citation-and-tutorials)

## Contents

- [Requirements and installation](#requirements-and-installation)
- [Quick start](#quick-start) and [input data](#input-data)
- [Growth models](#growth-models) and [configuration](#configuration)
- [Rolling windows and evaluation](#rolling-windows-and-evaluation)
- [Outputs](#outputs) and [epidemiological diagnostics](#epidemiological-diagnostics)
- [Implementation notes](#implementation-notes-and-limitations), [reproducibility](#reproducibility), and [troubleshooting](#troubleshooting)
- [Citation and tutorials](#citation-and-tutorials) · [Contributing](#contributing) · [License](#license)

## Requirements and installation

### MATLAB dependencies

The current implementation uses MATLAB together with the following products:

| Product | Used for |
| --- | --- |
| [Optimization Toolbox](https://www.mathworks.com/help/optim/ug/fmincon.html) | Constrained parameter estimation with `fmincon`. |
| [Global Optimization Toolbox](https://www.mathworks.com/help/gads/multistart.html) | `MultiStart`, optimization-problem construction, and starting-point sets. |
| [Statistics and Machine Learning Toolbox](https://www.mathworks.com/help/stats/nbinrnd.html) | Distribution functions and random draws for uncertainty calculations. |

A minimum supported MATLAB release and an Octave-compatible workflow have not been established for this documentation. Record your MATLAB release and installed toolbox versions with `ver` when reporting results or issues.

### Get the code

Clone the repository in a terminal:

```bash
git clone https://github.com/gchowell/GrowthPredict-Toolbox.git
```

Alternatively, use **Code → Download ZIP** on GitHub and extract the archive.

In MATLAB, set the **Current Folder** to the cloned or extracted repository root, then run:

```matlab
cd(fullfile(pwd, 'forecasting_growthmodels code'))
addpath(pwd)

if ~isfolder('input'),  mkdir('input');  end
if ~isfolder('output'), mkdir('output'); end
```

**The code-folder name contains a space:** `forecasting_growthmodels code`, not `forecasting_growthmodels_code`. Run the workflows from this folder because input and output paths are relative to the current working directory.

The repository also contains root-level copies of some functions. Use the implementations inside the code folder consistently; check MATLAB's function resolution with:

```matlab
which Run_Forecasting_GrowthModels -all
which plotForecast_GrowthModels -all
which computeQuantiles -all
```

The first result for each command should point to the intended code folder.

## Quick start

The bundled example is [`input/Most_Recent_Timeseries_US-CDC.txt`](forecasting_growthmodels%20code/input/Most_Recent_Timeseries_US-CDC.txt), used for the U.S. mpox example. The following workflow fits the first 20 observations and then evaluates a four-step forecast against subsequent observations in the file.

### 1. Configure the example

Open [`options_fit.m`](forecasting_growthmodels%20code/options_fit.m) and [`options_forecast.m`](forecasting_growthmodels%20code/options_forecast.m). Check the following assignments **inside both functions**, preserving their existing function definitions and estimator-selection logic:

```matlab
cadfilename1 = 'Most_Recent_Timeseries_US-CDC';  % No .txt extension needed
caddisease = 'Mpox';
datatype = 'cases';

method1 = 0;          % Unweighted nonlinear least squares
dist1 = 0;            % Gaussian bootstrap/predictive noise
numstartpoints = 10;  % Random starts for the original-data fit
B = 100;              % Bootstrap datasets

flag1 = 1;            % Generalized logistic growth model
model_name1 = 'GLM';
fixI0 = 1;            % Fix the initial state to the first observation
```

Also set `getperformance = 1` in `options_forecast.m` for this retrospective example.

These are demonstration settings, not universal recommendations. In particular, Gaussian draws can be negative; see the [uncertainty limitations](#implementation-notes-and-limitations). Assess the stability of uncertainty summaries before selecting `B` for a scientific analysis.

### 2. Fit and plot

```matlab
rng(1, 'twister');
Run_Fit_GrowthModels(1, 1, 20);
plotFit_GrowthModels(1, 1, 20);
```

The arguments specify the first window-start row, last window-start row, and number of calibration observations. Here, equal start and end indices request one window.

### 3. Forecast and plot

```matlab
rng(1, 'twister');
Run_Forecasting_GrowthModels(1, 1, 20, 4);
plotForecast_GrowthModels(1, 1, 20, 4);
```

The fourth argument is the number of forecast steps. Four steps correspond to four weeks for weekly observations or four days for daily observations.

Results are written to `output/`. The plotting functions load saved `.mat` results and write additional summaries; they do not replace the fitting step. **Keep the options and window arguments consistent between running and plotting.** Existing output files with matching names can be overwritten.

Calling these functions without arguments uses the window and horizon settings in the corresponding options function. Assigning similarly named variables only in the MATLAB base workspace does not override those options.

## Input data

Place each dataset in `input/` as a numeric, two-column `.txt` file with **no header**:

| Column | Meaning |
| --- | --- |
| 1 | A regularly spaced time index, preferably `0, 1, 2, …`, with one unit per observation interval. |
| 2 | Observed incidence: the count in each interval, not the cumulative total. |

For example, the beginning of the bundled file is:

```text
0   14
1   16
2   34
3   80
4   149
5   261
```

Set `cadfilename1` to the filename without its extension, for example `'my_incidence_series'`. Use descriptive names of at least ten characters to avoid the current short-filename prefix-check issue.

### Cumulative observations

A filename beginning with `cumulative`, ignoring case, triggers conversion to incidence:

```matlab
incidence = [cumulative(1); diff(cumulative)];
```

For example, use `cumulative_mpox_series.txt`. A hyphen is not required. Cumulative data without this prefix will be interpreted as incidence; incidence data with this prefix will be differenced incorrectly.

The first cumulative value is retained as the first incidence value. Check that this is appropriate for your series, especially when a file begins after an outbreak has already accumulated cases. Preparing a documented incidence series before loading is often clearer.

### Before fitting

Check that observations and times are finite, rows are ordered, and the time index has unit spacing. Missing intervals, irregular sampling, and negative reporting corrections require explicit preprocessing; do not silently replace missing observations with zeros or discard revisions. Keep the internal `DT=1` convention and express the observation interval through the time unit used in the analysis.

Poisson and negative-binomial likelihoods require **nonnegative integer counts**. Rates, normalized values, and smoothed fractional observations are not count data. Zero counts within a series are legitimate, but all-zero windows and windows beginning at zero need special attention because the current initial-state and negative-binomial boundary handling is incomplete.

## Growth models

The ODE evolves a cumulative state, $C(t)$. The fitting objective constructs the incidence vector as `[C(t_1); diff(C(t))]`; its first entry is an initialization convention rather than a difference across a preceding observed interval.

Set both `flag1` and `model_name1` to the matching values below.

| `flag1` | `model_name1` | Model | Implemented equation |
| ---: | --- | --- | --- |
| `-1` | `'EXP'` | Exponential growth | $dC/dt = rC$ |
| `0` | `'GGM'` | Generalized growth | $dC/dt = rC^p$ |
| `1` | `'GLM'` | Generalized logistic growth | $dC/dt = rC^p(1-C/K)$ |
| `2` | `'GRM'` | Generalized Richards | $dC/dt = rC^p[1-(C/K)^a]$ |
| `3` | `'LM'` | Logistic growth | $dC/dt = rC(1-C/K)$ |
| `4` | `'RICH'` | Richards | $dC/dt = rC[1-(C/K)^a]$ |
| `5` | `'GOM'` | Gompertz, time-dependent growth-rate form | $dC/dt = rC\exp(-at)$ |

Here, $r$ is a growth-scale parameter, $p$ controls departure from exponential growth, $K$ is a saturation parameter where present, and $a$ controls shape or growth-rate decay. Parameter meanings and units depend on the selected equation. The implemented Gompertz form does **not** use an independently estimated $K$.

`fixI0=1` anchors the initial ODE state to the first observation in each calibration window; `fixI0=0` estimates it within the bounds in `fit_model.m`. Exported parameter tables can contain fixed or unused entries, so interpret only the parameters active in the selected model.

Implementation: [`modifiedLogisticGrowth.m`](forecasting_growthmodels%20code/modifiedLogisticGrowth.m) and [`fit_model.m`](forecasting_growthmodels%20code/fit_model.m).

## Configuration

Use `options_fit.m` for fitting, `options_forecast.m` for forecasting, and `options_Rt.m` for generation-interval assumptions.

| Setting | Meaning |
| --- | --- |
| `cadfilename1` | Input filename in `input/`; use the base name without `.txt`. |
| `caddisease`, `datatype` | Labels used in plots and output filenames. |
| `flag1`, `model_name1`, `fixI0` | Growth model and initial-state treatment. |
| `method1`, `dist1` | Fitting criterion and bootstrap/predictive observation model; see below. |
| `numstartpoints` | Number of random starts for the original-data fit, in addition to an informed starting point. Current bootstrap refits use two random starts internally. |
| `B` | Number of simulated datasets to refit for bootstrap uncertainty. |
| `windowsize1` | Number of observations in each calibration window. |
| `tstart1`, `tend1` | First and last **window-start row indices**, inclusive; not the two endpoints of one calibration window. |
| `forecastingperiod` | Number of observations to forecast beyond each calibration window. Forecast options only. |
| `getperformance` | Controls some forecast reporting and plotting paths. It does not disable every evaluation call in the current forecasting runner. |

### Fitting criterion: `method1`

Let $\mu$ denote model-predicted incidence, $\alpha$ a dispersion parameter, and $d$ a variance-power parameter.

| `method1` | Estimation method | Observation variance for MLE |
| ---: | --- | --- |
| `0` | Unweighted nonlinear least squares (sum of squared errors) | Not selected by `dist1` in the fitting objective. |
| `1` | Poisson maximum likelihood | $\operatorname{Var}(Y)=\mu$ |
| `3` | Negative-binomial maximum likelihood | $\operatorname{Var}(Y)=\mu+\alpha\mu$ |
| `4` | Negative-binomial maximum likelihood | $\operatorname{Var}(Y)=\mu+\alpha\mu^2$ |
| `5` | Negative-binomial maximum likelihood | $\operatorname{Var}(Y)=\mu+\alpha\mu^d$ |

Use the method codes listed here; `method1=2` is not implemented in the active objective-function switch.

### Bootstrap and predictive noise: `dist1`

For `method1=0`, choose `dist1=0` for Gaussian noise, `dist1=1` for Poisson noise, or `dist1=2` for negative-binomial noise with an empirically estimated variance-to-mean factor. **Changing `dist1` does not turn least squares into weighted least squares.**

For MLE, the options functions set `dist1` automatically to match `method1`: `1`, `3`, `4`, or `5`. Retain that selection logic when editing the options.

The bootstrap refits simulated calibration datasets. Parameter-driven trajectories represent variability across those refits; additional observation-level draws are used for predictive intervals and quantiles. These are different uncertainty summaries. Increasing `B` improves Monte Carlo resolution but does not establish model adequacy or correct optimization and observation-model defects.

Implementation: [`plotModifiedLogisticGrowthMethods1.m`](forecasting_growthmodels%20code/plotModifiedLogisticGrowthMethods1.m), [`AddErrorStructure.m`](forecasting_growthmodels%20code/AddErrorStructure.m), and the two runner functions.

## Rolling windows and evaluation

For a window starting at row `i`, the calibration rows are:

```matlab
i : i + windowsize1 - 1
```

For example, the following requests five overlapping 20-observation windows and a four-step forecast from each:

```matlab
rng(1, 'twister');
Run_Forecasting_GrowthModels(1, 5, 20, 4);
plotForecast_GrowthModels(1, 5, 20, 4);
```

The first window uses rows `1:20`; the last uses rows `5:24`. For a dataset with `N` observations, choose:

```text
Fitting:                  tend1 + windowsize1 - 1 <= N
Fully scored forecasting: tend1 + windowsize1 + forecastingperiod - 1 <= N
```

These are **MATLAB row indices**, even when the first time label is zero. For live forecasting, future observations are unavailable and their forecast-error scores are undefined. Setting `getperformance=0` suppresses some reporting, but the current runner still calls evaluation helpers; check warnings and the output limitations below.

### What the scores mean

The standard performance CSVs report mean absolute error (MAE), mean squared error (MSE), empirical coverage of the nominal 95% prediction interval, and weighted interval score (WIS). Coverage is expressed as a percentage. RMSE and mean interval score are also computed internally and stored in performance `.mat` results.

Forecast metrics at horizon `h` summarize the first **`h` forecast observations together**, not just the error at the `h`-step-ahead observation. The runner-level forecast summary retains the selected maximum horizon for each window. Distinguish calibration scores from out-of-sample forecast scores.

MAE and MSE use the median of parameter-driven trajectories, whereas predictive quantiles and WIS use observation-level draws. These medians need not coincide. When comparing models, keep forecast targets, calibration windows, observation units, and score definitions aligned.

Implementation: [`computeforecastperformance.m`](forecasting_growthmodels%20code/computeforecastperformance.m), [`computeWIS.m`](forecasting_growthmodels%20code/computeWIS.m), and [`Run_Forecasting_GrowthModels.m`](forecasting_growthmodels%20code/Run_Forecasting_GrowthModels.m).

## Outputs

Outputs are saved under `output/`. The table lists the principal products, not every intermediate file. `*` represents the run-specific filename suffix.

| Files | Contents and interpretation |
| --- | --- |
| `Forecast-growthModel-*.mat` | Saved state for each calibration window, including fitted parameters, bootstrap results, and trajectories. Also used for fit-only runs with forecast horizon zero. Required by plotting functions. |
| `parameters-growthModel-*.mat`, `performanceCalibration-growthModel-*.mat`, `performanceForecasting-growthModel-*.mat` | Aggregated parameter and performance results, as applicable. |
| `QuantilesCalibration-growthModel-*.mat`, `QuantilesForecastingPerformance-growthModel-*.mat` | Calibration and forecast quantile arrays, as applicable. |
| `AICcs-rollingwindow-*.csv` | Window-start index (`time`), AICc, its objective and penalty components, and parameter count. Not a combined AIC/AICc/BIC export. |
| `parameters-rollingwindow-*.csv` | Bootstrap central estimates and 2.5th/97.5th percentile bounds. **Columns labeled `mean` currently contain medians.** |
| `MCSES-rollingwindow-*.csv` | Bootstrap standard deviation divided by $\sqrt{B}$. This is a mean-MCSE formula, not an MCSE for the reported median. |
| `SCIS-rollingwindow-*.csv` | The interval-span calculation $\log_{10}(UB/LB)$. Interpret only for appropriate positive bounds; it is not a formal identifiability test. |
| `Fit-*.csv`, `Forecast-*.csv` | `time`, `data`, `median`, `LB`, `UB`; exported by plotting functions. Forecast files contain calibration and future rows. |
| `quantile-fit-*.csv`, `quantile-forecast-*.csv` | Predictive quantiles exported by plotting functions. The 23 columns run from `Q_0.010` to `Q_0.990`, including `Q_0.025`, `Q_0.500`, and `Q_0.975`. |
| `performance-calibration-*.csv`, `performance-forecasting-*.csv` | Window-start index, calibration length or forecast horizon, MAE, MSE, 95% prediction-interval coverage, and WIS, where available. |
| `Rt-*.csv`, `doublingtimes-*.csv` | Reproduction-number or sequential doubling-time summaries from the relevant plotting routines. Subject to the diagnostic limitations below. |

The quantile CSVs do **not** include a time column. Their rows follow calibration observations and then forecast observations, where present. Align them with the matching trajectory file or saved `timevect2`; do not assume a standalone time-stamped forecast format.

Some forecast CSV paths write `NaN` for future time labels as well as unavailable observations. Recover target times from the saved forecast grid rather than treating missing labels as an absence of predictions.

Preserve filename capitalization when scripting imports. Some exports encode window starts as `time`, while trajectory files contain actual time labels. Output names do not capture every setting: changing `B`, the random seed, or some optimizer settings may reuse a filename. Archive each analysis separately, and do not mistake bundled historical outputs for results of your current run.

## Epidemiological diagnostics

Configure generation-interval assumptions in [`options_Rt.m`](forecasting_growthmodels%20code/options_Rt.m):

| Setting | Meaning |
| --- | --- |
| `type_GId1` | `1`: gamma; `2`: exponential; `3`: fixed-interval (delta) distribution. |
| `mean_GI1` | Generation-interval mean in the same time units as the input index. |
| `var_GI1` | Generation-interval variance in squared time units, used for the gamma option. |

For unit conversion, a mean of five days is `5/7` weeks; a standard deviation of eight days corresponds to variance `(8/7)^2` weeks squared. These values illustrate conversion, **not a pathogen-specific recommendation**.

After generating matching saved results, the diagnostic entry points are:

```matlab
% For the single-window examples above:
plotFit_ReproductionNumber(1, 1, 20);
plotForecast_ReproductionNumber(1, 1, 20, 4);
```

Growth-model plotting routines also calculate sequential doubling-time summaries. Reproduction numbers depend on the specified generation interval and the modeled incidence history; projected values are model-based extrapolations, not directly observed transmission measurements.

## Reproducibility

Set and record the random seed immediately before each analysis. Preserve the exact input data, preprocessing decisions, all three options files, window/horizon arguments, source revision, MATLAB/toolbox versions, warnings, and generated outputs. With Git, record the revision using `git rev-parse HEAD` from the repository.

Use a separate archived output directory or a separate working copy for each analysis. Compare results across seeds and optimizer settings, inspect boundary estimates, and check the stability of bootstrap summaries. A saved seed supports repeatability; it does not establish convergence, model validity, or equivalence across software versions.

## Troubleshooting

| Symptom | Check |
| --- | --- |
| `MultiStart`, `createOptimProblem`, `fmincon`, or a distribution function is unavailable | Confirm that the required MathWorks products are installed and licensed; inspect `ver` and `which`. |
| Input file is not found | Run from `forecasting_growthmodels code`, check `input/`, and verify `cadfilename1` and capitalization. |
| An unexpected function runs | Use `which functionName -all` and remove unintended duplicate copies from the active path. |
| Plotting cannot load a `.mat` file | Run the matching fitting/forecasting function first, and use identical model, estimator, initial-state, window, and horizon settings. Older bundled results may use different filenames. |
| A window is skipped or scores cannot be computed | Check row-index bounds, calibration length, number of free parameters, and availability of all requested future observations. |
| A fit fails, saturates at a bound, or gives implausible intervals | Inspect the input scale, zeros, model-specific bounds, observation model, and solver diagnostics. Do not treat a larger bootstrap size as a repair. |

## Citation and tutorials

Please cite the toolbox paper when using GrowthPredict in research or teaching, and identify the software revision used:

Chowell, G., Bleichrodt, A., Dahal, S., et al. (2024). **GrowthPredict: A toolbox and tutorial-based primer for fitting and forecasting growth trajectories using phenomenological growth models.** *Scientific Reports*, **14**, 1630. [doi:10.1038/s41598-024-51852-8](https://doi.org/10.1038/s41598-024-51852-8).

The [video tutorial series](https://www.youtube.com/watch?v=op93_wUeXXA&list=PLiMOXVNNZfvYLdwNKrIdBmH5NTvGk6IG2) provides additional demonstrations. Tutorial materials and bundled historical outputs may reflect earlier code revisions.

Related methodological references:

Chowell, G. (2017). Fitting dynamic models to epidemic outbreaks with quantified uncertainty: A primer for parameter uncertainty, identifiability, and forecasts. *Infectious Disease Modelling*, **2**(3), 379–398. [doi:10.1016/j.idm.2017.08.001](https://doi.org/10.1016/j.idm.2017.08.001).

Bürger, R., Chowell, G., & Lara-Díaz, L. Y. (2019). Comparative analysis of phenomenological growth models applied to epidemic outbreaks. *Mathematical Biosciences and Engineering*, **16**(5), 4250–4273. [doi:10.3934/mbe.2019212](https://doi.org/10.3934/mbe.2019212).

## Contributing

Report reproducible problems through [GitHub Issues](https://github.com/gchowell/GrowthPredict-Toolbox/issues). Include the code revision, MATLAB/toolbox versions, relevant options, exact commands, error or warning text, and a small shareable dataset or synthetic example. Do not upload sensitive or restricted data. For code changes, describe the expected behavior and include a regression test demonstrating it.

## License

The project declares the **GNU General Public License v3.0 (GPL-3.0)**. A standalone `LICENSE` file containing the full license text still needs to be included in the repository.

