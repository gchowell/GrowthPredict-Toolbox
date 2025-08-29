# GrowthPredict
**GrowthPredict** is a **user-friendly MATLAB toolbox** for fitting and forecasting time-series trajectories using **phenomenological dynamic growth models based on ordinary differential equations (ODEs)**. It is especially useful for modeling **epidemic outbreaks and other processes governed by growth dynamics**.
 
## 📚 Tutorial & Examples

For detailed examples and step-by-step tutorials:

- 📄 **Tutorial Paper**:  
  [GrowthPredict: A toolbox and tutorial-based primer for fitting and forecasting growth trajectories](https://www.nature.com/articles/s41598-024-51852-8)

- 🎥 **Video Tutorial**:  
  [YouTube Series on GrowthPredict Toolbox](https://www.youtube.com/watch?v=op93_wUeXXA&list=PLiMOXVNNZfvYLdwNKrIdBmH5NTvGk6IG2)

The tutorial showcases real-world applications using, for example, the 2022 U.S. monkeypox (mpox) epidemic dataset.

## ✨ Features

- **Fits a suite of phenomenological models**, including:
  - Exponential Growth
  - Generalized Growth Model (GGM)
  - Gompertz Model
  - Generalized Logistic Growth Model (GLM)
  - Richards Model
- **Flexible estimation methods**:
  - Nonlinear Least Squares (LSQ)
  - Maximum Likelihood Estimation (MLE) with support for **normal, Poisson, and negative binomial errors**
- **Forecasting with uncertainty quantification** via **parametric bootstrapping**
- **Epidemiological metrics calculation**, such as:
  - Doubling time
  - **Effective reproduction number (Rₜ)**
- **Forecast evaluation tools**:
  - Weighted Interval Score (WIS)
  - Mean Absolute Error (MAE)
  - Mean Squared Error (MSE)
- **Rolling window analysis** for time-varying parameter estimation and forecasting
- **High-quality visualizations** for:
  - Model fit and forecast curves
  - Prediction intervals
  - Forecast evaluation metrics

## 📦 Installation

```bash
git clone https://github.com/gchowell/GrowthPredict-Toolbox.git
cd GrowthPredict-Toolbox/forecasting_growthmodels_code
```

In MATLAB, add the directory to your path:

```matlab
addpath(genpath(pwd))
```

Make sure you have two folders in your working directory:

- `input` — **Store your data files here**
- `output` — **Toolbox results will be saved here**

## 📊 Usage

1. **Prepare input data**:  
   Create a `.txt` file with two columns:
   - Column 1: Time index (e.g., 0, 1, 2, ...)
   - Column 2: Observed incidence (e.g., case counts)
   - If using **cumulative incidence**, filename must start with `cumulative`.

2. **Configure settings**:  
   Modify the following files to specify model and estimation options:
   - `options_fit.m`
   - `options_forecast.m`


## Options files configurations
Quick reference for configuring GrowthPredict runs — edit the settings in options_fit.m, options_forecast.m, and options_Rt.m; these values are consumed by the toolbox’s fitting, forecasting, and Rt scripts.

| Setting                            | Where                                 | What it controls                                                                                                                                                                                                         | Typical values                                                                                                                                  |
| ---------------------------------- | ------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `cadfilename1`                     | `options_fit.m`, `options_forecast.m` | **Base name** of the input time series in `./input` (expects `<cadfilename1>.txt`, two columns, **no header**). If the series is cumulative, the filename must start with `cumulative-`.                                 | e.g., `Most_Recent_Timeseries_US-CDC`                                                                                                           |
| `caddisease`                       | `options_fit.m`, `options_forecast.m` | Disease/subject label used in outputs/filenames                                                                                                                                                                          | e.g., `Mpox`                                                                                                                                    |
| `datatype`                         | `options_fit.m`, `options_forecast.m` | Data type tag                                                                                                                                                                                                            | `cases`, `deaths`, `hospitalizations`, …                                                                                                        |
| `method1` / `dist1`                | **global** in options files / both    | Estimator & observation/error model. `method1=0` (LSQ) uses `dist1∈{0,1,2}` as weighting (Normal / Poisson-like / NegBin-like). `method1=1,3,4,5` auto-sets `dist1` to 1,3,4,5 respectively (Poisson / NegBin variants). | `method1`: `0`=LSQ, `1`=MLE Poisson, `3/4/5`=MLE NegBin; `dist1`: `0`=Normal, `1`=Poisson, `2`=NegBin (LSQ-like), `3/4/5`=NegBin (MLE variants) |
| `numstartpoints`                   | `options_fit.m`, `options_forecast.m` | MultiStart initial points for global optimization                                                                                                                                                                        | e.g., `10`                                                                                                                                      |
| `B`                                | `options_fit.m`, `options_forecast.m` | Bootstrap replicates for parameter/forecast uncertainty                                                                                                                                                                  | e.g., `100`                                                                                                                                     |
| `flag1` / `model_name1`            | `options_fit.m`, `options_forecast.m` | Growth model choice and its name (must correspond)                                                                                                                                                                       | `flag1`: `-1`=EXP, `0`=GGM, `1`=GLM, `2`=GRM, `3`=LM, `4`=RICH, `5`=GOM; `model_name1`: e.g., `GLM`                                             |
| `fixI0`                            | `options_fit.m`, `options_forecast.m` | Fix initial observed value to first datum (`1`) or estimate it (`0`)                                                                                                                                                     | `0` or `1` (often `1`)                                                                                                                          |
| `windowsize1`                      | `options_fit.m`, `options_forecast.m` | Rolling-window length (time steps)                                                                                                                                                                                       | Fit: e.g., `20`; Forecast: e.g., `10`                                                                                                           |
| `tstart1`, `tend1`                 | `options_fit.m`, `options_forecast.m` | Start/end indices of the first rolling window                                                                                                                                                                            | e.g., `1`, `1`                                                                                                                                  |
| `getperformance`                   | `options_forecast.m`                  | Compute forecast performance metrics                                                                                                                                                                                     | `0` or `1` (often `1`)                                                                                                                          |
| `forecastingperiod`                | `options_forecast.m`                  | Forecast horizon (steps ahead)                                                                                                                                                                                           | e.g., `4`                                                                                                                                       |
| `type_GId1`, `mean_GI1`, `var_GI1` | `options_Rt.m`                        | Generation-interval (GI) family & parameters for Rt (**same time units as your data/time step**)                                                                                                                         | `type_GId1`: `1`=Gamma, `2`=Exponential, `3`=Delta; e.g., `mean_GI1 = 5/7`, `var_GI1 = (8/7)^2`                                                 |


# Fitting the model to your data

To use the toolbox to fit a model to your data, you just need to:

<ul>
    <li>define the model parameter values and time series parameters by editing <code>options_fit.m</code> </li>
    <li>run the function <code>Run_Fit_GrowthModels.m</code> </li>
</ul>
  
# Plotting the best model fit and parameter estimates

After fitting the model to your data using the function <code>Run_Fit_GrowthModels.m</code>, you can use the toolbox to plot the best model fit and parameter estimates as follows:

<ul>
    <li>run the function <code>plotFit_GrowthModels.m</code> </li>
</ul>

The function also outputs files with parameter estimates, the best fit of the model, and the performance metrics for the calibration period.

# Generating forecasts

To use the toolbox to fit a model to your data and generate a forecast, you just need to:

<ul>
    <li>define the model parameter values and time series parameters by editing <code>options_forecast.m</code> </li>
    <li>run the function <code>Run_Forecasting_GrowthModels.m</code> </li>
</ul>
  
# Plotting forecasts and performance metrics based on the best-fit model

After running <code>Run_Forecasting_GrowthModels.m</code>, you can use the toolbox to plot forecasts and performance metrics derived from the best fit model as follows:

<ul>
    <li>run the function <code>plotForecast_GrowthModels.m</code></li>
</ul>

The function also outputs files with parameter estimates, the fit and forecast of the model, and the performance metrics for the calibration period and forecasting periods.


## 📄 License

This toolbox is licensed under the **GNU General Public License v3.0**.  
See the [LICENSE](LICENSE) file for more details.

# Publications

<ul>

<li>Chowell, G., Bleichrodt, A., Dahal, S. et al. GrowthPredict: A toolbox and tutorial-based primer for fitting and forecasting growth trajectories using phenomenological growth models. Sci Rep 14, 1630 (2024). https://doi.org/10.1038/s41598-024-51852-8 </li>
    
<li> Chowell, G. (2017). Fitting dynamic models to epidemic outbreaks with quantified uncertainty: A primer for parameter uncertainty, identifiability, and forecasts. Infectious Disease Modelling, 2(3), 379-398. </li>

<li> Bürger, R., Chowell, G., & Lara-Díıaz, L. Y. (2019). Comparative analysis of phenomenological growth models applied to epidemic outbreaks. Mathematical Biosciences and Engineering, 16(5), 4250-4273. </li>

</ul>
