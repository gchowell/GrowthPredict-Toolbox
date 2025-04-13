# GrowthPredict
**GrowthPredict** is a **user-friendly MATLAB toolbox** for fitting and forecasting time-series trajectories using **phenomenological dynamic growth models based on ordinary differential equations (ODEs)**. It is especially useful for modeling **epidemic outbreaks and other processes governed by growth dynamics**.
 
<p> <a href="https://www.nature.com/articles/s41598-024-51852-8" target="_blank">GrowthPredict Tutorial</a></p> 
<p><a href="https://www.youtube.com/watch?v=op93_wUeXXA&list=PLiMOXVNNZfvYLdwNKrIdBmH5NTvGk6IG2&index=4&t=18s" target="_blank">Video Tutorial</a></p>

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


# Fitting the model to your data

To use the toolbox to fit a model to your data, you just need to:

<ul>
    <li>download the code </li>
    <li>create 'input' folder in your working directory where your data is located </li>
    <li>create 'output' folder in your working directory where the output files will be stored</li>   
    <li>open a MATLAB session </li>
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

## 📚 Tutorial & Examples

For detailed examples and step-by-step tutorials:

- 📄 **Tutorial Paper**:  
  [GrowthPredict: A toolbox and tutorial-based primer for fitting and forecasting growth trajectories](https://www.nature.com/articles/s41598-024-51852-8)

- 🎥 **Video Tutorial**:  
  [YouTube Series on GrowthPredict Toolbox](https://www.youtube.com/watch?v=op93_wUeXXA&list=PLiMOXVNNZfvYLdwNKrIdBmH5NTvGk6IG2)

The tutorial showcases real-world applications using, for example, the 2022 U.S. monkeypox (mpox) epidemic dataset.

## 📄 License

This toolbox is licensed under the **GNU General Public License v3.0**.  
See the [LICENSE](LICENSE) file for more details.

# Publications

<ul>

<li>Chowell, G., Bleichrodt, A., Dahal, S. et al. GrowthPredict: A toolbox and tutorial-based primer for fitting and forecasting growth trajectories using phenomenological growth models. Sci Rep 14, 1630 (2024). https://doi.org/10.1038/s41598-024-51852-8 </li>
    
<li> Chowell, G. (2017). Fitting dynamic models to epidemic outbreaks with quantified uncertainty: A primer for parameter uncertainty, identifiability, and forecasts. Infectious Disease Modelling, 2(3), 379-398. </li>

<li> Bürger, R., Chowell, G., & Lara-Díıaz, L. Y. (2019). Comparative analysis of phenomenological growth models applied to epidemic outbreaks. Mathematical Biosciences and Engineering, 16(5), 4250-4273. </li>

</ul>
