"""GrowthPredict: fitting and forecasting growth trajectories with phenomenological growth models.

Python port of the GrowthPredict MATLAB toolbox (https://github.com/gchowell/GrowthPredict-Toolbox).
Functions keep their MATLAB names and argument order.

Quick start::

    import growthpredict as gp

    gp.rng(1)
    gp.Run_Forecasting_GrowthModels(1, 1, 20, 4)
    gp.plotForecast_GrowthModels(1, 1, 20, 4)

Reference: Chowell G, Bleichrodt A, Dahal S, et al. GrowthPredict: A toolbox and
tutorial-based primer for fitting and forecasting growth trajectories using
phenomenological growth models. Sci Rep 14, 1630 (2024).
https://doi.org/10.1038/s41598-024-51852-8
"""

from .aicc import get_nparams, getAICc, growthAICcCountAdjustment
from .doubling import getDoublingTimeCurve
from .errors import GrowthPredictError, GrowthPredictWarning
from .fitting import FitResult, fit_model, initialParams
from .make_readme_figures import make_readme_figures
from .metrics import boundHitFractions, computeforecastperformance, computeQuantiles, computeWIS
from .models import (
    EXP,
    GGM,
    GLM,
    GOM,
    GRM,
    LM,
    MODEL_NAMES,
    RICH,
    growthIncidenceCurve,
    growthODEOptions,
    modifiedLogisticGrowth,
)
from .noise import AddErrorStructure, estimatePearsonDispersion, getMeanVarianceRatio
from .objective import growthObjective, growthObjectiveValue, plotModifiedLogisticGrowthMethods1
from .options import (
    Options,
    options_fit,
    options_fit_cholera,
    options_fit_cholera_malawi,
    options_fit_EXP_cholera,
    options_fit_fluAlaska,
    options_fit_GLM,
    options_fit_GRM,
    options_fit_mpox,
    options_fit_mpoxCDC,
    options_fit_mpoxChile,
    options_fit_stochasticSIR,
    options_forecast,
    options_forecast_cholera,
    options_forecast_HIV,
    options_forecast_mpox,
    options_forecast_mpoxCDC,
    options_forecast_mpoxChile,
    options_forecast_nb,
    options_forecast_pois,
)
from .plotting import (
    plot_growth_illustration,
    plotFit_GrowthModels,
    plotFit_ReproductionNumber,
    plotForecast_GrowthModels,
    plotForecast_ReproductionNumber,
    plotGrowthModel,
    plotRt,
)
from .rt import get_Rt, options_Rt, options_Rt_cholera, options_Rt_mpox, selectRtRows
from .runners import Run_Fit_GrowthModels, Run_Forecasting_GrowthModels
from .titles import (
    GetAdjustedFontSize,
    addEstimateLine,
    fitParameterTitles,
    parameterTitleText,
    setParameterTitle,
)
from .utils import Struct, days2weeks, example_input_dir, load_timeseries, quantile, rng, smooth

__version__ = "1.1.1"

__all__ = [name for name in dir() if not name.startswith("_")]
