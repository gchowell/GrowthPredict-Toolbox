"""Options of the fit and forecast runners.

Python port of ``options_fit.m``, ``options_forecast.m`` and their bundled
variants.

In MATLAB an options file is a function that returns twelve (fit) or fourteen
(forecast) values and sets the global ``method1``. Here an options function
takes no input and returns an :class:`Options` object that holds the same
settings, ``method1`` included. Pass the function (or an ``Options`` object)
as the options argument of a runner or plotting function::

    Run_Forecasting_GrowthModels(1, 1, 20, 4, options_forecast_pois)

To keep a configuration of your own, write a function that returns a modified
copy::

    def options_fit_mydata():
        return options_fit().replace(cadfilename1='mydata', caddisease='Flu',
                                     flag1=GRM, model_name1='GRM', B=500)

or build one directly with ``Options(cadfilename1='mydata', ...)``.

Estimation and error models
    ``method1`` selects the estimator, ``dist1`` the observation noise used for
    the bootstrap and the predictive draws.

    method1 = 0  nonlinear least squares; choose dist1 in {0, 1, 2}:
                   0 = normal, 1 = Poisson, 2 = negative binomial with
                   VAR = factor1*mean (factor1 estimated from the residuals)
    method1 = 1  Poisson MLE                            -> dist1 = 1 (automatic)
    method1 = 3  negative binomial MLE, VAR = mean + alpha*mean     -> dist1 = 3
    method1 = 4  negative binomial MLE, VAR = mean + alpha*mean^2   -> dist1 = 4
    method1 = 5  negative binomial MLE, VAR = mean + alpha*mean^d   -> dist1 = 5

Growth models (``flag1``)
    EXP = -1, GGM = 0, GLM = 1, GRM = 2, LM = 3, RICH = 4, GOM = 5. Set
    ``flag1`` and ``model_name1`` together.

Input data (``./input``)
    Text file ``<cadfilename1>.txt`` with two columns and no header: the time
    index (0, 1, 2, ...) and the observed incidence. A file whose name begins
    with ``cumulative`` (case-insensitive) is converted to incidence on load.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

from .errors import GrowthPredictError
from .models import EXP, GGM, GLM, GOM, GRM, MODEL_NAMES

__all__ = [
    "Options",
    "resolve_options",
    "options_fit",
    "options_fit_EXP_cholera",
    "options_fit_GLM",
    "options_fit_GRM",
    "options_fit_cholera",
    "options_fit_cholera_malawi",
    "options_fit_fluAlaska",
    "options_fit_mpox",
    "options_fit_mpoxCDC",
    "options_fit_mpoxChile",
    "options_fit_stochasticSIR",
    "options_forecast",
    "options_forecast_HIV",
    "options_forecast_cholera",
    "options_forecast_mpox",
    "options_forecast_mpoxCDC",
    "options_forecast_mpoxChile",
    "options_forecast_nb",
    "options_forecast_pois",
]


@dataclass
class Options:
    """Settings of a fit or forecast analysis.

    Attributes
    ----------
    cadfilename1 : str
        Base name of the input file in ``./input`` (with or without ``.txt``).
    caddisease : str
        Disease/subject label for plots and CSV file names.
    datatype : str
        Data tag (e.g. 'cases', 'deaths', 'hospitalizations').
    method1 : int
        Estimation method (0 LSQ, 1 Poisson MLE, 3/4/5 negative binomial MLE).
    dist1 : int
        Error structure. Set automatically to ``method1`` for the likelihood
        methods (1, 3, 4, 5).
    numstartpoints : int
        Number of start points of the multi-start optimization.
    B : int
        Number of bootstrap realizations.
    flag1 : int
        Growth-model code.
    model_name1 : str
        Name of the model shown in plots (default: the abbreviation of ``flag1``).
    fixI0 : int
        1 = fix the initial value to the first data point, 0 = estimate it.
    getperformance : int
        1 = compute forecast performance metrics (forecast runner only).
    forecastingperiod : int
        Forecast horizon in time steps (forecast runner only).
    windowsize1 : int
        Length of the calibration window.
    tstart1, tend1 : int
        One-based row where the first and the last rolling window start.
    """

    cadfilename1: str = "Most_Recent_Timeseries_US-CDC"
    caddisease: str = "Mpox"
    datatype: str = "cases"
    method1: int = 0
    dist1: int = 0
    numstartpoints: int = 10
    B: int = 300
    flag1: int = GLM
    model_name1: str = None
    fixI0: int = 1
    getperformance: int = 1
    forecastingperiod: int = 0
    windowsize1: int = 20
    tstart1: int = 1
    tend1: int = 1

    def __post_init__(self):
        if self.method1 in (1, 3, 4, 5):
            self.dist1 = self.method1  # the error structure follows the likelihood
        if self.model_name1 is None:
            self.model_name1 = MODEL_NAMES.get(self.flag1, str(self.flag1))
        self.fixI0 = int(self.fixI0)

    def replace(self, **changes):
        """A copy with some settings changed (``dist1`` follows ``method1`` again)."""
        if "flag1" in changes and "model_name1" not in changes:
            changes["model_name1"] = None
        return dataclasses.replace(self, **changes)


def resolve_options(options_pass, default):
    """Turn the options argument of a runner into an :class:`Options` object.

    ``None`` -> ``default()``; a function -> its return value; an ``Options``
    object -> itself.
    """
    if options_pass is None:
        options_pass = default
    if isinstance(options_pass, Options):
        return options_pass
    if callable(options_pass):
        opts = options_pass()
        if not isinstance(opts, Options):
            raise GrowthPredictError(
                "GrowthPredict:options", "The options function must return a growthpredict.Options object."
            )
        return opts
    raise GrowthPredictError(
        "GrowthPredict:options",
        "options_pass must be an options function (e.g. options_fit) or an Options object.",
    )


# <============================================================================>
# <============================== Fit options =================================>
# <============================================================================>


def options_fit():
    """Default fit options (``options_fit.m``): U.S. mpox cases, GLM, least squares
    with normal noise, 10 start points, 300 bootstrap realizations, one window
    of 20 observations starting at row 1."""
    return Options(
        cadfilename1="Most_Recent_Timeseries_US-CDC",
        caddisease="Mpox",
        datatype="cases",
        method1=0,
        dist1=0,
        numstartpoints=10,
        B=300,
        flag1=GLM,
        model_name1="GLM",
        fixI0=1,
        windowsize1=20,
        tstart1=1,
        tend1=1,
    )


def options_fit_EXP_cholera():
    """Exponential growth fitted to the cholera epidemic (``options_fit_EXP_cholera.m``)."""
    return Options(
        cadfilename1="cholera-epidemic.txt", caddisease="Cholera", datatype="cases", method1=0, dist1=0,
        numstartpoints=10, B=300, flag1=EXP, model_name1="EXP", fixI0=1, windowsize1=30, tstart1=1, tend1=1,
    )  # fmt: skip


def options_fit_GLM():
    """Poisson MLE of the GLM on a simulated series with known parameters
    (r = 0.18, p = 0.94, K = 10000) (``options_fit_GLM.m``)."""
    return Options(
        cadfilename1="curve-GLM model-r-0.18-p-0.94-K-10000-M-1-dist1-1-factor1-1", caddisease="GLM",
        datatype="cases", method1=1, dist1=1, numstartpoints=10, B=300, flag1=GLM, model_name1="GLM",
        fixI0=1, windowsize1=120, tstart1=1, tend1=1,
    )  # fmt: skip


def options_fit_GRM():
    """Generalized Richards model fitted to NYC COVID-19 cases (``options_fit_GRM.m``)."""
    return Options(
        cadfilename1="curve-covid-nyc-cases", caddisease="Covid", datatype="cases", method1=0, dist1=0,
        numstartpoints=10, B=1000, flag1=GRM, model_name1="GRM", fixI0=1, windowsize1=100, tstart1=1, tend1=1,
    )  # fmt: skip


def options_fit_cholera():
    """Rolling GLM fits of 10 observations to the cholera epidemic (``options_fit_cholera.m``)."""
    return Options(
        cadfilename1="cholera-epidemic", caddisease="Cholera", datatype="cases", method1=0, dist1=0,
        numstartpoints=20, B=300, flag1=GLM, model_name1="GLM", fixI0=1, windowsize1=10, tstart1=1, tend1=80,
    )  # fmt: skip


def options_fit_cholera_malawi():
    """Rolling GLM fits of 7 observations to the series ``simulated`` (``options_fit_cholera_malawi.m``)."""
    return Options(
        cadfilename1="simulated", caddisease="cholera-malawi", datatype="cases", method1=0, dist1=0,
        numstartpoints=20, B=200, flag1=GLM, model_name1="GLM", fixI0=1, windowsize1=7, tstart1=1, tend1=80,
    )  # fmt: skip


def options_fit_fluAlaska():
    """Poisson MLE of the GGM on ``curve1`` (1918 influenza deaths) (``options_fit_fluAlaska.m``).

    The file's own window start (row 37) lies beyond the last possible window
    of this 50-row series; pass the window explicitly, e.g.
    ``Run_Fit_GrowthModels(31, 31, 20, options_fit_fluAlaska)``.
    """
    return Options(
        cadfilename1="curve1", caddisease="flu-1918", datatype="deaths", method1=1, dist1=1,
        numstartpoints=10, B=300, flag1=GGM, model_name1="GGM", fixI0=1, windowsize1=20, tstart1=37, tend1=37,
    )  # fmt: skip


def options_fit_mpox():
    """GLM, least squares, U.S. mpox cases (``options_fit_mpox.m``).

    Pass the window explicitly: the last 20 observations of this 32-row series
    are ``Run_Fit_GrowthModels(13, 13, 20, options_fit_mpox)``.
    """
    return Options(
        cadfilename1="Most_Recent_Timeseries_US-CDC", caddisease="Mpox", datatype="cases", method1=0, dist1=0,
        numstartpoints=10, B=300, flag1=GLM, model_name1="GLM", fixI0=1, windowsize1=20, tstart1=37, tend1=37,
    )  # fmt: skip


def options_fit_mpoxCDC():
    """GLM, least squares, U.S. mpox cases (``options_fit_mpoxCDC.m``); see :func:`options_fit_mpox`."""
    return Options(
        cadfilename1="Most_Recent_Timeseries_US-CDC", caddisease="Mpox", datatype="cases", method1=0, dist1=0,
        numstartpoints=10, B=300, flag1=GLM, model_name1="GLM", fixI0=1, windowsize1=20, tstart1=37, tend1=37,
    )  # fmt: skip


def options_fit_mpoxChile():
    """GLM, least squares, weekly mpox cases in Chile (``options_fit_mpoxChile.m``)."""
    return Options(
        cadfilename1="cases-weekly-mpox-Chile", caddisease="Mpox", datatype="cases", method1=0, dist1=0,
        numstartpoints=10, B=300, flag1=GLM, model_name1="GLM", fixI0=1, windowsize1=10, tstart1=1, tend1=1,
    )  # fmt: skip


def options_fit_stochasticSIR():
    """Poisson MLE of exponential growth on a stochastic SIR simulation (``options_fit_stochasticSIR.m``)."""
    return Options(
        cadfilename1="Stochastic SIR simulation results (n=2) Poisson R=10 TimeStep=0.1 Aggregated",
        caddisease="SEIRsim", datatype="cases", method1=1, dist1=1, numstartpoints=10, B=100, flag1=EXP,
        model_name1="EXP", fixI0=1, windowsize1=30, tstart1=1, tend1=1,
    )  # fmt: skip


# <============================================================================>
# <============================ Forecast options ==============================>
# <============================================================================>


def options_forecast():
    """Default forecast options (``options_forecast.m``): U.S. mpox cases, GLM,
    least squares with normal noise, 10 start points, 100 bootstrap
    realizations, one window of 10 observations and a 4-step forecast with
    performance metrics."""
    return Options(
        cadfilename1="Most_Recent_Timeseries_US-CDC",
        caddisease="Mpox",
        datatype="cases",
        method1=0,
        dist1=0,
        numstartpoints=10,
        B=100,
        flag1=GLM,
        model_name1="GLM",
        fixI0=1,
        getperformance=1,
        forecastingperiod=4,
        windowsize1=10,
        tstart1=1,
        tend1=1,
    )


def options_forecast_HIV():
    """Gompertz forecast of yearly HIV cases, initial value estimated (``options_forecast_HIV.m``)."""
    return Options(
        cadfilename1="yearly-HIV-cases-national-National", caddisease="HIV", datatype="cases", method1=0,
        dist1=0, numstartpoints=20, B=300, flag1=GOM, model_name1="GOM", fixI0=0, getperformance=0,
        forecastingperiod=11, windowsize1=12, tstart1=1, tend1=1,
    )  # fmt: skip


def options_forecast_cholera():
    """GLM forecast of the cholera epidemic (``options_forecast_cholera.m``)."""
    return Options(
        cadfilename1="cholera-epidemic", caddisease="Cholera", datatype="cases", method1=0, dist1=0,
        numstartpoints=10, B=300, flag1=GLM, model_name1="GLM", fixI0=1, getperformance=1,
        forecastingperiod=4, windowsize1=10, tstart1=1, tend1=1,
    )  # fmt: skip


def options_forecast_mpox():
    """GLM forecast of U.S. mpox cases beyond the data, initial value estimated (``options_forecast_mpox.m``)."""
    return Options(
        cadfilename1="Most_Recent_Timeseries_US-CDC", caddisease="monkeypox", datatype="cases", method1=0,
        dist1=0, numstartpoints=10, B=300, flag1=GLM, model_name1="GLM", fixI0=0, getperformance=0,
        forecastingperiod=7, windowsize1=20, tstart1=4, tend1=5,
    )  # fmt: skip


def options_forecast_mpoxCDC():
    """GLM forecast of U.S. mpox cases (``options_forecast_mpoxCDC.m``)."""
    return Options(
        cadfilename1="Most_Recent_Timeseries_US-CDC", caddisease="Mpox", datatype="cases", method1=0, dist1=0,
        numstartpoints=10, B=300, flag1=GLM, model_name1="GLM", fixI0=1, getperformance=1,
        forecastingperiod=4, windowsize1=10, tstart1=1, tend1=1,
    )  # fmt: skip


def options_forecast_mpoxChile():
    """Negative binomial (VAR = mean + alpha*mean) GLM forecast of mpox in Chile (``options_forecast_mpoxChile.m``)."""
    return Options(
        cadfilename1="cases-weekly-mpox-Chile", caddisease="Mpox", datatype="cases", method1=3, dist1=3,
        numstartpoints=10, B=300, flag1=GLM, model_name1="GLM", fixI0=0, getperformance=0,
        forecastingperiod=10, windowsize1=10, tstart1=1, tend1=1,
    )  # fmt: skip


def options_forecast_nb():
    """The default forecast with negative binomial MLE, VAR = mean + alpha*mean (``options_forecast_nb.m``)."""
    return options_forecast().replace(method1=3)


def options_forecast_pois():
    """The default forecast with Poisson MLE (``options_forecast_pois.m``)."""
    return options_forecast().replace(method1=1)


FIT_OPTIONS = {
    f.__name__: f
    for f in (
        options_fit, options_fit_EXP_cholera, options_fit_GLM, options_fit_GRM, options_fit_cholera,
        options_fit_cholera_malawi, options_fit_fluAlaska, options_fit_mpox, options_fit_mpoxCDC,
        options_fit_mpoxChile, options_fit_stochasticSIR,
    )
}  # fmt: skip

FORECAST_OPTIONS = {
    f.__name__: f
    for f in (
        options_forecast, options_forecast_HIV, options_forecast_cholera, options_forecast_mpox,
        options_forecast_mpoxCDC, options_forecast_mpoxChile, options_forecast_nb, options_forecast_pois,
    )
}  # fmt: skip
