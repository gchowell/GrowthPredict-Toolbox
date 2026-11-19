% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, flag1, model_name1, fixI0, getperformance, forecastingperiod, windowsize1, tstart1, tend1] = options_forecast

% <============================================================================>
% <=================== Declare Global Variables ==============================>
% <============================================================================>
% Global variables used throughout the function.
global method1; % Parameter estimation method.

% <============================================================================>
% <========================= Dataset Properties ==============================>
% <============================================================================>
% Description of the time-series data file:
% - The file contains the incidence curve of interest (e.g., new cases per unit of time).
% - The first column corresponds to the time index (e.g., 0, 1, 2, ...).
% - The second column contains the temporal incidence data.
% - If the time series file contains cumulative incidence count data, 
%   its name must start with "cumulative".

cadfilename1 = 'cholera-epidemic'; % Name of the time-series data file.
caddisease = 'Cholera';           % Name of the disease or subject of interest.
datatype = 'cases';               % Nature of the data (e.g., cases, deaths, hospitalizations).

% <============================================================================>
% <======================= Parameter Estimation ==============================>
% <============================================================================>
% Parameter estimation method:
% - 0: Nonlinear least squares (LSQ)
% - 1: Maximum Likelihood Estimation (MLE) Poisson
% - 3: MLE Negative Binomial (VAR = mean + alpha*mean)
% - 4: MLE Negative Binomial (VAR = mean + alpha*mean^2)
% - 5: MLE Negative Binomial (VAR = mean + alpha*mean^d)

method1 = 0; % Default parameter estimation method: Nonlinear least squares (LSQ).

% Error structure assumptions:
% - 0: Normal distribution (for method1 = 0)
% - 1: Poisson error structure (for method1 = 0 or 1)
% - 2: Negative Binomial (empirical variance = factor1 * mean)
% - 3: MLE Negative Binomial (VAR = mean + alpha*mean)
% - 4: MLE Negative Binomial (VAR = mean + alpha*mean^2)
% - 5: MLE Negative Binomial (VAR = mean + alpha*mean^d)

dist1 = 0; % Default error structure: Normal distribution.
switch method1
    case 1
        dist1 = 1; % Poisson error structure.
    case 3
        dist1 = 3; % Negative Binomial (VAR = mean + alpha*mean).
    case 4
        dist1 = 4; % Negative Binomial (VAR = mean + alpha*mean^2).
    case 5
        dist1 = 5; % Negative Binomial (VAR = mean + alpha*mean^d).
end

% Optimization settings:
numstartpoints = 10; % Number of initial guesses for the global optimization process (Multistart).
B = 300;             % Number of bootstrap realizations used to quantify parameter uncertainty.

% <============================================================================>
% <=========================== Growth Model ==================================>
% <============================================================================>
% Growth model options:
% -1: Exponential growth (EXP)
%  0: Generalized Growth Model (GGM)
%  1: Logistic Model (GLM)
%  2: Generalized Richards Model (GRM)
%  3: Linear Model (LM)
%  4: Richards Model (RICH)
%  5: Gompertz Model (GOM)

EXP = -1; GGM = 0; GLM = 1; GRM = 2; LM = 3; RICH = 4; GOM = 5;

flag1 = GLM;          % Selected growth model: Logistic Model (GLM).
model_name1 = 'GLM';  % Name of the selected model.
fixI0 = 1;            % Fix initial value to the first data point (1 = true, 0 = false).

% <============================================================================>
% <======================= Forecasting Parameters ============================>
% <============================================================================>
% Parameters related to forecasting:
getperformance = 1;      % Calculate forecasting performance metrics (1 = true, 0 = false).
forecastingperiod = 4;   % Forecast horizon: Number of time units ahead to predict.

% <============================================================================>
% <============ Rolling Window Analysis Parameters ===========================>
% <============================================================================>
% Parameters for the rolling window analysis:
% - windowsize1: Size of the moving window.
% - tstart1: Starting time point for rolling window analysis.
% - tend1: Ending time point for rolling window analysis.

windowsize1 = 10; % Size of the rolling window (e.g., 10 days).
tstart1 = 1;      % Starting time point for rolling window analysis.
tend1 = 1;        % Ending time point for rolling window analysis.

end

