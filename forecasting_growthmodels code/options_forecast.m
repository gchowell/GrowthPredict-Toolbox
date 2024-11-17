% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, flag1, model_name1, fixI0, getperformance, forecastingperiod, windowsize1, tstart1, tend1] = options_forecast

% <============================================================================>
% <=================== Declare Global Variables ==============================>
% <============================================================================>
% Global variable used to define the parameter estimation method.
global method1; % Parameter estimation method

% <============================================================================>
% <========================= Dataset Properties ==============================>
% <============================================================================>
% The time series data file contains the incidence curve of interest (e.g., new cases per unit of time).
% - The first column corresponds to the time index (e.g., 0, 1, 2, ...).
% - The second column contains the temporal incidence data.
% Note:
% - If the file contains cumulative incidence data, its name must start with "cumulative".

cadfilename1 = 'Most_Recent_Timeseries_US-CDC'; % Name of the time-series data file
caddisease = 'Mpox';                            % Name of the disease or subject related to the data
datatype = 'cases';                             % Type of data (e.g., cases, deaths, hospitalizations)

% <============================================================================>
% <======================= Parameter Estimation ==============================>
% <============================================================================>
% Estimation method options:
% 0 - Nonlinear least squares (LSQ)
% 1 - Maximum Likelihood Estimation (MLE) Poisson
% 3 - MLE Negative Binomial (VAR = mean + alpha*mean)
% 4 - MLE Negative Binomial (VAR = mean + alpha*mean^2)
% 5 - MLE Negative Binomial (VAR = mean + alpha*mean^d)

method1 = 0; % Default estimation method: Nonlinear least squares (LSQ)

% Error structure options based on method1:
% 0 - Normal distribution (method1 = 0)
% 1 - Poisson error structure (method1 = 0 or 1)
% 2 - Negative Binomial (VAR = factor1 * mean, empirically estimated)
% 3 - MLE Negative Binomial (VAR = mean + alpha*mean)
% 4 - MLE Negative Binomial (VAR = mean + alpha*mean^2)
% 5 - MLE Negative Binomial (VAR = mean + alpha*mean^d)

dist1 = 0; % Default error structure: Normal distribution
switch method1
    case 1
        dist1 = 1; % Poisson error structure
    case 3
        dist1 = 3; % Negative Binomial (VAR = mean + alpha*mean)
    case 4
        dist1 = 4; % Negative Binomial (VAR = mean + alpha*mean^2)
    case 5
        dist1 = 5; % Negative Binomial (VAR = mean + alpha*mean^d)
end

% Optimization settings:
numstartpoints = 10; % Number of initial guesses for global optimization (Multistart)
B = 300;             % Number of bootstrap realizations for parameter uncertainty characterization

% <============================================================================>
% <========================== Growth Model ===================================>
% <============================================================================>
% Growth model options:
% -1: Exponential growth (EXP)
%  0: Generalized Growth Model (GGM)
%  1: Logistic Model (GLM)
%  2: Generalized Richards Model (GRM)
%  3: Linear Model (LM)
%  4: Richards Model (RICH)
%  5: Gompertz Model (GOM)

EXP = -1;  GGM = 0;  GLM = 1;  GRM = 2;  LM = 3;  RICH = 4;  GOM = 5;

flag1 = GLM;         % Selected growth model: Logistic Model (GLM)
model_name1 = 'GLM'; % Name of the selected model
fixI0 = 1;           % Boolean: Fix the initial value to the first data point (true) or estimate it (false)

% <============================================================================>
% <====================== Forecasting Parameters =============================>
% <============================================================================>
% Parameters for forecasting analysis:
% - getperformance: Boolean to enable/disable forecasting performance metrics
% - forecastingperiod: Time horizon for forecasting (number of time units ahead)

getperformance = 1;    % Enable forecasting performance metrics (1 = yes, 0 = no)
forecastingperiod = 4; % Forecast horizon: Number of time units ahead

% <============================================================================>
% <======= Parameters for Rolling Window Analysis ===========================>
% <============================================================================>
% Parameters for rolling window analysis:
% - windowsize1: Size of the moving window
% - tstart1: Start time point for rolling window analysis
% - tend1: End time point for rolling window analysis

windowsize1 = 10; % Size of the rolling window (e.g., 10 time units)
tstart1 = 1;      % Start time point for rolling window analysis
tend1 = 1;        % End time point for rolling window analysis

end
