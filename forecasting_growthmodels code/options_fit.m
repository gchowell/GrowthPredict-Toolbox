% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>
function [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, flag1, model_name1, fixI0, windowsize1, tstart1, tend1] = options_fit

% <============================================================================>
% <=================== Declare Global Variables ==============================>
% <============================================================================>
% Global variables used throughout the function.
global method1; % Parameter estimation method

% <============================================================================>
% <========================= Dataset Properties ==============================>
% <============================================================================>
% The time series data file is a text file (*.txt) located in the input folder. 
% This file contains the incidence curve of interest (e.g., new cases per unit of time).
% - The first column corresponds to the time index: 0, 1, 2, ...
% - The second column contains the temporal incidence data.
% Note: If the time series file contains cumulative incidence count data, 
%       its name must start with "cumulative".

cadfilename1 = 'Most_Recent_Timeseries_US-CDC'; % Name of the data file containing the time-series data.
caddisease = 'Mpox';                            % Name of the disease or subject related to the time series.
datatype = 'cases';                             % Nature of the data (e.g., cases, deaths, hospitalizations).

% <============================================================================>
% <======================= Parameter Estimation ==============================>
% <============================================================================>
% Method used for parameter estimation:
% 0 - Nonlinear least squares (LSQ)
% 1 - Maximum Likelihood Estimation (MLE) Poisson
% 3 - MLE Negative Binomial (VAR = mean + alpha*mean)
% 4 - MLE Negative Binomial (VAR = mean + alpha*mean^2)
% 5 - MLE Negative Binomial (VAR = mean + alpha*mean^d)

method1 = 0; % Default estimation method: Nonlinear least squares (LSQ).

% Error structure assumptions:
% 0 - Normal distribution (for method1 = 0)
% 1 - Poisson error structure (for method1 = 0 or 1)
% 2 - Negative Binomial (VAR = factor1 * mean, empirically estimated)
% 3 - MLE Negative Binomial (VAR = mean + alpha*mean)
% 4 - MLE Negative Binomial (VAR = mean + alpha*mean^2)
% 5 - MLE Negative Binomial (VAR = mean + alpha*mean^d)

dist1 = 0; % Default error structure: Normal distribution.
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
numstartpoints = 10; % Number of initial guesses for global optimization (Multistart).
B = 300;             % Number of bootstrap realizations for parameter uncertainty characterization.

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

flag1 = GLM;         % Selected growth model: Logistic Model (GLM).
model_name1 = 'GLM'; % Name of the selected model.
fixI0 = 1;           % Boolean: Fix initial value to the first data point (true) or estimate it (false).

% <============================================================================>
% <=========== Parameters for Rolling Window Analysis =======================>
% <============================================================================>
% Settings for rolling window analysis:
% - windowsize1: Size of the moving window.
% - tstart1: Time point where rolling window analysis starts.
% - tend1: Time point where rolling window analysis ends.

windowsize1 = 20; % Size of the rolling window (e.g., 20 days).
tstart1 = 37;     % Start time point for rolling window analysis.
tend1 = 37;       % End time point for rolling window analysis.

end
