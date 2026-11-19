% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>
function [cadfilename1,caddisease,datatype, dist1, numstartpoints,B,flag1,model_name1,fixI0,windowsize1,tstart1,tend1]=options_fit

% ==============================================================================
% ========================== Declare Global Variables ==========================
% ==============================================================================

global method1; % Parameter estimation method

% ==============================================================================
% ========================== Dataset Properties ================================
% ==============================================================================

% Time series data file properties
% The data file is a text file (*.txt) located in the input folder. It contains
% the incidence curve of the epidemic. The first column corresponds to the 
% time index (e.g., 0, 1, 2, ...) and the second column corresponds to the 
% temporal incidence data. 
% If the data file contains cumulative incidence counts, its name must start
% with "cumulative".

cadfilename1 = 'cholera-epidemic'; % Name of the time-series data file
caddisease = 'Cholera';           % Name of the disease
datatype = 'cases';               % Nature of the data (e.g., cases, deaths, etc.)

% ==============================================================================
% ======================== Parameter Estimation Settings =======================
% ==============================================================================

% Parameter estimation method
% Available methods:
% 0: Nonlinear least squares (LSQ)
% 1: Maximum Likelihood Estimation (MLE) Poisson
% 3: MLE Negative Binomial (VAR = mean + alpha*mean)
% 4: MLE Negative Binomial (VAR = mean + alpha*mean^2)
% 5: MLE Negative Binomial (VAR = mean + alpha*mean^d)

method1 = 0; % Set the parameter estimation method

% Error structure assumptions
% Available structures:
% 0: Normal distribution (for method1 = 0)
% 1: Poisson error structure (for method1 = 0 or 1)
% 2: Negative Binomial (var = factor1*mean, estimated empirically)
% 3: MLE Negative Binomial (VAR = mean + alpha*mean)
% 4: MLE Negative Binomial (VAR = mean + alpha*mean^2)
% 5: MLE Negative Binomial (VAR = mean + alpha*mean^d)

dist1 = 0; % Default error structure for LSQ

if method1 > 0
    dist1 = method1; % Automatically match the error structure to the method
end

% Optimization settings
numstartpoints = 20; % Number of initial guesses for global optimization
B = 300;             % Number of bootstrap realizations for parameter uncertainty

% ==============================================================================
% ============================== Growth Model ==================================
% ==============================================================================

% Growth model types
EXP  = -1; % Exponential growth
GGM  = 0;  % Generalized Growth Model (GGM)
GLM  = 1;  % Logistic Model (GLM)
GRM  = 2;  % Generalized Richards Model (GRM)
LM   = 3;  % Linear Model (LM)
RICH = 4;  % Richards Model
GOM  = 5;  % Gompertz Model

flag1 = GLM;            % Set the growth model
model_name1 = 'GLM';    % Name of the selected growth model
fixI0 = true;           % Fix initial value to the first data point (true/false)

% ==============================================================================
% ==================== Rolling Window Analysis Parameters ======================
% ==============================================================================

windowsize1 = 10;  % Size of the moving window
tstart1 = 1;       % Start time point for rolling window analysis
tend1 = 80;        % End time point for rolling window analysis


