function [tend1,adjusted]=adjustLastWindowStart(tend1,windowsize1,nObservations)

% ADJUSTLASTWINDOWSTART  Last window start that fits in the series.
%
%   [tend1,adjusted] = adjustLastWindowStart(tend1,windowsize1,nObservations)
%
% Inputs
%   tend1          row where the last rolling window starts, as requested
%   windowsize1    number of observations in a window
%   nObservations  number of rows of the series
%
% Outputs
%   tend1     the requested value, or nObservations-windowsize1+1 when a window
%             starting at the requested row would run past the end of the series
%   adjusted  true when the value was changed
%
% Why
%   Run_Fit_GrowthModels.m and Run_Forecasting_GrowthModels.m apply this rule
%   before their rolling-window loop and name every file they save with the
%   adjusted value ('...-tend-<tend1>-...'). The plotting functions call this
%   function so that they look for the same file names. Until this version they
%   used the requested value: after Run_Fit_GrowthModels(1,20,20) on a series
%   of 32 rows, which saves its files with tend-13, plotFit_GrowthModels(1,20,20)
%   looked for files with tend-20 and found none.
%
%   The result can be below the first window start (for example when the
%   window is longer than the series). The caller then has no window to
%   process, as in the runners.

if ~isnumeric(tend1) || ~isscalar(tend1) || ~isnumeric(windowsize1) || ~isscalar(windowsize1) || ...
        ~isnumeric(nObservations) || ~isscalar(nObservations)
    error('GrowthPredict:adjustLastWindowStart:input', ...
        'tend1, windowsize1 and nObservations must be numbers.');
end

adjusted=false;

if (tend1+windowsize1-1) > nObservations

    tend1=nObservations-windowsize1+1;

    adjusted=true;

end
