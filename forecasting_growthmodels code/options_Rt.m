function [type_GId1, mean_GI1, var_GI1] = options_Rt

% OPTIONS_RT  GrowthPredict options for effective reproduction number (Rt)
%
% Overview
%   Specifies the generation-interval (GI) distribution and its parameters
%   used to compute Rt from the fitted growth model (and, optionally, to
%   project Rt over short horizons).
%
% Usage
%   [type_GId1, mean_GI1, var_GI1] = options_Rt;
%
% Returns
%   type_GId1  int      GI distribution family:
%                        1=Gamma, 2=Exponential, 3=Delta (fixed interval)
%   mean_GI1   double   Mean of the GI (use the SAME time units as the data’s time step)
%   var_GI1    double   Variance of the GI (units^2, consistent with mean_GI1)
%
% Examples (unit consistency)
%   • If the data/time step are in weeks but literature GI is in days:
%       mean 5 days  →  5/7 weeks
%       sd   8 days  →  (8/7) weeks  ⇒  variance = (8/7)^2 weeks^2
%
% Notes
%   • Ensure GI units are consistent with your data’s time index (e.g., days vs weeks).
%   • Pick GI parameters from literature appropriate to your pathogen/context.


% <=======================================================================================>
% <========================== Reproduction Number Parameters =============================>
% <=======================================================================================>
% This function defines parameters related to the generation interval distribution,
% which is essential for estimating the reproduction number (Rt).

type_GId1 = 1; % Type of generation interval distribution:
               % 1 = Gamma distribution (flexible, commonly used for infectious diseases).
               % 2 = Exponential distribution (simpler, assumes memoryless property).
               % 3 = Delta distribution (fixed generation interval, no variability).

mean_GI1 = 5 / 7; % Mean generation interval in the same units as the data time step.
                  % For weekly data, 5/7 is approximately 0.71 weeks (5 days).

var_GI1 = (8 / 7)^2; % Generation-interval variance in squared data-time units.
                     % For weekly data, (8/7)^2 is approximately 1.31 weeks^2.
