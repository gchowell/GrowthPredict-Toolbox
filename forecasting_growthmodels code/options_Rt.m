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
%   mean_GI1   double   Mean of the GI (use the SAME time units as the data's time step)
%   var_GI1    double   Variance of the GI (units^2, consistent with mean_GI1)
%
% Default values
%   The default data file is the weekly U.S. mpox series, so the default is the
%   serial interval of mpox estimated from U.S. case pairs by symptom onset:
%   mean 8.5 days (95% CrI 7.3-9.9), SD 5.0 days (95% CrI 4.0-6.4).
%
%       Madewell ZJ, Charniga K, Masters NB, et al. Serial interval and
%       incubation period estimates of monkeypox virus infection in 12
%       jurisdictions, United States, May-August 2022. Emerg Infect Dis.
%       2023;29(4). doi:10.3201/eid2904.221622
%
%   The same study gives 7.0 days (SD 4.2) by rash onset. A longer alternative
%   is the generation time of 12.5 days (95% CI 7.5-17.3) estimated in Italy
%   (Guzzetta G, et al. Emerg Infect Dis. 2022;28(10). doi:10.3201/eid2810.221126).
%   A serial interval is used here in place of the generation interval: the two
%   have the same mean in theory but can differ in variance.
%
%   Until this version the default was mean 5 days, SD 8 days, which is the
%   serial interval commonly assumed for cholera; see options_Rt_cholera.m.
%
% Examples (unit consistency)
%   If the data/time step are in weeks but the literature GI is in days:
%       mean 8.5 days  ->  8.5/7 weeks
%       sd   5.0 days  ->  5.0/7 weeks  =>  variance = (5.0/7)^2 weeks^2
%
% Notes
%   Ensure GI units are consistent with your data's time index (e.g., days vs weeks).
%   Pick GI parameters from literature appropriate to your pathogen/context: copy
%   this file, change the values and pass the handle of the copy to the plotting
%   function.


% <=======================================================================================>
% <========================== Reproduction Number Parameters =============================>
% <=======================================================================================>
% This function defines parameters related to the generation interval distribution,
% which is essential for estimating the reproduction number (Rt).

type_GId1 = 1; % Type of generation interval distribution:
               % 1 = Gamma distribution (flexible, commonly used for infectious diseases).
               % 2 = Exponential distribution (simpler, assumes memoryless property).
               % 3 = Delta distribution (fixed generation interval, no variability).

mean_GI1 = 8.5 / 7; % Mean generation interval in the same units as the data time step.
                    % For weekly data, 8.5/7 is approximately 1.21 weeks (8.5 days).

var_GI1 = (5.0 / 7)^2; % Generation-interval variance in squared data-time units.
                       % For weekly data, (5.0/7)^2 is approximately 0.51 weeks^2 (SD 5.0 days).
