function [type_GId1, mean_GI1, var_GI1] = options_Rt_cholera

% OPTIONS_RT_CHOLERA  Generation-interval options for cholera, weekly data.
%
%   Serial interval commonly assumed for cholera: gamma distribution with mean
%   5 days and SD 8 days. For example, the analysis of the 2016-18 epidemic in
%   Yemen used a gamma distribution shifted by 1 day with mean 5 days and SD
%   8 days (appendix of the reference below).
%
%       Camacho A, Bouhenia M, Alyusfi R, et al. Cholera epidemic in Yemen,
%       2016-18: an analysis of surveillance data. Lancet Glob Health.
%       2018;6(6):e680-e690. doi:10.1016/S2214-109X(18)30230-4
%
%   get_Rt.m has no shift, so the distribution is applied from lag zero. With
%   mean 5 days and SD 8 days its shape is 0.39 and about 41% of its mass lies
%   within one day, which pulls Rt toward 1.
%
%   These were the values of options_Rt.m until its default became the mpox
%   serial interval. For daily data use mean_GI1=5 and var_GI1=8^2.

% <=======================================================================================>
% <========================== Reproduction Number Parameters =============================>
% <=======================================================================================>
% This function defines parameters related to the generation interval distribution,
% which is essential for estimating the reproduction number (Rt).

type_GId1 = 1; % Type of generation interval distribution:
               % 1 = Gamma distribution (flexible, commonly used for infectious diseases).
               % 2 = Exponential distribution (simpler, assumes memoryless property).
               % 3 = Delta distribution (fixed generation interval, no variability).

mean_GI1 = 5 / 7; % Mean of the generation interval distribution (weeks): 5 days.

var_GI1 = (8 / 7)^2; % Variance of the generation interval distribution (weeks^2): SD 8 days.
