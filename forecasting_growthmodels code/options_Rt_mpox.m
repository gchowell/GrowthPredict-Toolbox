function [type_GId1,mean_GI1,var_GI1]=options_Rt_mpox

% OPTIONS_RT_MPOX  Generation-interval options for mpox, weekly data.
%
%   Serial interval of mpox estimated from U.S. case pairs by symptom onset:
%   mean 8.5 days (95% CrI 7.3-9.9), SD 5.0 days (95% CrI 4.0-6.4).
%
%       Madewell ZJ, Charniga K, Masters NB, et al. Serial interval and
%       incubation period estimates of monkeypox virus infection in 12
%       jurisdictions, United States, May-August 2022. Emerg Infect Dis.
%       2023;29(4). doi:10.3201/eid2904.221622
%
%   These are the values of options_Rt.m. For daily data use mean_GI1=8.5 and
%   var_GI1=5.0^2.
%
%   Until this version the file held mean_GI1=3.8 with var_GI1=2*(1/3.8)^2.
%   That variance belongs to a gamma distribution with shape 2 and rate 3.8,
%   whose mean is 2/3.8, so the two values did not describe one distribution.

% <=======================================================================================>
% <========================== Reproduction number parameters =============================>
% <=======================================================================================>

type_GId1=1; % type of Generation interval distribution (1=Gamma, 2=Exponential, 3=Delta)

mean_GI1=8.5/7;  % mean of the generation interval distribution (weeks)

var_GI1=(5.0/7)^2; % variance of the generation interval distribution (weeks^2)
