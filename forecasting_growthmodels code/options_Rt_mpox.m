function [type_GId1,mean_GI1,var_GI1]=options_Rt


% <=======================================================================================>
% <========================== Reproduction number number parameters =======================>
% <========================================================================================>

type_GId1=1; % type of Generation interval distribution (1=Gamma, 2=Exponential, 3=Delta)

mean_GI1=3.8;  % mean of the generation interval distribution

var_GI1=2*(1//3.8)^2; % variance of the generation interval distribution

