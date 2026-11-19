function [type_GId1, mean_GI1, var_GI1] = options_Rt

% <=======================================================================================>
% <========================== Reproduction Number Parameters =============================>
% <=======================================================================================>
% This function defines parameters related to the generation interval distribution,
% which is essential for estimating the reproduction number (Rt).

type_GId1 = 1; % Type of generation interval distribution:
               % 1 = Gamma distribution (flexible, commonly used for infectious diseases).
               % 2 = Exponential distribution (simpler, assumes memoryless property).
               % 3 = Delta distribution (fixed generation interval, no variability).

mean_GI1 = 5 / 7; % Mean of the generation interval distribution (in time units, e.g., days).
                  % Example: 5/7 corresponds to approximately 0.71 days.

var_GI1 = (8 / 7)^2; % Variance of the generation interval distribution (in time units squared).
                     % Example: (8/7)^2 corresponds to approximately 1.31 days^2.
