function options = growthODEOptions()
%GROWTHODEOPTIONS  One numerical policy for every ODE solve in GrowthPredict.
%
%   options = growthODEOptions()
%
%   Used by growthIncidenceCurve for the fitting objective, the final fitted
%   curve and the forecast, so the curve that is scored is the curve that is
%   reported. (Ported from quantdiffODEOptions in QuantDiffForecast.)
%
%   RelTol 1e-6 / AbsTol 1e-8 replace the ode15s defaults (1e-3 / 1e-6): the
%   optimizer differentiates the objective by finite differences, and solver
%   error at the default tolerances is of the same order as those differences.
%   The cumulative curve C(t) is constrained to be nonnegative.

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'NonNegative', 1);

end
