function [incidence, C] = growthIncidenceCurve(z, tvec, flag1)
%GROWTHINCIDENCECURVE  Solve the growth model and return the incidence curve.
%
%   [incidence, C] = growthIncidenceCurve(z, tvec, flag1)
%
%   z      parameter vector; z(1:5) = [r p a K I0] (longer vectors are accepted,
%          the error-model entries alpha and d are ignored)
%   tvec   time points (at least 3, increasing)
%   flag1  growth-model code (-1 EXP, 0 GGM, 1 GLM, 2 GRM, 3 LM, 4 RICH, 5 GOM)
%
%   incidence  abs([C(1); diff(C)]), a column the length of tvec
%   C          the cumulative curve
%
%   This is the single place where the model is integrated: the objective
%   (growthObjective), the final fit and the forecast (fit_model) all call it,
%   with the tolerances of growthODEOptions. The function throws when the solve
%   is incomplete or returns non-finite / complex values; growthObjective turns
%   that into Inf for the optimizer, and fit_model lets it surface for the
%   final curve, so a failed solve can never be saved as a fitted curve.

tvec = tvec(:);
if numel(tvec) < 3
    % With a two-element time span ode15s returns every internal step instead
    % of the two requested points.
    error('GrowthPredict:growthIncidenceCurve:TooFewTimes', ...
        'At least 3 time points are needed (got %d).', numel(tvec));
end

r = z(1); p = z(2); a = z(3); K = z(4); I0 = z(5);

[tout, C] = ode15s(@modifiedLogisticGrowth, tvec, I0, growthODEOptions(), r, p, a, K, flag1);

if numel(tout) ~= numel(tvec) || size(C, 1) ~= numel(tvec) || ...
        ~isreal(C) || any(~isfinite(C(:)))
    error('GrowthPredict:growthIncidenceCurve:InvalidTrajectory', ...
        'The ODE solve must return complete finite real values at the requested times.');
end

C = C(:, 1);
incidence = abs([C(1); diff(C)]);

end
