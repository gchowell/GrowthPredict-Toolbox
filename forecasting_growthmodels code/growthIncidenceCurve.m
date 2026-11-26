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
%
%   Gompertz model. Its growth rate r*exp(-a*t) depends on time, unlike the other
%   models, so the time origin matters. Time is measured from tvec(1): r is the
%   growth rate at the first point of the curve, whatever its time label. Until
%   this version the time labels were used as they are, so the same observations
%   labelled from 12 instead of 0 needed an r larger by exp(12*a), which pushed r
%   onto its upper bound for windows starting later in a series.

tvec = tvec(:);
if numel(tvec) < 3
    % With a two-element time span ode15s returns every internal step instead
    % of the two requested points.
    error('GrowthPredict:growthIncidenceCurve:TooFewTimes', ...
        'At least 3 time points are needed (got %d).', numel(tvec));
end

if flag1 == 5
    tvec = tvec - tvec(1);   % elapsed time; the other models do not depend on t
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
