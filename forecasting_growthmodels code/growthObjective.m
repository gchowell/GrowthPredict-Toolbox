function [objfunction, fitcurve, C] = growthObjective(z, tvec, ydata, flag1, method1)
%GROWTHOBJECTIVE  Objective of the growth-model fit, without global variables.
%
%   objfunction = growthObjective(z, tvec, ydata, flag1, method1)
%   [objfunction, fitcurve, C] = growthObjective(...)
%
%   z = [r p a K I0 alpha d].
%
%   Called with one output (the optimizer interface), a trial whose ODE solve
%   fails or whose score is not a finite real number returns Inf. fit_model
%   maps Inf to a value above every valid objective value, so a failed solve
%   can never look better than a valid point. (The former fixed penalty, 1e10,
%   can be smaller than legitimate sums of squares.)
%
%   Called with more outputs (the final evaluation in fit_model), the same
%   failures throw, so a penalty value can never be saved as an accepted fit.
%   The returned fitcurve is exactly the curve that was scored.
%
%   Because every input is an argument, the function can run on parallel
%   workers (MultiStart 'UseParallel').

objfunction = Inf;
fitcurve = [];
C = [];

try
    [fitcurve, C] = growthIncidenceCurve(z, tvec, flag1);
    objfunction = growthObjectiveValue(fitcurve, ydata, method1, z(6), z(7));
catch ME
    if nargout > 1
        rethrow(ME);
    end
    objfunction = Inf;
end

end
