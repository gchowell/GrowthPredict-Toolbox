function [adj, info] = growthAICcCountAdjustment(method1, ydata, fixI0, alpha, d)
%GROWTHAICCCOUNTADJUSTMENT  Complete the count objective to a full log-likelihood for AICc.
%
%   [adj, info] = growthAICcCountAdjustment(method1, ydata, fixI0, alpha, d)
%
%   For the Poisson and negative binomial methods (method1 = 1, 3, 4, 5) the
%   fitting objective fval is the negative log-likelihood WITHOUT the term that
%   depends only on the data, sum(log(y!)). That is harmless for fitting and
%   for the bootstrap, but it makes the AICc of these methods incomparable with
%   the AICc of least squares (method1 = 0), whose formula in getAICc is a
%   complete Gaussian log-likelihood.
%
%   adj is the amount to add to fval so that fval + adj is the full negative
%   log-likelihood of the observations that count toward the sample size:
%
%       adj = sum(gammaln(y + 1)) over the informative observations
%             - (objective of the fixed first point)            [fixI0 = 1 only]
%
%   With fixI0 = 1 the first fitted value equals the first observation by
%   construction; that point is not counted in n (see getAICc), so its
%   likelihood term is taken out as well. The term is evaluated with
%   growthObjectiveValue at mean = observation with the fitted alpha and d,
%   i.e. exactly what fval contains for that point.
%
%   adj is 0 for method1 = 0. fval itself is never changed.
%
%   (Ported from quantdiffAICcCountAdjustment in QuantDiffForecast.)

adj = 0;
info = struct('logFactorialTerm', 0, 'firstPointTerm', 0, 'nInformative', numel(ydata));

if ~any(method1 == [1 3 4 5])
    return
end
if nargin < 4 || isempty(alpha), alpha = 1; end
if nargin < 5 || isempty(d),     d = 1;     end

y = double(ydata(:));

pinned = false(size(y));
if ~isempty(fixI0) && fixI0 == 1
    pinned(1) = true;
end

logFact = sum(gammaln(max(y(~pinned), 0) + 1));

firstTerm = 0;
if any(pinned)
    y1 = y(pinned);
    mu1 = y1;
    mu1(mu1 == 0) = 0.001;   % the objective's own replacement for a zero mean
    firstTerm = growthObjectiveValue(mu1, y1, method1, alpha, d);
end

adj = logFact - firstTerm;
info = struct('logFactorialTerm', logFact, 'firstPointTerm', firstTerm, 'nInformative', nnz(~pinned));

end
