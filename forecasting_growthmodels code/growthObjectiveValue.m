function objfunction = growthObjectiveValue(fitcurve, observations, method1, alpha, d)
%GROWTHOBJECTIVEVALUE  Score a fitted incidence curve against the data.
%
%   objfunction = growthObjectiveValue(fitcurve, observations, method1, alpha, d)
%
%   method1 = 0  sum of squared residuals
%           = 1  Poisson negative log-likelihood
%           = 3  negative binomial NLL, VAR = mean + alpha*mean
%           = 4  negative binomial NLL, VAR = mean + alpha*mean^2
%           = 5  negative binomial NLL, VAR = mean + alpha*mean^d
%
%   The likelihoods leave out the term that depends only on the data,
%   sum(log(y!)); growthAICcCountAdjustment adds it back for AICc.
%
%   Same formulas as the former plotModifiedLogisticGrowthMethods1, with two
%   changes (ported from quantdiffObjectiveValue in QuantDiffForecast):
%     * the negative binomial sums  sum_{j=0}^{y-1} log(j+m)  are evaluated as
%       gammaln(y+m) - gammaln(m). For integer y the value is identical; the
%       cost no longer grows with the number of cases. For non-integer y the
%       gamma form is the continuous extension of the likelihood (the loop
%       silently truncated y).
%     * invalid input or a non-finite / complex score raises an error instead
%       of returning a number; growthObjective maps that to Inf.
%   Exact zeros of the fitted curve are replaced by 0.001, as before.

if ~isvector(fitcurve) || ~isvector(observations) || isempty(fitcurve) || ...
        numel(fitcurve) ~= numel(observations) || ...
        ~isreal(fitcurve) || ~isreal(observations) || ...
        any(~isfinite(fitcurve(:))) || any(~isfinite(observations(:)))
    error('GrowthPredict:growthObjectiveValue:InvalidInput', ...
        'A score requires matching finite real observation and curve vectors.');
end

yfit  = fitcurve(:);
ydata = observations(:);

if any(yfit < 0)
    error('GrowthPredict:growthObjectiveValue:NegativeCurve', ...
        'The fitted incidence curve must be nonnegative.');
end
if sum(yfit) == 0
    error('GrowthPredict:growthObjectiveValue:ZeroCurve', ...
        'The fitted curve is identically zero.');
end

yfit(yfit == 0) = 0.001;   % allows the logarithms below
yg = max(ydata, 0);        % argument of the gamma terms

switch method1

    case 0  % least squares

        objfunction = sum((ydata - yfit).^2);

    case 1  % Poisson

        objfunction = -sum(ydata .* log(yfit) - yfit);

    case 3  % negative binomial, VAR = mean + alpha*mean

        m = yfit ./ alpha;
        objfunction = -sum(gammaln(yg + m) - gammaln(m) ...
            + ydata .* log(alpha) - (ydata + m) .* log(1 + alpha));

    case 4  % negative binomial, VAR = mean + alpha*mean^2

        m = 1 / alpha;
        objfunction = -sum(gammaln(yg + m) - gammaln(m) ...
            + ydata .* log(alpha .* yfit) - (ydata + m) .* log(1 + alpha .* yfit));

    case 5  % negative binomial, VAR = mean + alpha*mean^d

        m = (1 ./ alpha) .* yfit.^(2 - d);
        objfunction = -sum(gammaln(yg + m) - gammaln(m) ...
            + ydata .* log(alpha .* yfit.^(d - 1)) ...
            - (ydata + m) .* log(1 + alpha .* yfit.^(d - 1)));

    otherwise
        error('GrowthPredict:growthObjectiveValue:UnsupportedMethod', ...
            'method1 = %g is not supported; use 0, 1, 3, 4 or 5.', method1);

end

if ~isscalar(objfunction) || ~isreal(objfunction) || ~isfinite(objfunction)
    error('GrowthPredict:growthObjectiveValue:InvalidObjective', ...
        'The parameters do not give a finite real objective.');
end

end
