function [AICc,part1,part2,numparams,nEff]=getAICc(method1,dist1,flag1,fixI0,fval,n,fvalAdjust)
%GETAICC  Corrected Akaike information criterion of a growth-model fit.
%
%   [AICc, part1, part2, numparams, nEff] = getAICc(method1, dist1, flag1, fixI0, fval, n)
%   [...] = getAICc(method1, dist1, flag1, fixI0, fval, n, fvalAdjust)
%
%   fval        objective value at the fit (sum of squares for method1 = 0,
%               negative log-likelihood without the sum(log(y!)) term otherwise)
%   n           number of observations in the calibration window
%   fvalAdjust  (count likelihoods) amount that completes fval to the full
%               negative log-likelihood, from growthAICcCountAdjustment.
%               Omitted or empty -> 0, i.e. the value reported before
%               November 2026.
%
%   AICc = part1 + part2, part1 = -2 log L, part2 = 2k + 2k(k+1)/(nEff-k-1).
%
%   Revised November 2026 (following QuantDiffForecast):
%     * Least squares uses the complete Gaussian log-likelihood,
%           -2 log L = nEff*log(SSE/nEff) + nEff*log(2*pi) + nEff.
%       The values are now comparable with the count likelihoods (given
%       fvalAdjust).
%     * nEff = n - 1 when fixI0 = 1: the first fitted value equals the first
%       observation by construction, so that point carries no information.
%     * Effect on least-squares model ranking within one window: with
%       fixI0 = 0 the new value differs from the former n*log(SSE) + ... by
%       a constant, so differences between models are unchanged. With
%       fixI0 = 1 the sample size is n - 1 instead of n, so differences
%       change by a factor (n-1)/n in the fit term and slightly in the
%       small-sample term; close calls can change order.
%     * The number of parameters comes from get_nparams (one definition).
%     * Too few data points is an error at the call site of the runners; here
%       AICc = Inf with a warning, as before.

if nargin < 7 || isempty(fvalAdjust)
    fvalAdjust = 0;
end

numparams=get_nparams(method1,dist1,flag1,fixI0);

nEff=n;
if ~isempty(fixI0) && fixI0==1
    nEff=n-1;
end

if nEff <= numparams + 1

    AICc = Inf;
    part1 = NaN;
    part2 = Inf;

    warning('GrowthPredict:getAICc:InsufficientSample', ...
        ['AICc is not defined for n=%g informative observations and k=%g; ' ...
        'returning Inf (requires n>k+1).'], ...
        nEff,numparams);

    return
end

part2=2*numparams + (2*numparams*(numparams+1))/(nEff-numparams-1);

switch method1

    case 0

        part1=nEff*log(fval/nEff) + nEff*log(2*pi) + nEff;

    otherwise

        part1=2*(fval+fvalAdjust);

end

AICc=part1+part2;

end
