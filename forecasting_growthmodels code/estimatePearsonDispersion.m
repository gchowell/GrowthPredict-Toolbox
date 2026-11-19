function [phi, info] = estimatePearsonDispersion(ydata, fitcurve, numModelParams, nSeries, fixI0, minMean)
%ESTIMATEPEARSONDISPERSION  Variance-to-mean ratio of the data around the fitted curve.
%
%   [phi, info] = estimatePearsonDispersion(ydata, fitcurve, numModelParams, nSeries, fixI0)
%   [phi, info] = estimatePearsonDispersion(..., minMean)
%
%   Pearson (quasi-Poisson) dispersion estimate for the noise model
%   VAR = phi * mean, used to scale the negative binomial bootstrap noise of
%   AddErrorStructure (dist1 = 2) when the model is fitted by least squares
%   (method1 = 0):
%
%       phi = sum( (y_i - mu_i)^2 / mu_i ) / (n_used - numModelParams)
%
%   with mu = fitcurve. The fitted curve removes the epidemic trend, so phi
%   measures the scatter around the curve only. (The former estimate, the
%   var/mean ratio of 7-point bins of the raw data, also counted the trend
%   inside each bin and was inflated, e.g. about 5.6 for Poisson data.)
%
%   ydata, fitcurve  stacked observations and fitted curve (same length; one
%                    block of equal length per fitted series, as in the runners)
%   numModelParams   estimated quantities of the CURVE (model parameters and free
%                    initial values), not counting the Gaussian variance of the
%                    least-squares fit; the runners pass numparams - 1
%   nSeries          number of fitted series stacked in ydata (default 1)
%   fixI0            1 -> the first point of each series is fixed to the data
%                    by construction (zero residual) and is excluded (default 0)
%   minMean          points with fitted mean below this are excluded, since
%                    (y-mu)^2/mu is unstable for tiny mu (default 0.5). The
%                    default assumes COUNT data; for rates or proportions the
%                    VAR = phi*mean model itself does not apply
%
%   phi   the estimate; NaN when fewer than numModelParams+1 points are usable
%   info  struct: nUsed, dof, minMean, excludedSmallMean, excludedFixed
%
%   The caller decides what to do with phi <= 1 or NaN (the runners fall
%   back to Poisson noise).
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if nargin < 4 || isempty(nSeries), nSeries = 1; end
if nargin < 5 || isempty(fixI0),   fixI0 = 0;   end
if nargin < 6 || isempty(minMean), minMean = 0.5; end

y  = double(ydata(:));
mu = double(fitcurve(:));
if numel(y) ~= numel(mu)
    error('estimatePearsonDispersion:size', ...
        'ydata (%d) and fitcurve (%d) must have the same number of elements.', numel(y), numel(mu));
end
if mod(numel(y), nSeries) ~= 0
    error('estimatePearsonDispersion:series', ...
        '%d observations cannot be split into %d series of equal length.', numel(y), nSeries);
end
nT = numel(y) / nSeries;

use = isfinite(y) & isfinite(mu);
fixedMask = false(size(y));
if fixI0
    fixedMask(1:nT:end) = true;          % first point of every series
end
smallMask = mu < minMean;
use = use & ~fixedMask & ~smallMask;

dof = nnz(use) - numModelParams;
if dof < 1
    phi = NaN;
else
    phi = sum((y(use) - mu(use)).^2 ./ mu(use)) / dof;
end

info = struct('nUsed', nnz(use), 'dof', dof, 'minMean', minMean, ...
    'excludedSmallMean', nnz(smallMask & ~fixedMask), 'excludedFixed', nnz(fixedMask));
end
