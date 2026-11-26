% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>
%
function [fracLB, fracUB, hitLB, hitUB] = boundHitFractions(Phatss, LB, UB, relTol)
% BOUNDHITFRACTIONS  Share of bootstrap parameter estimates sitting on the optimizer box.
%
%   [fracLB, fracUB, hitLB, hitUB] = boundHitFractions(Phatss, LB, UB, relTol)
%
%   Bound-censoring diagnostic that complements PII_box / PII_ratio. The 2.5% and
%   97.5% quantiles behind the PII can hide a pile-up of estimates at ONE bound
%   (e.g., 80% of draws at LB with the rest interior gives a moderate-looking
%   PII_ratio although the estimate is entirely bound-determined). This function
%   reports, for every column of Phatss, the fraction of bootstrap draws that lie
%   within a relative tolerance of the lower and of the upper bound.
%
%   Inputs
%     Phatss  M-by-d matrix of bootstrap estimates (rows = draws, columns = parameters).
%     LB, UB  Optimizer box. Either 1-by-d vectors (one common box for all draws) or
%             M-by-d matrices (the per-draw box; needed when a bound is data-dependent,
%             e.g. the X0 upper bound sum(abs(data)) when params.fixI0==0).
%     relTol  Relative tolerance (default 1e-3, i.e. 0.1% of the bound value).
%             A draw P counts as hitting LB when  P <= LB + relTol*|LB|
%             and as hitting UB when              P >= UB - relTol*|UB|.
%             For a bound that is exactly 0 the tolerance is 0 (exact hit only).
%
%   Outputs
%     fracLB, fracUB  1-by-d fractions in [0,1]; NaN for pinned columns (LB==UB, i.e.
%                     fixed parameters or fixI0 initial conditions) and for columns
%                     whose draws are all NaN.
%     hitLB, hitUB    M-by-d logical masks of the individual hits.
%
%   Reading the output
%     fracLB+fracUB ~ 0   : the bootstrap distribution is interior; the box is not
%                           shaping the CI.
%     fracLB or fracUB large (e.g. > 0.05) while PII_ratio is moderate: the quantile-
%                           based PII understates bound censoring; the estimate is
%                           being pinned on one side by the box, not by the data.
%     both large          : draws split between the two bounds, typical of a flat or
%                           bimodal likelihood over the whole box.

if nargin < 4 || isempty(relTol)
    relTol = 1e-3;
end

[M, d] = size(Phatss);

LB = double(LB); UB = double(UB);

% Broadcast a common box (given as a row or column vector) to M-by-d. A vector whose
% length along dim 1 already equals M is a per-draw box with d==1 and is kept as is.
if isvector(LB) && size(LB,1) ~= M, LB = repmat(LB(:)', M, 1); end
if isvector(UB) && size(UB,1) ~= M, UB = repmat(UB(:)', M, 1); end

if size(LB,1) ~= M || size(UB,1) ~= M || size(LB,2) < d || size(UB,2) < d
    error('boundHitFractions:sizeMismatch', ...
        'LB/UB must be 1-by-d (or longer) or M-by-d (or wider) to match Phatss (%d-by-%d).', M, d);
end

% Trim a box longer than Phatss (the box also carries X0 / error-model entries).
LB = LB(:, 1:d);
UB = UB(:, 1:d);

% Guard against flipped boxes (fit_model already swaps them, but be safe).
flipMask = LB > UB;
if any(flipMask(:))
    tmp = LB(flipMask); LB(flipMask) = UB(flipMask); UB(flipMask) = tmp;
end

tolLB = relTol * abs(LB);
tolUB = relTol * abs(UB);

hitLB = Phatss <= LB + tolLB;
hitUB = Phatss >= UB - tolUB;

valid = ~isnan(Phatss);
hitLB = hitLB & valid;
hitUB = hitUB & valid;

nValid = sum(valid, 1);
fracLB = sum(hitLB, 1) ./ max(nValid, 1);
fracUB = sum(hitUB, 1) ./ max(nValid, 1);

% Pinned columns (LB==UB in every draw) and all-NaN columns carry no information.
pinned = all(LB == UB, 1);
fracLB(pinned | nValid == 0) = NaN;
fracUB(pinned | nValid == 0) = NaN;

end
