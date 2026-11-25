function [idxRows,reliable]=selectRtRows(giCoverage,minCoverage,stepsPerUnit)

% SELECTRTROWS  Reproduction-number estimates backed by enough incidence history.
%
%   [idxRows,reliable] = selectRtRows(giCoverage)
%   [idxRows,reliable] = selectRtRows(giCoverage,minCoverage,stepsPerUnit)
%
% Inputs
%   giCoverage    third output of get_Rt.m: for each estimate, the fraction of
%                 the generation-interval distribution covered by the incidence
%                 history available at that time
%   minCoverage   smallest fraction accepted (default 0.95)
%   stepsPerUnit  grid steps of get_Rt.m per unit of time (default 100)
%
% Outputs
%   reliable      logical column like giCoverage: giCoverage >= minCoverage.
%                 Use it to draw the curve: plot(ts(find(reliable)+1), Rs(reliable)).
%   idxRows       row vector of indices into the estimates for a table with one
%                 row per unit of time after the start of the curve, ending at
%                 the last time, restricted to reliable estimates. Estimate
%                 idxRows(k) belongs to time ts(idxRows(k)+1). When no estimate
%                 is reliable, idxRows is the last index alone, so that a table
%                 is never empty; its coverage shows that it is not reliable.
%
% Why
%   get_Rt.m uses only the incidence inside the curve. Near its start the
%   renewal denominator misses the infections that happened earlier, and the
%   estimate is too high: with the mpox default of options_Rt.m, about 2.5
%   times the true value one week in and 20% above it two weeks in when
%   incidence falls by 20% per week.

if nargin<2 || isempty(minCoverage)
    minCoverage=0.95;
end

if nargin<3 || isempty(stepsPerUnit)
    stepsPerUnit=100;
end

if ~isscalar(minCoverage) || ~isfinite(minCoverage) || minCoverage<=0 || minCoverage>1
    error('GrowthPredict:selectRtRows:minCoverage', ...
        'minCoverage must be a number above 0 and at most 1.');
end

if ~isscalar(stepsPerUnit) || ~isfinite(stepsPerUnit) || stepsPerUnit<1 || stepsPerUnit~=round(stepsPerUnit)
    error('GrowthPredict:selectRtRows:stepsPerUnit', ...
        'stepsPerUnit must be a positive integer.');
end

giCoverage=giCoverage(:);

n=numel(giCoverage);

reliable=giCoverage>=minCoverage;

idxRows=zeros(1,0);

if n==0
    return
end

% One index per unit of time after the start, and the last index.
idxAll=stepsPerUnit:stepsPerUnit:n;

if isempty(idxAll) || idxAll(end)~=n
    idxAll=[idxAll n];
end

keep=reliable(idxAll);

idxRows=idxAll(keep(:)');

if isempty(idxRows)
    idxRows=n;
end
