function [rows,cols]=forecastPanelGrid(nWindows)

% FORECASTPANELGRID  Rows and columns of the panel with one forecast per window.
%
%   [rows,cols] = forecastPanelGrid(nWindows)
%
% Input
%   nWindows  number of rolling windows (a positive integer)
%
% Outputs
%   rows, cols  layout for subplot(rows,cols,k), k = 1..nWindows, with
%               rows*cols = nWindows, so that every window has a panel and no
%               panel is left empty
%
% Rule
%   The layout comes from the prime factors of nWindows, f1 <= f2 <= ...:
%
%     one factor (1 or a prime)   nWindows-by-1
%     two factors                 f1-by-f2              (6 windows: 2-by-3)
%     three factors               (f1*f2)-by-f3         (12 windows: 4-by-3)
%     four or more factors        d-by-(nWindows/d), d the largest divisor of
%                                 nWindows not above its square root
%                                 (16: 4-by-4, 24: 4-by-6, 32: 4-by-8)
%
% History
%   plotForecast_GrowthModels.m held the first three cases inline and used
%   f1-by-f2 for every other count. With four or more prime factors (16, 24,
%   32, 36, 40, 48, ... windows) that grid has fewer panels than windows and
%   subplot stopped with an index error. Layouts for counts with at most three
%   prime factors are unchanged.

if ~isnumeric(nWindows) || ~isscalar(nWindows) || ~isreal(nWindows) || ~isfinite(nWindows) || nWindows<1 || nWindows~=round(nWindows)
    error('GrowthPredict:forecastPanelGrid:nWindows', ...
        'nWindows must be a positive integer.');
end

nWindows=double(nWindows);

factors=factor(nWindows);

if length(factors)==1

    rows=factors;
    cols=1;

elseif length(factors)==2

    rows=factors(1);
    cols=factors(2);

elseif length(factors)==3

    rows=factors(1)*factors(2);
    cols=factors(3);

else

    % Four or more prime factors: the pair of divisors closest to a square.
    rows=1;

    for d=2:nWindows

        if d*d>nWindows
            break
        end

        if mod(nWindows,d)==0
            rows=d;
        end

    end

    cols=nWindows/rows;

end
