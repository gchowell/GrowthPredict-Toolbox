function h = addEstimateLine(estimate, ax)
%ADDESTIMATELINE  Mark a point estimate on a histogram with a vertical dashed line.
%
%   h = addEstimateLine(estimate)
%   h = addEstimateLine(estimate, ax)
%
%   Draws a red dashed vertical line at x = estimate spanning the full
%   height of the axes (default: gca) and returns its handle. Used on the
%   bootstrap histograms of the parameters, where the red dashed horizontal
%   segment shows the 95% interval and this line shows the fit to the
%   observed data.
%
%   The y-limits are frozen at their value before the line is drawn, so the
%   line spans the axes without rescaling them. If the estimate lies outside
%   the current x-limits (possible when the fit is outside its own bootstrap
%   distribution), the x-limits are widened to include it with a 5% margin.
%   A non-finite estimate draws nothing and returns an empty handle.

if nargin < 2 || isempty(ax)
    ax = gca;
end

h = gobjects(0);

if ~isscalar(estimate) || ~isreal(estimate) || ~isfinite(estimate)
    return
end

holdState = ishold(ax);
hold(ax, 'on');

yl = ylim(ax);
xl = xlim(ax);

if estimate < xl(1) || estimate > xl(2)
    lo = min(xl(1), estimate);
    hi = max(xl(2), estimate);
    pad = 0.05 * max(hi - lo, eps(max(abs([lo hi]))));
    xlim(ax, [lo - pad, hi + pad]);
end

h = plot(ax, [estimate estimate], yl, 'r--', 'LineWidth', 2);

ylim(ax, yl);

if ~holdState
    hold(ax, 'off');
end

end
