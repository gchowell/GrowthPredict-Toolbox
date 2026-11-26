function h = setParameterTitle(name, values, ax)
%SETPARAMETERTITLE  Title of a parameter histogram: estimate and 95% CI on two lines.
%
%   h = setParameterTitle(name, values)
%   h = setParameterTitle(name, values, ax)
%
%   name    parameter name, e.g. 'r'
%   values  [estimate lower upper]
%   ax      axes (default: current axes)
%
%   Writes, for example,
%
%       r = 1.73
%       95% CI (1.49, 2.08)
%
%   above the axes, in normal weight, and sizes the text to the panel with
%   fitParameterTitles. The figure's SizeChangedFcn is set to
%   fitParameterTitles, unless the figure already has one, so the titles
%   keep fitting when the figure is resized.
%
%   Call it after the font size of the axes is set: setting the font size
%   of the axes afterwards resets the size of the title.
%
%   Before v1.1.0 the title was one bold line, e.g.
%   'r=1.7 (95% CI:1.5,2.1)', which is wider than one panel of a row of
%   four and ran into the titles next to it.
%
%   See also parameterTitleText, fitParameterTitles.

if nargin < 3 || isempty(ax)
    ax = gca;
end

h = title(ax, parameterTitleText(name, values));
set(h, 'FontWeight', 'normal', 'Tag', 'GrowthPredictParameterTitle');

fig = ancestor(ax, 'figure');

fitParameterTitles(fig);

try
    if isempty(get(fig, 'SizeChangedFcn'))
        set(fig, 'SizeChangedFcn', @fitParameterTitles);
    end
catch
    % Releases without SizeChangedFcn: the titles are sized once, above.
end

end
