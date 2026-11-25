function fontSize = fitParameterTitles(fig, varargin)
%FITPARAMETERTITLES  Size the parameter-histogram titles of a figure to their panels.
%
%   fontSize = fitParameterTitles(fig)
%
%   Finds the titles created by setParameterTitle in figure fig (default:
%   the current figure) and gives them one common font size, in points, at
%   which the longest title is no wider than the space of its panel. The
%   function is also the SizeChangedFcn of those figures, so the titles
%   follow the figure when it is resized; extra inputs (the callback's
%   event data) are ignored. fontSize is empty when the figure has no such
%   title.
%
%   Font size, for each title, with W the width of its axes in points,
%   n the number of characters of its longest line, F the font size of its
%   axes and H the height of the figure in points:
%
%       min( 0.75*F,  1.25*W/(0.6*n),  0.031*H ),   not below 6 points
%
%   The smallest value over the titles is used for all of them. The second
%   term keeps a line of n characters, at 0.6 of the font size per
%   character, within 1.25 axes widths; neighbouring panels of
%   subplot(2,4,k) are 1.32 axes widths apart. The third term keeps two
%   lines within the margin above the top row. The rendered width is then
%   measured and the size reduced if a title is wider than 1.15 axes
%   widths, which covers fonts wider than Helvetica or Arial.
%
%   See also setParameterTitle, parameterTitleText.

fontSize = [];

if nargin < 1 || isempty(fig)
    fig = get(0, 'CurrentFigure');
end

if isempty(fig) || ~all(ishghandle(fig))
    return
end

fig = fig(1);

titles = findall(fig, 'Type', 'text', 'Tag', 'GrowthPredictParameterTitle');

if isempty(titles)
    return
end

capFraction = 0.75;      % of the axes font size
roomFraction = 1.25;     % of the axes width, by character count
charWidth = 0.6;         % of the font size, per character (an upper value for digits)
heightFraction = 0.031;  % of the figure height
maxFraction = 1.15;      % of the axes width, as rendered
minFontSize = 6;

try

    pointsPerPixel = 72 / get(0, 'ScreenPixelsPerInch');
    figPos = getpixelposition(fig);
    figHeight = figPos(4) * pointsPerPixel;

    sizes = zeros(numel(titles), 1);
    axesOf = cell(numel(titles), 1);

    for k = 1:numel(titles)

        ax = ancestor(titles(k), 'axes');
        axesOf{k} = ax;

        axPos = getpixelposition(ax);
        axWidth = axPos(3) * pointsPerPixel;

        nchar = max(1, longestLine(get(titles(k), 'String')));

        sizes(k) = min([capFraction * get(ax, 'FontSize'), ...
            roomFraction * axWidth / (charWidth * nchar), ...
            heightFraction * figHeight]);

    end

    fontSize = max(minFontSize, min(sizes));

    set(titles, 'FontUnits', 'points', 'FontSize', fontSize);

catch

    % A figure without a usable size (e.g. while it is being deleted):
    % leave the titles as they are.
    fontSize = [];
    return

end

% Rendered width, as a fraction of the axes width.
try

    widest = 0;

    for k = 1:numel(titles)

        ax = axesOf{k};

        if strcmp(get(ax, 'XScale'), 'linear') && strcmp(get(titles(k), 'Units'), 'data')

            extent = get(titles(k), 'Extent');
            xl = get(ax, 'XLim');
            widest = max(widest, extent(3) / (xl(2) - xl(1)));

        end

    end

    if isfinite(widest) && widest > maxFraction

        fontSize = max(minFontSize, fontSize * maxFraction / widest);
        set(titles, 'FontSize', fontSize);

    end

catch
    % Keep the size from the character count.
end

end


function n = longestLine(str)
%LONGESTLINE  Number of characters of the longest line of a title string.

if iscell(str)
    n = max(cellfun(@(s) numel(char(s)), str));
elseif ischar(str)
    n = size(str, 2);
else
    n = max(strlength(string(str)));
end

if isempty(n)
    n = 0;
end

end
