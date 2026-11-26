function tests = test_parameterTitle
%TEST_PARAMETERTITLE  Titles of the parameter histograms (r, p, a, K).
%
%   The titles are written on two lines and sized to their panels by
%   parameterTitleText, setParameterTitle and fitParameterTitles. Before
%   v1.1.0 each title was one bold line that ran into its neighbours.
%
%   Run from the code folder:  runtests('tests')

tests = functiontests(localfunctions);

end


function setupOnce(testCase)

testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename('fullpath')))));

end


% ------------------------------------------------------------------ text

function testThreeSignificantDigits(testCase)

inputs   = [1.7302 0.8526 1   28547.1 801.8 0.17743 0.0012345 999.6  9.996 123456.7 10   0.5   -0.52];
expected = {'1.73' '0.853' '1' '28500' '802' '0.177' '0.00123' '1000' '10' '123000' '10' '0.5' '-0.52'};

for k = 1:numel(inputs)
    lines = parameterTitleText('r', [inputs(k) 0 1]);
    verifyEqual(testCase, lines{1}, ['r = ' expected{k}], sprintf('input %g', inputs(k)));
end

end


function testExponentOutsideRange(testCase)

lines = parameterTitleText('K', [1234567 2.85e8 3.2e-5]);

verifyEqual(testCase, lines{1}, 'K = 1.23e6');
verifyEqual(testCase, lines{2}, '95% CI (2.85e8, 3.20e-5)');

end


function testQuickStartTitles(testCase)

verifyEqual(testCase, parameterTitleText('r', [1.7302 1.49 2.08]), {'r = 1.73', '95% CI (1.49, 2.08)'});
verifyEqual(testCase, parameterTitleText('p', [0.8526 0.828 0.874]), {'p = 0.853', '95% CI (0.828, 0.874)'});
verifyEqual(testCase, parameterTitleText('a', [1 1 1]), {'a = 1', '95% CI (1, 1)'});
verifyEqual(testCase, parameterTitleText('K', [28547 27295 29915]), {'K = 28500', '95% CI (27300, 29900)'});

end


function testSpecialValues(testCase)

lines = parameterTitleText('a', [0 NaN Inf]);

verifyEqual(testCase, lines, {'a = 0', '95% CI (NaN, Inf)'});

lines = parameterTitleText('a', [0 -Inf 0]);

verifyEqual(testCase, lines{2}, '95% CI (-Inf, 0)');

end


function testColumnInputAndErrors(testCase)

verifyEqual(testCase, parameterTitleText('r', [1.7302; 1.49; 2.08]), {'r = 1.73', '95% CI (1.49, 2.08)'});

verifyError(testCase, @() parameterTitleText('r', [1 2]), 'GrowthPredict:parameterTitleText:values');
verifyError(testCase, @() parameterTitleText('r', 'abc'), 'GrowthPredict:parameterTitleText:values');
verifyError(testCase, @() parameterTitleText(1, [1 2 3]), 'GrowthPredict:parameterTitleText:name');

end


% ---------------------------------------------------------------- layout

function testTitlesDoNotOverlap(testCase)

% Size of docs/images/fit.png and forecast.png (12 x 6.6 in), the default
% figure size, and a large figure.
sizesPt = [864 475; 560 420; 1200 700];

sets = {[1.7302 1.49 2.08; 0.8526 0.828 0.874; 1 1 1; 28547 27295 29915], ...
    [801.8 571 1210; 0.178 0.128 0.223; 1 1 1; 18153 18026 18470], ...
    [0.00123 0.00101 0.00187; 0.999 0.951 1; 4.07 0.0123 9.87; 2.85e8 1.02e8 9.91e8]};

for s = 1:size(sizesPt, 1)
    for v = 1:numel(sets)

        [fig, ax] = histogramFigure(testCase, sizesPt(s, :));
        titleHandles = newTitles(ax, sets{v});
        drawnow

        [left, right] = titleEdges(ax, titleHandles);
        label = sprintf('figure %d x %d pt, value set %d', sizesPt(s, 1), sizesPt(s, 2), v);

        verifyGreaterThan(testCase, left(2:4) - right(1:3), 0, ['titles overlap: ' label]);
        verifyLessThan(testCase, right(4), 1, ['last title leaves the figure: ' label]);
        verifyGreaterThan(testCase, left(1), 0, ['first title leaves the figure: ' label]);

        close(fig);

    end
end

end


function testOldTitlesOverlapped(testCase)

% The check above must be able to see an overlap: the one-line bold titles
% used before v1.1.0 overlap at the size of the README figures.
[~, ax] = histogramFigure(testCase, [864 475]);

old = {'r=1.7 (95% CI:1.5,2.1)', 'p=0.85 (95% CI:0.83,0.87)', 'a=1 (95% CI:1,1)', 'K=2.9e+04 (95% CI:2.7e+04,3e+04)'};
titleHandles = gobjects(1, 4);

for k = 1:4
    titleHandles(k) = title(ax(k), old{k});
end

drawnow

[left, right] = titleEdges(ax, titleHandles);

verifyLessThan(testCase, left(2) - right(1), 0);
verifyLessThan(testCase, left(4) - right(3), 0);

end


function testOneFontSizeWithinLimits(testCase)

[fig, ax] = histogramFigure(testCase, [864 475]);
titleHandles = newTitles(ax, [1.7302 1.49 2.08; 0.8526 0.828 0.874; 1 1 1; 28547 27295 29915]);

drawnow

fontSize = fitParameterTitles(fig);
sizes = arrayfun(@(h) get(h, 'FontSize'), titleHandles);

verifyEqual(testCase, sizes, repmat(fontSize, 1, 4), 'AbsTol', 1e-9, 'the four titles share one font size');
verifyLessThanOrEqual(testCase, fontSize, 0.75 * 24 + 1e-9);
verifyGreaterThanOrEqual(testCase, fontSize, 6);

% Rule by character count for this figure: axes 0.1566*864 pt wide,
% longest line '95% CI (0.828, 0.874)' and '95% CI (27300, 29900)', 21 characters.
axPos = get(ax(1), 'Position');
expected = min([0.75 * 24, 1.25 * axPos(3) * 864 / (0.6 * 21), 0.031 * 475]);

verifyLessThanOrEqual(testCase, fontSize, expected + 0.05);
verifyGreaterThan(testCase, fontSize, 0.6 * expected, 'the size is close to the rule, not at the floor');

verifyTrue(testCase, all(strcmp(get(titleHandles, 'FontWeight'), 'normal')));

end


function testTitlesFollowFigureSize(testCase)

[fig, ax] = histogramFigure(testCase, [1440 800]);
titleHandles = newTitles(ax, [801.8 571 1210; 0.178 0.128 0.223; 1 1 1; 18153 18026 18470]);

drawnow

large = fitParameterTitles(fig);

set(fig, 'Position', [20 20 480 400]);
drawnow

small = fitParameterTitles(fig);

verifyLessThan(testCase, small, large);

[left, right] = titleEdges(ax, titleHandles);

verifyGreaterThan(testCase, left(2:4) - right(1:3), 0);

end


function testFigureWithoutParameterTitles(testCase)

fig = figure('Visible', 'off');
testCase.addTeardown(@() closeIfOpen(fig));

axes('Parent', fig);
title('another title');

verifyEmpty(testCase, fitParameterTitles(fig));

end


% ---------------------------------------------------------------- helpers

function [fig, ax] = histogramFigure(testCase, sizePt)
%HISTOGRAMFIGURE  Invisible figure with the top row of the toolbox figures.

fig = figure('Visible', 'off', 'Units', 'points', 'Position', [20 20 sizePt]);
testCase.addTeardown(@() closeIfOpen(fig));

values = [0.9 1.0 1.0 1.1 1.1 1.1 1.2 1.2 1.3 1.4];
names = {'r', 'p', 'a', 'K'};
ax = gobjects(1, 4);

for k = 1:4
    ax(k) = subplot(2, 4, k, 'Parent', fig);
    hist(ax(k), values);
    xlabel(ax(k), names{k});
    set(ax(k), 'FontSize', 24);
end

subplot(2, 4, [5 6 7 8], 'Parent', fig);

end


function titleHandles = newTitles(ax, values)

names = {'r', 'p', 'a', 'K'};
titleHandles = gobjects(1, 4);

for k = 1:4
    titleHandles(k) = setParameterTitle(names{k}, values(k, :), ax(k));
end

end


function [left, right] = titleEdges(ax, titleHandles)
%TITLEEDGES  Left and right edge of each title, as fractions of the figure width.

left = zeros(1, numel(ax));
right = zeros(1, numel(ax));

for k = 1:numel(ax)

    set(ax(k), 'Units', 'normalized');
    axPos = get(ax(k), 'Position');
    xl = get(ax(k), 'XLim');
    extent = get(titleHandles(k), 'Extent');          % data units of the axes

    left(k) = axPos(1) + (extent(1) - xl(1)) / (xl(2) - xl(1)) * axPos(3);
    right(k) = axPos(1) + (extent(1) + extent(3) - xl(1)) / (xl(2) - xl(1)) * axPos(3);

end

end


function closeIfOpen(fig)

if ishghandle(fig)
    close(fig);
end

end
