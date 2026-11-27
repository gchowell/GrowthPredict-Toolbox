function tests = test_forecastPanelGrid
%TEST_FORECASTPANELGRID  Regression tests for forecastPanelGrid.m.
%
%   Run from the code folder:
%       results = runtests('tests/test_forecastPanelGrid.m');
%       table(results)
%
%   forecastPanelGrid.m gives the rows and columns of figure 400 of
%   plotForecast_GrowthModels.m, which has one forecast panel per rolling
%   window. Reference layouts marked 'Python' come from a Python
%   re-implementation (docs/validation/python/plot_port.py).
%
%   No toolbox is required.

tests = functiontests(localfunctions);

end


function setupOnce(testCase)

codeDir = fileparts(fileparts(mfilename('fullpath')));
addpath(codeDir);

testCase.TestData.codeDir = codeDir;

end


function teardownOnce(testCase)

rmpath(testCase.TestData.codeDir);

end


% <============================================================================>
% <================================ Tests =====================================>
% <============================================================================>

function testEveryWindowHasAPanel(testCase)
% rows*cols equals the number of windows, so subplot(rows,cols,k) exists for
% every window and no panel is left empty.

for n = 1:200

    [rows, cols] = forecastPanelGrid(n);

    verifyEqual(testCase, rows*cols, n, sprintf('%d windows', n));
    verifyEqual(testCase, [rows cols], round([rows cols]), sprintf('%d windows', n));
    verifyGreaterThanOrEqual(testCase, min(rows, cols), 1, sprintf('%d windows', n));

end

end


function testCountsThatUsedToFail(testCase)
% Four or more prime factors: the former grid was f1-by-f2 = 2-by-2 for 16, 24
% and 32 windows and 2-by-3 for 54, fewer panels than windows.

counts = [16 24 32 36 40 48 54 56 60 64 72 80 81 96 100];
expected = [4 4; 4 6; 4 8; 6 6; 5 8; 6 8; 6 9; 7 8; 6 10; 8 8; 8 9; 8 10; 9 9; 8 12; 10 10];    % Python

for k = 1:numel(counts)

    [rows, cols] = forecastPanelGrid(counts(k));

    verifyEqual(testCase, [rows cols], expected(k, :), sprintf('%d windows', counts(k)));

    f = factor(counts(k));
    verifyLessThan(testCase, f(1)*f(2), counts(k), 'the former grid was too small for this count');

end

end


function testOtherCountsAreUnchanged(testCase)
% One, two or three prime factors: the layout plotForecast_GrowthModels.m used
% before. The README example with six windows stays 2-by-3.

for n = 1:200

    f = factor(n);

    if numel(f) == 1
        expected = [f 1];
    elseif numel(f) == 2
        expected = [f(1) f(2)];
    elseif numel(f) == 3
        expected = [f(1)*f(2) f(3)];
    else
        continue
    end

    [rows, cols] = forecastPanelGrid(n);

    verifyEqual(testCase, [rows cols], expected, sprintf('%d windows', n));

end

[rows, cols] = forecastPanelGrid(6);
verifyEqual(testCase, [rows cols], [2 3]);

[rows, cols] = forecastPanelGrid(1);
verifyEqual(testCase, [rows cols], [1 1]);

[rows, cols] = forecastPanelGrid(7);
verifyEqual(testCase, [rows cols], [7 1]);

[rows, cols] = forecastPanelGrid(12);
verifyEqual(testCase, [rows cols], [4 3]);

end


function testNewLayoutsAreTheMostSquare(testCase)
% For four or more prime factors rows is the largest divisor not above the
% square root, so rows <= cols and no other divisor pair is closer to a square.

for n = 1:400

    if numel(factor(n)) < 4
        continue
    end

    [rows, cols] = forecastPanelGrid(n);

    verifyLessThanOrEqual(testCase, rows, cols, sprintf('%d windows', n));

    divisors = find(mod(n, 1:n) == 0);
    verifyEqual(testCase, rows, max(divisors(divisors.^2 <= n)), sprintf('%d windows', n));

end

end


function testSubplotAcceptsEveryPanel(testCase)
% The calls plotForecast_GrowthModels.m makes: one panel per window, then the
% bottom row and the first panel again for the axis labels.

for n = [6 12 16 24 32]

    [rows, cols] = forecastPanelGrid(n);

    fig = figure('Visible', 'off');
    testCase.addTeardown(@() closeIfOpen(fig));

    for k = 1:n
        subplot(rows, cols, k, 'Parent', fig);
    end

    for c = 1:cols
        subplot(rows, cols, (rows - 1)*cols + c, 'Parent', fig);
    end

    subplot(rows, cols, 1, 'Parent', fig);

    verifyEqual(testCase, numel(findobj(fig, 'Type', 'axes')), n, sprintf('%d windows', n));

    close(fig);

end

end


function testInputChecks(testCase)

verifyError(testCase, @() forecastPanelGrid(0), 'GrowthPredict:forecastPanelGrid:nWindows');
verifyError(testCase, @() forecastPanelGrid(-3), 'GrowthPredict:forecastPanelGrid:nWindows');
verifyError(testCase, @() forecastPanelGrid(2.5), 'GrowthPredict:forecastPanelGrid:nWindows');
verifyError(testCase, @() forecastPanelGrid([4 6]), 'GrowthPredict:forecastPanelGrid:nWindows');
verifyError(testCase, @() forecastPanelGrid(NaN), 'GrowthPredict:forecastPanelGrid:nWindows');
verifyError(testCase, @() forecastPanelGrid('6'), 'GrowthPredict:forecastPanelGrid:nWindows');

% An integer type is accepted.
[rows, cols] = forecastPanelGrid(int32(16));
verifyEqual(testCase, [rows cols], [4 4]);

end


% <============================================================================>
% <============================== Local helper ================================>
% <============================================================================>

function closeIfOpen(fig)

if ishghandle(fig)
    close(fig);
end

end
