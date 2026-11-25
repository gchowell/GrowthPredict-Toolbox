function tests = test_doublingTime
%TEST_DOUBLINGTIME  Regression tests for getDoublingTimeCurve.m.
%
%   Run from the code folder:
%       results = runtests('tests/test_doublingTime.m');
%       table(results)
%
%   Reference values marked 'Python' come from an independent Python
%   re-implementation of getDoublingTimeCurve.m. Times are found on a grid of
%   0.01, so they are compared with a tolerance of one grid step.

tests = functiontests(localfunctions);

end


function setupOnce(testCase)

codeDir = fileparts(fileparts(mfilename('fullpath')));
addpath(codeDir);

testCase.TestData.codeDir = codeDir;

data = load(fullfile(codeDir, 'input', 'Most_Recent_Timeseries_US-CDC.txt'));
testCase.TestData.mpox = data(:, 2);

end


function teardownOnce(testCase)

rmpath(testCase.TestData.codeDir);

end


% <============================================================================>
% <================================ Tests =====================================>
% <============================================================================>

function testNoPriorCasesIsTheFormerBehaviour(testCase)
% A window that starts at the first row of the series has no earlier cases:
% three inputs, and a fourth input of zero, give the former result.

y = testCase.TestData.mpox(1:20);

[tds3, C03, curve3, dt3] = getDoublingTimeCurve(y, 1, 0);
[tds4, C04, curve4, dt4] = getDoublingTimeCurve(y, 1, 0, 0);

verifyEqual(testCase, tds4, tds3);
verifyEqual(testCase, C04, C03);
verifyEqual(testCase, curve4, curve3);
verifyEqual(testCase, dt4, dt3);

verifyEqual(testCase, C03, 14);
verifyEqual(testCase, dt3, [0.87 0.89 0.84 0.93 1.06 1.12 1.04 1.17 1.54 2.40]', 'AbsTol', 0.011);    % Python

end


function testExponentialGrowthInsideAnEpidemic(testCase)
% Incidence exp(0.2*k), k = 0..39, window k = 20..39. With the earlier cases
% the count doubles every ln(2)/0.2 = 3.47 time units; when the count restarts
% at the window, the first 'doubling' takes 0.81.

I = exp(0.2*(0:39))';

window1 = I(21:40);
prior1 = sum(I(1:20));

[~, ~, ~, dtNew] = getDoublingTimeCurve(window1, 1, 0, prior1);
[~, ~, ~, dtOld] = getDoublingTimeCurve(window1, 1, 0);

verifyEqual(testCase, dtNew, [3.40 3.46 3.44 3.47 3.46]', 'AbsTol', 0.011);    % Python
verifyLessThan(testCase, max(abs(dtNew - log(2)/0.2)), 0.08);

verifyEqual(testCase, dtOld(1:3), [0.81 1.34 1.93]', 'AbsTol', 0.011);         % Python

% Same result as following the count of the whole series from the window start.
[~, ~, ~, dtWhole] = getDoublingTimeCurve([sum(I(1:21)); I(22:40)], 1, 0);

verifyEqual(testCase, dtNew, dtWhole, 'AbsTol', 1e-9);

end


function testFallingIncidenceHasNoDoubling(testCase)
% Logistic epidemic of 20000 cases, window starting at its midpoint: the count
% goes from 10000 to 19999 and never doubles. Restarting the count at the
% window reported doublings after 1.00 and 2.78 time units.

k = (0:40)';
C = 20000./(1 + exp(-0.5*(k - 12)));
I = [C(1); diff(C)];

window1 = I(13:32);
prior1 = sum(I(1:12));

[tds, C0data, ~, dtNew] = getDoublingTimeCurve(window1, 1, 0, prior1);
[~, ~, ~, dtOld] = getDoublingTimeCurve(window1, 1, 0);

verifyEmpty(testCase, dtNew);
verifyEqual(testCase, C0data, 10000, 'RelTol', 1e-9);
verifyEqual(testCase, size(tds, 1), 1);

verifyEqual(testCase, dtOld, [1.00 2.78]', 'AbsTol', 0.011);    % Python

% A window in the rising phase still doubles, more slowly than its own running total.
[~, ~, ~, dtNew] = getDoublingTimeCurve(I(7:26), 1, 0, sum(I(1:6)));

verifyEqual(testCase, dtNew, [1.43 1.65 1.93 3.31]', 'AbsTol', 0.011);    % Python

end


function testMpoxWindows(testCase)
% Observed U.S. mpox counts. Rows 13-32: 11725 cases before the window; the
% count goes from 14752 to 29751 and doubles once, after 15.09 weeks.

y = testCase.TestData.mpox;

prior1 = sum(y(1:12));

verifyEqual(testCase, prior1, 11725);

[tds, C0data, ~, dtNew] = getDoublingTimeCurve(y(13:32), 1, 0, prior1);
[~, ~, ~, dtOld] = getDoublingTimeCurve(y(13:32), 1, 0);

verifyEqual(testCase, C0data, 14752);
verifyEqual(testCase, dtNew, 15.09, 'AbsTol', 0.011);            % Python
verifyEqual(testCase, tds(end, 1), 15.09, 'AbsTol', 0.011);
verifyLessThanOrEqual(testCase, tds(end, 2), 2*14752);

verifyEqual(testCase, dtOld, [1.14 3.09]', 'AbsTol', 0.011);     % Python

% Rows 8-27: 1033 cases before the window.
[~, ~, ~, dtNew] = getDoublingTimeCurve(y(8:27), 1, 0, sum(y(1:7)));

verifyEqual(testCase, dtNew, [1.17 1.67 2.76]', 'AbsTol', 0.011);    % Python

end


function testInputHandling(testCase)

% One point: nothing to follow.
[tds, C0data, ~, doublingtimes] = getDoublingTimeCurve(5, 1, 0);

verifyEmpty(testCase, tds);
verifyEmpty(testCase, C0data);
verifyEmpty(testCase, doublingtimes);

% A row vector is accepted.
[~, ~, ~, dtRow] = getDoublingTimeCurve(exp(0.2*(0:19)), 1, 0);
[~, ~, ~, dtColumn] = getDoublingTimeCurve(exp(0.2*(0:19))', 1, 0);

verifyEqual(testCase, dtRow, dtColumn);

% The count before the window cannot be negative or missing.
verifyError(testCase, @() getDoublingTimeCurve(ones(10, 1), 1, 0, -1), ...
    'GrowthPredict:getDoublingTimeCurve:priorCumulative');
verifyError(testCase, @() getDoublingTimeCurve(ones(10, 1), 1, 0, NaN), ...
    'GrowthPredict:getDoublingTimeCurve:priorCumulative');

end
