function tests = test_selectRtRows
%TEST_SELECTRTROWS  Regression tests for selectRtRows.m (coverage filter for Rt).
%
%   Run from the code folder:
%       results = runtests('tests/test_selectRtRows.m');
%       table(results)
%
%   Reference values marked 'Python' come from an independent Python
%   re-implementation of get_Rt.m and selectRtRows.m.
%
%   Requires the Statistics and Machine Learning Toolbox.

tests = functiontests(localfunctions);

end


function setupOnce(testCase)

codeDir = fileparts(fileparts(mfilename('fullpath')));
addpath(codeDir);

testCase.TestData.codeDir = codeDir;

% Mpox serial interval of options_Rt.m and cholera serial interval, in weeks.
testCase.TestData.mpox = [8.5/7, (5.0/7)^2];
testCase.TestData.cholera = [5/7, (8/7)^2];

end


function teardownOnce(testCase)

rmpath(testCase.TestData.codeDir);

end


% <============================================================================>
% <================================ Tests =====================================>
% <============================================================================>

function testMpoxDefaultDropsFirstTwoWeeks(testCase)
% Window labelled 12..31. Coverage reaches 95% at elapsed time 2.58, so the
% rows for t = 13 and t = 14 are left out.

gi = testCase.TestData.mpox;

t = 12:31;
I = 100*exp(-0.2*(0:19));

[Rs, ts, giCoverage] = get_Rt(t, I, 1, gi(1), gi(2));

[idxRows, reliable] = selectRtRows(giCoverage);

verifySize(testCase, reliable, size(Rs));
verifyEqual(testCase, find(reliable, 1), 258);                        % Python
verifyEqual(testCase, giCoverage(200), 0.867348, 'AbsTol', 1e-5);     % Python
verifyEqual(testCase, giCoverage(300), 0.976684, 'AbsTol', 1e-5);     % Python

verifyEqual(testCase, idxRows, 300:100:1900);
verifyEqual(testCase, ts(idxRows + 1), 15:31, 'AbsTol', 1e-9);

% The first reliable index agrees with the distribution's 95th percentile.
b = gi(2)/gi(1);
a = gi(1)/b;

verifyEqual(testCase, find(reliable, 1), ceil(gaminv(0.95, a, b)/0.01));

end


function testDroppedRowsAreTheBiasedOnes(testCase)
% Incidence falling by 20% per week: exact R = (1 + r*b)^a = 0.7759. The rows
% left out are 150% and 20% too high; every row kept is within 4%.

gi = testCase.TestData.mpox;

b = gi(2)/gi(1);
a = gi(1)/b;
exact = (1 - 0.2*b)^a;

t = 12:31;
I = 100*exp(-0.2*(0:19));

[Rs, ~, giCoverage] = get_Rt(t, I, 1, gi(1), gi(2));

idxRows = selectRtRows(giCoverage);

verifyEqual(testCase, Rs(100)/exact, 2.496, 'AbsTol', 0.01);     % Python: +149.6%
verifyEqual(testCase, Rs(200)/exact, 1.204, 'AbsTol', 0.01);     % Python: +20.4%

verifyLessThan(testCase, max(abs(Rs(idxRows)/exact - 1)), 0.04);
verifyEqual(testCase, Rs(300), 0.802078, 'RelTol', 1e-5);        % Python

end


function testOtherGenerationIntervals(testCase)

t = 12:31;
I = 100*exp(0.1*(0:19));

% Cholera interval: 95% coverage at elapsed time 2.99.
gi = testCase.TestData.cholera;

[~, ts, giCoverage] = get_Rt(t, I, 1, gi(1), gi(2));
[idxRows, reliable] = selectRtRows(giCoverage);

verifyEqual(testCase, find(reliable, 1), 300);                   % Python
verifyEqual(testCase, ts(idxRows + 1), 15:31, 'AbsTol', 1e-9);

% Fixed interval of 1: estimates exist, and are reliable, from elapsed time 1.
[~, ts, giCoverage] = get_Rt(t, I, 3, 1, 0);
[idxRows, reliable] = selectRtRows(giCoverage);

verifyEqual(testCase, find(reliable, 1), 100);
verifyEqual(testCase, ts(idxRows + 1), 13:31, 'AbsTol', 1e-9);

end


function testThresholdAndStepArguments(testCase)

gi = testCase.TestData.mpox;

[~, ~, giCoverage] = get_Rt(12:31, 100*exp(0.1*(0:19)), 1, gi(1), gi(2));

% Defaults are 0.95 and 100.
[idxDefault, reliableDefault] = selectRtRows(giCoverage);
[idxRows, reliable] = selectRtRows(giCoverage, 0.95, 100);

verifyEqual(testCase, idxRows, idxDefault);
verifyEqual(testCase, reliable, reliableDefault);

% A stricter threshold drops one more week.
idxRows = selectRtRows(giCoverage, 0.99);

verifyEqual(testCase, idxRows, 400:100:1900);

% Two rows per unit of time.
idxRows = selectRtRows(giCoverage, 0.95, 50);

verifyEqual(testCase, idxRows, 300:50:1900);

end


function testLastTimeIsAlwaysTheLastRow(testCase)
% A curve whose length is not a whole number of time units still ends with
% its last time.

giCoverage = linspace(0, 1, 250)';

[idxRows, reliable] = selectRtRows(giCoverage);

verifyEqual(testCase, find(reliable, 1), 238);
verifyEqual(testCase, idxRows, 250);

giCoverage = [linspace(0, 0.99, 100)'; ones(150, 1)];

idxRows = selectRtRows(giCoverage);

verifyEqual(testCase, idxRows, [100 200 250]);

end


function testWindowTooShortKeepsLastRowOnly(testCase)
% Generation interval of mean 10 and SD 5 on a window of 5 time units: coverage
% ends at 0.14. The table keeps the last time, flagged by its coverage.

[~, ~, giCoverage] = get_Rt(0:5, exp(0.2*(0:5)), 1, 10, 25);

[idxRows, reliable] = selectRtRows(giCoverage);

verifyEqual(testCase, giCoverage(end), 0.1429, 'AbsTol', 1e-3);    % Python
verifyFalse(testCase, any(reliable));
verifyEqual(testCase, idxRows, 500);

end


function testInputChecks(testCase)

giCoverage = linspace(0, 1, 300)';

verifyError(testCase, @() selectRtRows(giCoverage, 0), 'GrowthPredict:selectRtRows:minCoverage');
verifyError(testCase, @() selectRtRows(giCoverage, 1.2), 'GrowthPredict:selectRtRows:minCoverage');
verifyError(testCase, @() selectRtRows(giCoverage, 0.95, 0), 'GrowthPredict:selectRtRows:stepsPerUnit');
verifyError(testCase, @() selectRtRows(giCoverage, 0.95, 2.5), 'GrowthPredict:selectRtRows:stepsPerUnit');

% No estimates, no rows.
[idxRows, reliable] = selectRtRows([]);

verifyEmpty(testCase, idxRows);
verifyEmpty(testCase, reliable);

% A row vector is accepted.
idxRows = selectRtRows(linspace(0, 1, 300));

verifyEqual(testCase, idxRows, 300);

end
