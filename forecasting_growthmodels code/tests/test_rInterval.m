function tests = test_rInterval
%TEST_RINTERVAL  Regression tests for the interval of r in fit_model.m.
%
%   Run from the code folder:
%       results = runtests('tests/test_rInterval.m');
%       table(results)
%
%   The series is the bundled U.S. mpox file, GLM, least squares, fixI0 = 1.
%   Reference values come from an independent Python re-implementation of the
%   fit (same model, objective and parameter box; scipy least squares), which
%   reproduces the README quick-start estimates. Tolerances are loose enough
%   for differences between optimizers.
%
%   Requires the Optimization and Global Optimization Toolboxes.

tests = functiontests(localfunctions);

end


function setupOnce(testCase)

codeDir = fileparts(fileparts(mfilename('fullpath')));
addpath(codeDir);

testCase.TestData.codeDir = codeDir;
testCase.TestData.data = load(fullfile(codeDir, 'input', 'Most_Recent_Timeseries_US-CDC.txt'));

end


function teardownOnce(testCase)

rmpath(testCase.TestData.codeDir);

end


function setup(testCase)

global method1 %#ok<GVMIS>

testCase.TestData.method1Before = method1;
method1 = 0;                 % least squares
rng(1, 'twister');

end


function teardown(testCase)

global method1 %#ok<GVMIS>

method1 = testCase.TestData.method1Before;

end


% <============================================================================>
% <================================ Tests =====================================>
% <============================================================================>

function testQuickStartWindowIsUnchanged(testCase)
% Rows 1-20: the fit was interior before the change and must stay the same.

[P, fval, diagnostics] = localFit(testCase, 1, 20, 1, []);

verifyEqual(testCase, P(1), 1.7302, 'RelTol', 0.02);
verifyEqual(testCase, P(2), 0.8526, 'AbsTol', 0.005);
verifyEqual(testCase, P(4), 28547, 'RelTol', 0.01);
verifyEqual(testCase, fval, 3.1779e5, 'RelTol', 0.01);
verifyFalse(testCase, diagnostics.atLB(1) || diagnostics.atUB(1));

end


function testDecliningWindowIsNotLimitedByR(testCase)
% Rows 13-32 (README: 'the last 20 observations'). With r capped at 10 the
% sum of squares was 9.51e5; the data support r of about 800 and 2.32e4.

[P, fval, diagnostics] = localFit(testCase, 13, 20, 1, []);

verifyLessThan(testCase, fval, 3.0e4);
verifyGreaterThan(testCase, P(1), 10);
verifyFalse(testCase, diagnostics.atLB(1) || diagnostics.atUB(1));
verifyFalse(testCase, diagnostics.rBoundsWidened);

end


function testMidEpidemicWindowIsNotLimitedByR(testCase)
% Rows 8-27. With r capped at 10 the sum of squares was 5.49e5; reference 9.24e4.

[P, fval, diagnostics] = localFit(testCase, 8, 20, 1, []);

verifyLessThan(testCase, fval, 1.0e5);
verifyEqual(testCase, P(1), 26.97, 'RelTol', 0.10);
verifyFalse(testCase, diagnostics.atLB(1) || diagnostics.atUB(1));

end


function testIntervalScalesWithTheData(testCase)
% GGM, GLM and GRM: r in [1e-4, 10*max(y)]/DT.

data = testCase.TestData.data;
ymax = max(data(13:32, 2));

for flag1 = [0 1 2]
    [~, ~, diagnostics] = localFit(testCase, 13, 20, flag1, []);
    verifyEqual(testCase, diagnostics.LB(1), 1e-4, 'RelTol', 1e-12);
    verifyEqual(testCase, diagnostics.UB(1), 10*ymax, 'RelTol', 1e-12);
end

end


function testFixedExponentModelsKeepTheirInterval(testCase)
% Logistic model (p fixed at 1): r is a rate per unit time and its interval
% still comes from the early slope of log incidence, [slope/20, 20*slope].

data = testCase.TestData.data;
b = polyfit(data(1:7, 1), log(max(data(1:7, 2), 0.5)), 1);

[~, ~, diagnostics] = localFit(testCase, 1, 20, 3, []);

verifyEqual(testCase, diagnostics.LB(1), b(1)/20, 'RelTol', 1e-10);
verifyEqual(testCase, diagnostics.UB(1), b(1)*20, 'RelTol', 1e-10);

end


function testUserIntervalIsHonoured(testCase)
% fitopts.rBounds replaces the interval and switches off the retry.

fitopts = struct('rBounds', [1e-3 10]);

[P, fval, diagnostics] = localFit(testCase, 13, 20, 1, fitopts);

verifyEqual(testCase, diagnostics.LB(1), 1e-3, 'RelTol', 1e-12);
verifyEqual(testCase, diagnostics.UB(1), 10, 'RelTol', 1e-12);
verifyLessThanOrEqual(testCase, P(1), 10*(1 + 1e-9));
verifyTrue(testCase, diagnostics.atUB(1));
verifyFalse(testCase, diagnostics.rBoundsWidened);
verifyEqual(testCase, fval, 9.51e5, 'RelTol', 0.03);

end


function testInvalidUserIntervalIsRejected(testCase)

data = testCase.TestData.data;
data1 = data(1:20, :);
params0 = initialParams(data1(:, 2), 1);

verifyError(testCase, ...
    @() fit_model(data1, params0, 1, 2, 1, 1, 0, struct('rBounds', [5 1])), ...
    'GrowthPredict:fit_model:rBounds');

verifyError(testCase, ...
    @() fit_model(data1, params0, 1, 2, 1, 1, 0, struct('rBounds', [0 1])), ...
    'GrowthPredict:fit_model:rBounds');

end


% <============================================================================>
% <============================== Local helper ================================>
% <============================================================================>

function [P, fval, diagnostics] = localFit(testCase, firstRow, windowsize1, flag1, fitopts)

data = testCase.TestData.data;
data1 = data(firstRow:firstRow + windowsize1 - 1, :);

params0 = initialParams(data1(:, 2), flag1);

fixI0 = 1;
numstartpoints = 10;
DT = 1;

[P, ~, ~, ~, ~, ~, fval, ~, diagnostics] = fit_model(data1, params0, fixI0, numstartpoints, DT, flag1, 0, fitopts);

end
