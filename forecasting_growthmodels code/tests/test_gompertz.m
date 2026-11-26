function tests = test_gompertz
%TEST_GOMPERTZ  Regression tests for the time origin of the Gompertz model.
%
%   Run from the code folder:
%       results = runtests('tests/test_gompertz.m');
%       table(results)
%
%   The Gompertz model dC/dt = r*C*exp(-a*t) is the only model whose growth rate
%   depends on time. Time is measured from the first point of the curve, so the
%   result does not depend on the time label of the first observation.
%
%   Reference values marked 'Python' come from an independent Python
%   re-implementation of the fit (same model, objective, bounds, seed and start
%   points; scipy least squares in place of fmincon).
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

function testCurveMatchesClosedForm(testCase)
% C(tau) = C0*exp((r/a)*(1 - exp(-a*tau))), tau = time since the first point.

r = 0.8;
a = 0.15;
C0 = 20;

tau = (0:9)';
C = C0*exp((r/a)*(1 - exp(-a*tau)));
expected = [C(1); diff(C)];

z = [r 1 a 1 C0];

incidence = growthIncidenceCurve(z, tau, 5);

verifyEqual(testCase, incidence, expected, 'AbsTol', 0.05);
verifyEqual(testCase, incidence(2), 22.0401, 'AbsTol', 0.01);     % Python

end


function testCurveDoesNotDependOnTimeLabels(testCase)
% The same parameters give the same curve for labels 0..19, 12..31 and 100..119.
% Before the change the curve for 12..31 was the curve of r*exp(-12*a).

z = [0.8 1 0.15 1 20];

base = growthIncidenceCurve(z, (0:19)', 5);

for shift1 = [1 12 100]
    shifted = growthIncidenceCurve(z, (0:19)' + shift1, 5);
    verifyEqual(testCase, shifted, base, 'RelTol', 1e-9);
end

% The other models never depended on the labels.
zGLM = [0.9 0.8 1 5000 20];

base = growthIncidenceCurve(zGLM, (0:19)', 1);
shifted = growthIncidenceCurve(zGLM, (12:31)', 1);

verifyEqual(testCase, shifted, base, 'AbsTol', 0.05);

end


function testFitDoesNotDependOnTimeLabels(testCase)
% U.S. mpox, rows 1-20, labelled 0..19, 12..31 and 100..119. Before the change
% the last two ended with r on its upper bound of 10 and sums of squares of
% 3.94e6 and 4.11e7.

data = testCase.TestData.data;
y = data(1:20, 2);

for shift1 = [0 12 100]

    rng(1, 'twister');

    [P, fval, diagnostics] = localFit([(0:19)' + shift1, y], 1);

    verifyEqual(testCase, P(1), 1.4738, 'RelTol', 0.02);      % Python
    verifyEqual(testCase, P(3), 0.18763, 'RelTol', 0.02);     % Python
    verifyEqual(testCase, fval, 2.71708e6, 'RelTol', 0.01);   % Python
    verifyFalse(testCase, diagnostics.atLB(1) || diagnostics.atUB(1));

end

end


function testWindowsInsideTheSeries(testCase)
% Rows 8-27 (labels 7..26): the fit was the same before, but r was reported
% as 8.22 = 0.998*exp(7*0.301). Rows 13-32 (labels 12..31): r ended on its
% upper bound of 8.33 with a sum of squares of 2.81e6.

data = testCase.TestData.data;

[P, fval, diagnostics] = localFit(data(8:27, :), 1);

verifyEqual(testCase, P(1), 0.99825, 'RelTol', 0.02);     % Python
verifyEqual(testCase, P(3), 0.30114, 'RelTol', 0.02);     % Python
verifyEqual(testCase, fval, 185716, 'RelTol', 0.01);      % Python
verifyFalse(testCase, diagnostics.atLB(1) || diagnostics.atUB(1));

[P, fval, diagnostics] = localFit(data(13:32, :), 1);

verifyEqual(testCase, P(1), 0.6732, 'RelTol', 0.02);      % Python
verifyEqual(testCase, P(3), 0.3939, 'RelTol', 0.02);      % Python
verifyEqual(testCase, fval, 452472, 'RelTol', 0.01);      % Python
verifyFalse(testCase, diagnostics.atLB(1) || diagnostics.atUB(1));

end


function testTimeStepRescalesTheRates(testCase)
% With DT = 7 the same observations give the same curve and rates 7 times
% smaller. The bounds of r and a are per unit of time.

data = testCase.TestData.data;
data1 = [(0:19)', data(1:20, 2)];

[P, fval, diagnostics] = localFit(data1, 7);

verifyEqual(testCase, 7*P(1), 1.4738, 'RelTol', 0.02);
verifyEqual(testCase, 7*P(3), 0.18763, 'RelTol', 0.02);
verifyEqual(testCase, fval, 2.71708e6, 'RelTol', 0.01);
verifyEqual(testCase, 7*diagnostics.UB(1), 10.005, 'RelTol', 1e-3);    % Python

end


function testForecastContinuesTheFit(testCase)
% The forecast curve is solved from the same origin as the fitted curve.

data = testCase.TestData.data;
data1 = data(13:32, :);

params0 = initialParams(data1(:, 2), 5);

[~, ~, fitcurve, forecastcurve] = fit_model(data1, params0, 1, 10, 1, 5, 4, []);

forecastcurve = forecastcurve(:);

verifyEqual(testCase, numel(forecastcurve), 24);
verifyEqual(testCase, forecastcurve(1:20), fitcurve(:), 'AbsTol', 0.5);
verifyEqual(testCase, forecastcurve(21:24), [5.22 3.52 2.38 1.60]', 'AbsTol', 0.3);    % Python

end


% <============================================================================>
% <============================== Local helper ================================>
% <============================================================================>

function [P, fval, diagnostics] = localFit(data1, DT)

params0 = initialParams(data1(:, 2), 5);

fixI0 = 1;
numstartpoints = 10;

[P, ~, ~, ~, ~, ~, fval, ~, diagnostics] = fit_model(data1, params0, fixI0, numstartpoints, DT, 5, 0, []);

end
