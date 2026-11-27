function tests = test_adjustLastWindowStart
%TEST_ADJUSTLASTWINDOWSTART  Regression tests for adjustLastWindowStart.m.
%
%   Run from the code folder:
%       results = runtests('tests/test_adjustLastWindowStart.m');
%       table(results)
%
%   adjustLastWindowStart.m gives the plotting functions the last window start
%   that the runners use when the requested one does not fit in the series.
%   The expected values are those of the rule in Run_Fit_GrowthModels.m and
%   Run_Forecasting_GrowthModels.m, written out in localRunnerRule below.
%
%   No toolbox is required.

tests = functiontests(localfunctions);

end


function setupOnce(testCase)

codeDir = fileparts(fileparts(mfilename('fullpath')));
addpath(codeDir);

testCase.TestData.codeDir = codeDir;

data = load(fullfile(codeDir, 'input', 'Most_Recent_Timeseries_US-CDC.txt'));
testCase.TestData.nMpox = size(data, 1);

end


function teardownOnce(testCase)

rmpath(testCase.TestData.codeDir);

end


% <============================================================================>
% <================================ Tests =====================================>
% <============================================================================>

function testWindowThatFitsIsKept(testCase)
% U.S. mpox series, 32 rows: the README examples are not adjusted.

n = testCase.TestData.nMpox;

verifyEqual(testCase, n, 32);

[tend1, adjusted] = adjustLastWindowStart(1, 20, n);       % Run_Fit_GrowthModels(1,1,20)
verifyEqual(testCase, tend1, 1);
verifyFalse(testCase, adjusted);

[tend1, adjusted] = adjustLastWindowStart(6, 20, n);       % six rolling windows
verifyEqual(testCase, tend1, 6);
verifyFalse(testCase, adjusted);

[tend1, adjusted] = adjustLastWindowStart(13, 20, n);      % rows 13-32, the last window that fits
verifyEqual(testCase, tend1, 13);
verifyFalse(testCase, adjusted);

end


function testWindowPastTheEndIsMovedBack(testCase)
% Run_Fit_GrowthModels(1,20,20) saves its files with tend-13; the plotting
% functions looked for tend-20.

n = testCase.TestData.nMpox;

[tend1, adjusted] = adjustLastWindowStart(14, 20, n);
verifyEqual(testCase, tend1, 13);
verifyTrue(testCase, adjusted);

[tend1, adjusted] = adjustLastWindowStart(20, 20, n);
verifyEqual(testCase, tend1, 13);
verifyTrue(testCase, adjusted);

% options_fit_mpox.m as shipped asks for row 37 of 32: the last window that
% fits starts at row 13, below the first window start of 37.
[tend1, adjusted] = adjustLastWindowStart(37, 20, n);
verifyEqual(testCase, tend1, 13);
verifyTrue(testCase, adjusted);

% options_fit_cholera.m: 80 windows of 10 in a series of 107 rows fit.
[tend1, adjusted] = adjustLastWindowStart(80, 10, 107);
verifyEqual(testCase, tend1, 80);
verifyFalse(testCase, adjusted);

end


function testSameAsTheRunnerRule(testCase)
% Every combination of a series of 5 to 60 rows, a window of 2 to 40 and a
% requested start of 1 to 70.

for n = [5 16 32 50 60]
    for w = [2 7 10 20 30 40]
        for tend1 = 1:70

            [expected, expectedAdjusted] = localRunnerRule(tend1, w, n);

            [actual, adjusted] = adjustLastWindowStart(tend1, w, n);

            label = sprintf('tend1 = %d, window %d, %d rows', tend1, w, n);

            verifyEqual(testCase, actual, expected, label);
            verifyEqual(testCase, adjusted, expectedAdjusted, label);

            if adjusted
                % the adjusted window ends exactly at the last row
                verifyEqual(testCase, actual + w - 1, n, label);
            else
                verifyLessThanOrEqual(testCase, actual + w - 1, n, label);
            end

        end
    end
end

end


function testWindowLongerThanTheSeries(testCase)
% The result is below 1: no window fits, and the caller has nothing to process.

[tend1, adjusted] = adjustLastWindowStart(1, 40, 32);

verifyEqual(testCase, tend1, -7);
verifyTrue(testCase, adjusted);
verifyLessThan(testCase, tend1, 1);

end


function testInputChecks(testCase)

verifyError(testCase, @() adjustLastWindowStart([1 2], 20, 32), 'GrowthPredict:adjustLastWindowStart:input');
verifyError(testCase, @() adjustLastWindowStart(1, [], 32), 'GrowthPredict:adjustLastWindowStart:input');
verifyError(testCase, @() adjustLastWindowStart(1, 20, 'n'), 'GrowthPredict:adjustLastWindowStart:input');

end


% <============================================================================>
% <============================== Local helper ================================>
% <============================================================================>

function [tend1, adjusted] = localRunnerRule(tend1, windowsize1, nRows)
% The lines of Run_Fit_GrowthModels.m and Run_Forecasting_GrowthModels.m.

adjusted = false;

if (tend1+windowsize1-1) > nRows

    tend1 = nRows-windowsize1+1;
    adjusted = true;

end

end
