function tests = test_getRt
%TEST_GETRT  Regression tests for get_Rt.m.
%
%   Run from the code folder:
%       results = runtests('tests/test_getRt.m');
%       table(results)
%
%   Reference values marked 'Python' come from an independent Python
%   re-implementation of get_Rt.m (scipy gamma and exponential distributions,
%   numpy convolution). The other targets are exact: for incidence growing as
%   exp(r*t) the reproduction number is (1 + r*b)^a for a gamma generation
%   interval with shape a and scale b, and exp(r*g) for a fixed interval g.
%
%   Requires the Statistics and Machine Learning Toolbox.

tests = functiontests(localfunctions);

end


function setupOnce(testCase)

codeDir = fileparts(fileparts(mfilename('fullpath')));
addpath(codeDir);

testCase.TestData.codeDir = codeDir;

% Cholera serial interval (options_Rt_cholera.m): mean 5 days, SD 8 days, in
% weeks. Its shape is below 1, the hardest case for the discretization.
testCase.TestData.meanGI = 5/7;
testCase.TestData.varGI = (8/7)^2;

end


function teardownOnce(testCase)

rmpath(testCase.TestData.codeDir);

end


% <============================================================================>
% <================================ Tests =====================================>
% <============================================================================>

function testResultDoesNotDependOnFirstTimeLabel(testCase)
% The same curve labelled 0..19, 1..20 and 12..31 must give the same values.
% Before the change the gamma estimate at the end of the 12..31 curve was 3431
% instead of 0.81.

t = 0:19;
I = 50*exp(0.25*t)./(1 + exp(0.5*(t - 10)));      % rises, then falls

settings = {1, testCase.TestData.meanGI, testCase.TestData.varGI; ...
            2, 1, 0; ...
            3, 1, 0};

for k = 1:size(settings, 1)

    base = get_Rt(t, I, settings{k, :});

    for shift1 = [1 12 1000]

        shifted = get_Rt(t + shift1, I, settings{k, :});

        verifyEqual(testCase, isnan(shifted), isnan(base));

        ok = isfinite(base);
        verifyEqual(testCase, shifted(ok), base(ok), 'RelTol', 1e-9);

    end

end

end


function testGammaAgreesWithExactValue(testCase)
% Exponential growth: exact R = (1 + r*b)^a. Weekly samples are interpolated
% linearly, which costs a few tenths of a percent.

meanGI = testCase.TestData.meanGI;
varGI = testCase.TestData.varGI;

b = varGI/meanGI;
a = meanGI/b;

for r = [-0.2 0.1 0.3]
    for start1 = [0 1 12]

        t = start1 + (0:19);
        I = 100*exp(r*(0:19));

        Rs = get_Rt(t, I, 1, meanGI, varGI);

        verifyEqual(testCase, Rs(end), (1 + r*b)^a, 'RelTol', 0.01);

    end
end

end


function testGammaAndExponentialMatchPythonReference(testCase)
% t = 0..19, incidence 100*exp(0.3*t). A curve starting at zero gave these
% values before the change as well.

t = 0:19;
I = 100*exp(0.3*t);

[Rs, ts] = get_Rt(t, I, 1, testCase.TestData.meanGI, testCase.TestData.varGI);

verifyEqual(testCase, numel(Rs), 1900);
verifyEqual(testCase, ts(501), 5, 'AbsTol', 1e-9);
verifyEqual(testCase, Rs(500), 1.1850088598, 'RelTol', 1e-6);    % Python
verifyEqual(testCase, Rs(1000), 1.1822348870, 'RelTol', 1e-6);   % Python
verifyEqual(testCase, Rs(end), 1.1822069459, 'RelTol', 1e-6);    % Python

Rs = get_Rt(t, I, 2, 1, 0);

verifyEqual(testCase, Rs(500), 1.2945025419, 'RelTol', 1e-6);    % Python
verifyEqual(testCase, Rs(end), 1.2925563361, 'RelTol', 1e-6);    % Python

end


function testFixedIntervalIsExact(testCase)
% R = I(t)/I(t-g) = exp(r*g). Before the change the lag was g + 0.01.

t = 0:19;
I = 100*exp(0.3*t);

Rs = get_Rt(t, I, 3, 1, 0);

verifyEqual(testCase, Rs(end), exp(0.3), 'RelTol', 1e-9);
verifyEqual(testCase, find(isfinite(Rs), 1), 100);    % first estimate at elapsed time 1

Rs = get_Rt(t, I, 3, 2, 0);

verifyEqual(testCase, Rs(end), exp(0.6), 'RelTol', 1e-9);

end


function testFixedIntervalsThatUsedToFail(testCase)
% 0.14, 0.28, 0.29 and 0.57 stopped with an indexing error; 5/7 is not on the
% grid and is rounded to 0.71.

t = 0:0.01:19;                 % dense input: no interpolation error
I = exp(0.3*t);

for g = [0.14 0.28 0.29 0.57 5/7]

    [Rs, ts] = get_Rt(t, I, 3, g, 0);

    verifyEqual(testCase, numel(Rs), numel(ts) - 1);
    verifyEqual(testCase, Rs(end), exp(0.3*round(g/0.01)*0.01), 'RelTol', 1e-9);

end

end


function testFixedIntervalLimits(testCase)

t = 0:19;
I = 100*exp(0.3*t);

% Longer than the curve: no lagged incidence, so no estimate.
Rs = get_Rt(t, I, 3, 25, 0);

verifyEqual(testCase, numel(Rs), 1900);
verifyTrue(testCase, all(isnan(Rs)));

% Shorter than half a grid step.
verifyError(testCase, @() get_Rt(t, I, 3, 0.004, 0), 'GrowthPredict:get_Rt:lag');

end


function testOutputsAreAligned(testCase)

meanGI = testCase.TestData.meanGI;
varGI = testCase.TestData.varGI;

t = 12:31;
I = 100*exp(0.1*(0:19));

[Rs, ts, giCoverage] = get_Rt(t, I, 1, meanGI, varGI);

verifyEqual(testCase, ts(1), 12);
verifyEqual(testCase, ts(end), 31, 'AbsTol', 1e-9);
verifyEqual(testCase, numel(Rs), numel(ts) - 1);
verifySize(testCase, Rs, [numel(ts) - 1, 1]);
verifySize(testCase, giCoverage, size(Rs));

% Coverage is the generation-interval distribution function of elapsed time.
b = varGI/meanGI;
a = meanGI/b;

verifyGreaterThanOrEqual(testCase, min(diff(giCoverage)), 0);
verifyEqual(testCase, giCoverage(100), gamcdf(1, a, b), 'AbsTol', 1e-10);
verifyEqual(testCase, giCoverage(end), gamcdf(19, a, b), 'AbsTol', 1e-10);

end


function testInputChecks(testCase)

t = 0:19;
I = 100*exp(0.1*t);

verifyError(testCase, @() get_Rt(t, I, 4, 1, 1), 'GrowthPredict:get_Rt:type');
verifyError(testCase, @() get_Rt(t, I, 1, -1, 1), 'GrowthPredict:get_Rt:mean');
verifyError(testCase, @() get_Rt(t, I, 1, 1, 0), 'GrowthPredict:get_Rt:variance');

% A single observation gives no estimate.
[Rs, ts] = get_Rt(5, 10, 1, 1, 1);

verifyEmpty(testCase, Rs);
verifyEqual(testCase, ts, 5);

end
