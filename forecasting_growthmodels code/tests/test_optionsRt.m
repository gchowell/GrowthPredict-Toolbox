function tests = test_optionsRt
%TEST_OPTIONSRT  Regression tests for the generation-interval option files.
%
%   Run from the code folder:
%       results = runtests('tests/test_optionsRt.m');
%       table(results)
%
%   Reference values marked 'Python' come from an independent Python
%   re-implementation of get_Rt.m.
%
%   Requires the Statistics and Machine Learning Toolbox.

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

function testDefaultIsMpoxSerialInterval(testCase)
% U.S. case pairs by symptom onset: mean 8.5 days, SD 5.0 days, in weeks.

[type_GId1, mean_GI1, var_GI1] = options_Rt();

verifyEqual(testCase, type_GId1, 1);
verifyEqual(testCase, mean_GI1, 8.5/7, 'RelTol', 1e-12);
verifyEqual(testCase, var_GI1, (5.0/7)^2, 'RelTol', 1e-12);

end


function testMpoxFileMatchesDefault(testCase)
% options_Rt_mpox.m held mean 3.8 with variance 2*(1/3.8)^2, which do not
% describe one distribution, and a syntax error.

[typeDefault, meanDefault, varDefault] = options_Rt();
[type_GId1, mean_GI1, var_GI1] = options_Rt_mpox();

verifyEqual(testCase, type_GId1, typeDefault);
verifyEqual(testCase, mean_GI1, meanDefault, 'RelTol', 1e-12);
verifyEqual(testCase, var_GI1, varDefault, 'RelTol', 1e-12);

end


function testCholeraFileKeepsFormerDefault(testCase)
% Mean 5 days, SD 8 days, in weeks: the default of options_Rt.m until this version.

[type_GId1, mean_GI1, var_GI1] = options_Rt_cholera();

verifyEqual(testCase, type_GId1, 1);
verifyEqual(testCase, mean_GI1, 5/7, 'RelTol', 1e-12);
verifyEqual(testCase, var_GI1, (8/7)^2, 'RelTol', 1e-12);

end


function testImpliedGammaDistributions(testCase)
% Shape = mean^2/variance. The mpox interval has almost no mass within one
% day; the cholera interval, applied without its one-day shift, has 41%.

[~, m, v] = options_Rt();
b = v/m;
a = m/b;

verifyEqual(testCase, a, 2.89, 'RelTol', 1e-10);
verifyEqual(testCase, gamcdf(1/7, a, b), 0.0066, 'AbsTol', 0.0005);

[~, m, v] = options_Rt_cholera();
b = v/m;
a = m/b;

verifyEqual(testCase, a, 25/64, 'RelTol', 1e-10);
verifyEqual(testCase, gamcdf(1/7, a, b), 0.4072, 'AbsTol', 0.0005);

end


function testRtWithDefaultOnExponentialGrowth(testCase)
% t = 0..19, incidence 100*exp(0.3*t). Exact R = (1 + r*b)^a = 1.4093.

[type_GId1, mean_GI1, var_GI1] = options_Rt();

t = 0:19;
I = 100*exp(0.3*t);

Rs = get_Rt(t, I, type_GId1, mean_GI1, var_GI1);

b = var_GI1/mean_GI1;
a = mean_GI1/b;

verifyEqual(testCase, Rs(end), (1 + 0.3*b)^a, 'RelTol', 0.01);
verifyEqual(testCase, Rs(500), 1.4006215372, 'RelTol', 1e-6);    % Python
verifyEqual(testCase, Rs(end), 1.4004398175, 'RelTol', 1e-6);    % Python

% The cholera interval gives a value much closer to 1 for the same curve.
[type_GId1, mean_GI1, var_GI1] = options_Rt_cholera();

Rs = get_Rt(t, I, type_GId1, mean_GI1, var_GI1);

verifyEqual(testCase, Rs(end), 1.1822069459, 'RelTol', 1e-6);    % Python

end


function testEveryRtOptionsFileRuns(testCase)
% Each options_Rt*.m in the code folder must return a usable setting.

files = dir(fullfile(testCase.TestData.codeDir, 'options_Rt*.m'));

verifyGreaterThanOrEqual(testCase, numel(files), 3);

for k = 1:numel(files)

    [~, name1] = fileparts(files(k).name);

    optionsFunction = str2func(name1);

    [type_GId1, mean_GI1, var_GI1] = optionsFunction();

    verifyTrue(testCase, isscalar(type_GId1) && any(type_GId1 == [1 2 3]), name1);
    verifyTrue(testCase, isscalar(mean_GI1) && isfinite(mean_GI1) && mean_GI1 > 0, name1);
    verifyTrue(testCase, isscalar(var_GI1) && isfinite(var_GI1) && var_GI1 >= 0, name1);

    if type_GId1 == 1
        verifyGreaterThan(testCase, var_GI1, 0, name1);
    end

end

end
