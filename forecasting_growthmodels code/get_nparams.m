function numparams=get_nparams(method1,dist1,flag1,fixI0)
%GET_NPARAMS  Number of estimated quantities of a growth-model fit.
%
%   numparams = get_nparams(method1, dist1, flag1, fixI0)
%
%   Model parameters (by flag1) + the initial value when fixI0 = 0 + the
%   parameters of the error model. getAICc calls this function, so the count
%   is defined in one place only.
%
%   Least squares (method1 = 0) counts one parameter for the residual
%   variance whatever dist1 is: the AICc of a least-squares fit is a Gaussian
%   likelihood with an estimated variance, and dist1 only selects the noise
%   used for the bootstrap. (Before November 2026 the variance was counted
%   only for dist1 = 0, so the same fit had two different AICc values.)
%
%   dist1 is kept in the argument list for compatibility; it is not used.

%flag -1=EXP 0=GGM 1=GLM 2=GRM 3=Logistic 4=RIC 5=Gompertz

switch flag1 % model indicator

    case -1
        numparams=1;

    case 0
        numparams=2;

    case 1
        numparams=3;

    case 2
        numparams=4;

    case 3
        numparams=2;

    case 4
        numparams=3;

    case 5
        numparams=2;

    otherwise
        error('GrowthPredict:get_nparams:UnknownModel', ...
            'flag1 = %g is not a supported growth model (-1 to 5).', flag1);

end

if fixI0==0 % initial value is estimated

    numparams=numparams+1;

end

switch method1

    case 0 % least squares: residual variance

        numparams=numparams+1;

    case 1 % Poisson: no extra parameter

    case {3,4} % negative binomial: alpha

        numparams=numparams+1;

    case 5 % negative binomial: alpha and d

        numparams=numparams+2;

    otherwise
        error('GrowthPredict:get_nparams:UnsupportedMethod', ...
            'method1 = %g is not supported; use 0, 1, 3, 4 or 5.', method1);

end

end
