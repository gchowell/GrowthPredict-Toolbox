function [Rs,ts,giCoverage]=get_Rt(timevect,incidence1,type_GId1,mean1,var1)

% GET_RT  Reproduction number from an incidence curve with the renewal equation.
%
%   [Rs,ts] = get_Rt(timevect,incidence1,type_GId1,mean1,var1)
%   [Rs,ts,giCoverage] = get_Rt(...)
%
% Inputs
%   timevect    time labels of the incidence curve (increasing)
%   incidence1  incidence at those times
%   type_GId1   generation-interval distribution: 1 = gamma, 2 = exponential,
%               3 = fixed interval
%   mean1       mean of the generation interval, in the units of timevect
%               (for type 3: the fixed interval, rounded to the nearest 0.01)
%   var1        variance of the generation interval (gamma only)
%
% Outputs
%   ts          grid timevect(1):0.01:timevect(end)
%   Rs          column with numel(ts)-1 values; Rs(k) is the estimate at ts(k+1),
%               so plot it against ts(2:end). NaN where the denominator is zero.
%   giCoverage  column like Rs: fraction of the generation-interval distribution
%               that the incidence history available at ts(k+1) covers.
%
% Method
%   Incidence is interpolated linearly on the grid, and
%
%       R(t) = I(t) / sum_s I(t-s) w(s),
%
%   where w(s) is the probability that the generation interval lies in the grid
%   step ending at s. The interval is measured from the start of the curve:
%   s is elapsed time, so the result does not depend on the value of timevect(1).
%   For a fixed interval g the estimate is I(t)/I(t-g).
%
% Note
%   Only the incidence inside the curve is used. Until the history covers most
%   of the generation-interval distribution the denominator is too small and
%   the estimate too high; giCoverage shows when that is the case.
%
% History
%   Before this version the gamma and exponential distributions were evaluated
%   at the time labels themselves, which is right only when timevect(1) = 0; a
%   curve starting at any other time gave values that were far too large. The
%   fixed interval was placed one grid step late (g + 0.01), and intervals such
%   as 0.29 or 0.57 stopped with an indexing error because mean1/0.01 is not an
%   integer in floating point.

dt1=0.01;

if ~isscalar(type_GId1) || ~any(type_GId1==[1 2 3])
    error('GrowthPredict:get_Rt:type', ...
        'type_GId1 must be 1 (gamma), 2 (exponential) or 3 (fixed interval).');
end

if ~isscalar(mean1) || ~isfinite(mean1) || mean1<=0
    error('GrowthPredict:get_Rt:mean', ...
        'The mean of the generation interval must be a positive number.');
end

if type_GId1==1 && (nargin<5 || ~isscalar(var1) || ~isfinite(var1) || var1<=0)
    error('GrowthPredict:get_Rt:variance', ...
        'The variance of the gamma generation interval must be a positive number.');
end

timevect=timevect(:)';
incidence1=incidence1(:)';

Rs=[];
giCoverage=[];

if numel(timevect)<2
    ts=timevect;
    return
end

timevect2=timevect(1):dt1:timevect(end);

incidence2=interp1(timevect,incidence1,timevect2);

ts=timevect2;

n=numel(incidence2);

if n<2
    return
end

% Elapsed time since the start of the curve: the lag axis of the generation interval.
tau=(0:n-1)*dt1;

switch type_GId1

    case 1 % gamma

        b=var1/mean1;
        a=mean1/b;

        prob1=diff(gamcdf(tau,a,b));

    case 2 % exponential

        prob1=diff(expcdf(tau,mean1));

    case 3 % fixed interval

        lag1=round(mean1/dt1); % number of grid steps

        if lag1<1
            error('GrowthPredict:get_Rt:lag', ...
                'A fixed generation interval must be at least %g (the grid step).', dt1);
        end

        prob1=zeros(1,n-1);

        if lag1<=n-1
            prob1(lag1)=1;
        end

end

% Renewal denominator by full convolution: for the estimate at grid index i,
% c(i-1) = sum_{j=1}^{i-1} I(i-j) * prob1(j), where prob1(j) belongs to lag j*dt1.

I=incidence2(:);
w=prob1(:);

c=conv(I,w,'full');
denominator=c(1:n-1);

Rs=NaN(n-1,1);

currentIncidence=I(2:n);

valid=denominator>0;

Rs(valid)=currentIncidence(valid)./denominator(valid);

giCoverage=cumsum(w);
