function [Rs,ts]=get_Rt(timevect,incidence1,type_GId1,mean1,var1)

% type of GI distribution 1=Gamma, 2=Exponential, 3=Delta

dt1 = 0.01;
timevect2 = timevect(1):dt1:timevect(end);
incidence1 = interp1(timevect, incidence1, timevect2, 'linear');
timevect = timevect2;
timevect3 = timevect - min(timevect);
ts = timevect;

% Compute generation interval distribution
if type_GId1 == 1
    b = var1 / mean1;
    a = mean1 / b;
    prob1 = gamcdf(timevect3, a, b);
    prob1 = diff(prob1)';
elseif type_GId1 == 2
    prob1 = expcdf(timevect3, mean1);
    prob1 = diff(prob1)';
elseif type_GId1 == 3
    prob1 = zeros(length(timevect), 1);
    index = round(mean1 / dt1) + 1;
    if index <= length(prob1)
        prob1(index) = 1;
    else
        error('mean1/dt1 exceeds the range of the time vector.');
    end
else
    error('Invalid type_GId1. Must be 1 (Gamma), 2 (Exponential), or 3 (Delta).');
end

% Normalize prob1
prob1 = prob1 / sum(prob1);

Rs=[];

for i=2:length(incidence1)

    sum1=0;

    for j=1:i-1
        sum1=sum1+(incidence1(i-j))*prob1(j,1);
        %'prob1='
        %prob1(j,1)

    end

    if sum1>0
        R1=incidence1(i)/sum1;

    else
        R1=NaN;

    end

    Rs=[Rs;R1];

end

