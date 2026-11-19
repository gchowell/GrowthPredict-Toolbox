function curves=AddErrorStructure(yi,M,dist1,factor1,d)
%ADDERRORSTRUCTURE  Draw M noisy incidence curves around a model curve.
%
%   curves = AddErrorStructure(yi, M, dist1, factor1, d)
%
%   yi       CUMULATIVE model curve; the mean of observation t is yi(t)-yi(t-1)
%   M        number of realizations (columns of curves)
%   dist1    0 normal, sd = factor1
%            1 Poisson
%            2 negative binomial, VAR = factor1*mean
%            3 negative binomial, VAR = mean + factor1*mean
%            4 negative binomial, VAR = mean + factor1*mean^2
%            5 negative binomial, VAR = mean + factor1*mean^d
%
%   The first point of every realization is yi(1), as before.
%
%   Revised November 2026 (following QuantDiffForecast): a zero model mean is
%   replaced by 0.001 for the negative binomial cases (it gave p = 0/0 and NaN
%   draws), a negative binomial with no over-dispersion (VAR <= mean) is drawn
%   as Poisson (it gave an infinite size parameter), the output is
%   preallocated, and non-finite draws stop with an error instead of entering
%   the bootstrap data.

yi=yi(:);

n=numel(yi);

curves=zeros(n,M);

means=abs(diff(yi));   % model mean of observations 2..n

for realization=1:M

    yirData=zeros(n,1);

    yirData(1)=yi(1);

    switch dist1

        case 0  % normal

            for t=2:n
                yirData(t,1)=normrnd(yi(t)-yi(t-1),factor1);
            end

        case 1  % Poisson

            for t=2:n
                yirData(t,1)=poissrnd(means(t-1),1,1);
            end

        case 2  % negative binomial, VAR = factor1*mean

            for t=2:n
                mean1=max(means(t-1),0.001);
                yirData(t,1)=localNegBin(mean1,mean1*factor1);
            end

        case 3  % negative binomial, VAR = mean + alpha*mean

            for t=2:n
                mean1=max(means(t-1),0.001);
                yirData(t,1)=localNegBin(mean1,mean1+mean1*factor1);
            end

        case 4  % negative binomial, VAR = mean + alpha*mean^2

            for t=2:n
                mean1=max(means(t-1),0.001);
                yirData(t,1)=localNegBin(mean1,mean1+factor1*mean1^2);
            end

        case 5  % negative binomial, VAR = mean + alpha*mean^d

            for t=2:n
                mean1=max(means(t-1),0.001);
                yirData(t,1)=localNegBin(mean1,mean1+factor1*mean1^d);
            end

        otherwise
            error('GrowthPredict:AddErrorStructure:UnsupportedDist', ...
                'dist1 = %g is not supported; use 0 to 5.', dist1);

    end

    curves(:,realization)=yirData;

end

if any(~isfinite(curves(:)))
    error('GrowthPredict:AddErrorStructure:InvalidDraws', ...
        ['Non-finite draws were generated (dist1 = %g, factor1 = %g). ' ...
        'Check the fitted curve and the noise parameters.'], dist1, factor1);
end

end


function x=localNegBin(mean1,var1)
% One negative binomial draw with the given mean and variance. Without
% over-dispersion (VAR <= mean, up to rounding) the limit is Poisson.
if ~(var1 > mean1*(1+1e-10))
    x=poissrnd(mean1,1,1);
else
    p1=mean1/var1;
    r1=mean1*p1/(1-p1);
    x=nbinrnd(r1,p1,1,1);
end
end
