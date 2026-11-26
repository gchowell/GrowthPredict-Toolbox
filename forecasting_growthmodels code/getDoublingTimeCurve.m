function [tds,C0data,curve,doublingtimes]=getDoublingTimeCurve(curve,DT,print1,priorCumulative)

% GETDOUBLINGTIMECURVE  Successive doubling times of cumulative incidence.
%
%   [tds,C0data,curve,doublingtimes] = getDoublingTimeCurve(curve,DT,print1)
%   [tds,C0data,curve,doublingtimes] = getDoublingTimeCurve(curve,DT,print1,priorCumulative)
%
% Inputs
%   curve            incidence at successive time points (column)
%   DT               time between points
%   print1           1 to plot the cumulative curve and the doubling times
%   priorCumulative  cases counted before the first point of the curve
%                    (default 0). For a window that starts inside a series,
%                    pass the sum of the observations before the window.
%
% Outputs
%   tds            rows [time, cumulative count] at the start and at each doubling
%   C0data         cumulative count at the first point
%   curve          cumulative curve at the time points, [time, count]
%   doublingtimes  time taken by each successive doubling; empty when the
%                  count does not double within the curve
%
% The count that doubles is priorCumulative + cumsum(curve). Without
% priorCumulative the count restarts at the first point of the curve, so a
% window taken from the middle of an epidemic shows short 'doublings' of its
% own running total, even when incidence is falling.

dt1=0.01;

if nargin<4 || isempty(priorCumulative)
    priorCumulative=0;
end

if ~isscalar(priorCumulative) || ~isfinite(priorCumulative) || priorCumulative<0
    error('GrowthPredict:getDoublingTimeCurve:priorCumulative', ...
        'priorCumulative must be a number that is zero or positive.');
end

tds=[];
C0data=[];
doublingtimes=[];

curve=curve(:);

data=[(0:1:length(curve)-1)' curve];

if length(data(:,1))<=1
    return
end

I0=data(1,2); % initial condition


[max1,index1]=max(data(:,2));

cumdata=priorCumulative+cumsum(data(:,2));

cum1=cumdata(end);
%
%

tb=data(:,1);

tb=tb*DT;

cumulative1=priorCumulative+cumsum(data(:,2));

C0data=cumulative1(1);

XX=(tb(1):dt1:tb(end))';

cumulative1=interp1(tb,cumulative1,XX);

index1=find(cumulative1<=cum1);

XX=XX(index1);

cumulative1=cumulative1(index1);

curve=[XX cumulative1];

curve=curve((1:(1./dt1):end),:);

if print1
    subplot(1,2,1)
    line1=plot(XX(1:index1(end)),cumulative1(1:index1(end)),'-')
    xlabel('\fontsize{16}Time')
    ylabel('\fontsize{16}Cumulative cases')

end


%

x=cumulative1;
timevect=XX;

tds=[];

C0=x(1);

tds=[0 C0];

for i=1:1:25

    C0i1=2*C0;

    index1=find(x<=C0i1);

    %[timevect(index1) x(index1)]

    tds=[tds;[timevect(index1(end)) x(index1(end))]];

    C0=C0i1;

    %pause

end

tds=unique(tds,'rows');

tds=tds(1:end-1,:);


if print1
    subplot(1,2,2)

    %semilogx(diff(tds(1:end,2)),diff(tds(:,1)),code1(j,:))
    plot(tds(1:end,2),[tds(1,1);diff(tds(:,1))],'ro--')

    ylabel('\fontsize{16}Doubling time')
    xlabel('\fontsize{16}Doubling cumulative incidence')

end


if length(tds(:,1))<=1 % redo realization
    doublingtimes=[];
    return
end


maxdoublings=30;

doublingtimes=[diff(tds(:,1))];

doublingtimes=doublingtimes(1:min(length(doublingtimes),maxdoublings),:);

window1=length(doublingtimes);

doublingtimes=doublingtimes(end-(window1-1):1:end);

cumcases=tds(end,2);

meanDoublingTime=mean(doublingtimes);


