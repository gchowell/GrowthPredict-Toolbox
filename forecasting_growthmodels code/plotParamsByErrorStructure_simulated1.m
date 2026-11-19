
clear
close all

params1=[];
params2=[];
params3=[];
params4=[];
params5=[];
params6=[];

calibrationperiods1=30:1:100;

for i=calibrationperiods1

    i

    table1=readmatrix(strcat('./output/parameters-rollingwindow-model_name-SEIR_plain_unrep_LSQ-fixI0-1-method-0-dist-0-tstart-1-tend-1-calibrationperiod-',num2str(i),'-horizon-30-simulated-Cases.csv'));
    params1=[params1;table1];

    table1=readmatrix(strcat('./output/parameters-rollingwindow-model_name-SEIR_plain_unrep_1-fixI0-1-method-1-dist-1-tstart-1-tend-1-calibrationperiod-',num2str(i),'-horizon-30-simulated-Cases.csv'));
    params2=[params2;table1];

    table1=readmatrix(strcat('./output/parameters-rollingwindow-model_name-SEIR_plain_unrep_3-fixI0-1-method-3-dist-3-tstart-1-tend-1-calibrationperiod-',num2str(i),'-horizon-30-simulated-Cases.csv'));
    params3=[params3;table1];

    %GLM

    table1=readmatrix(strcat('./output/parameters-rollingwindow-model_name-GLM model_1-fixI0-1-method-0-dist-0-tstart-1-tend-1-calibrationperiod-',num2str(i),'-horizon-30-simulated-Cases.csv'));
    params4=[params4;table1];

    table1=readmatrix(strcat('./output/parameters-rollingwindow-model_name-GLM model_1-fixI0-1-method-1-dist-1-tstart-1-tend-1-calibrationperiod-',num2str(i),'-horizon-30-simulated-Cases.csv'));
    params5=[params5;table1];
   
    table1=readmatrix(strcat('./output/parameters-rollingwindow-model_name-GLM model_1-fixI0-1-method-3-dist-3-tstart-1-tend-1-calibrationperiod-',num2str(i),'-horizon-30-simulated-Cases.csv'));
    params6=[params6;table1];


end


%%  parameter estimates

names = {'Normal'; 'Poisson';'Negative Binomial'};

figure
%\beta
subplot(3,3,1)
line1=plot(calibrationperiods1,params1(:,2),'bo-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params1(:,3),'b--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params1(:,4),'b--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.4,'k:')
set(line1,'LineWidth',2)

ylabel('\beta estimate')
title(names(1))
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0 5])

subplot(3,3,2)
line1=plot(calibrationperiods1,params2(:,2),'go-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params2(:,3),'g--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params2(:,4),'g--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.4,'k:')
set(line1,'LineWidth',2)

title(names(2))
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0 5])

subplot(3,3,3)
line1=plot(calibrationperiods1,params3(:,2),'ro-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params3(:,3),'r--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params3(:,4),'r--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.4,'k:')
set(line1,'LineWidth',2)

title(names(3))
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0 5])

% 
% %\kappa
% 
% subplot(3,3,4)
% line1=plot(calibrationperiods1,params1(:,5),'bo-')
% set(line1,'LineWidth',2)
% hold on
% line1=plot(calibrationperiods1,params1(:,6),'b--')
% set(line1,'LineWidth',2)
% line1=plot(calibrationperiods1,params1(:,7),'b--')
% set(line1,'LineWidth',2)
% line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.3,'k:')
% set(line1,'LineWidth',2)
% 
% ylabel('\kappa estimate')
% set(gca,'FontSize', 18);
% set(gcf,'color','white')
% axis([calibrationperiods1(1) calibrationperiods1(end) 0 4])
% 
% subplot(3,3,5)
% line1=plot(calibrationperiods1,params2(:,5),'go-')
% set(line1,'LineWidth',2)
% hold on
% line1=plot(calibrationperiods1,params2(:,6),'g--')
% set(line1,'LineWidth',2)
% line1=plot(calibrationperiods1,params2(:,7),'g--')
% set(line1,'LineWidth',2)
% line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.3,'k:')
% set(line1,'LineWidth',2)
% 
% set(gca,'FontSize', 18);
% set(gcf,'color','white')
% axis([calibrationperiods1(1) calibrationperiods1(end) 0 4])
% 
% subplot(3,3,6)
% line1=plot(calibrationperiods1,params3(:,5),'ro-')
% set(line1,'LineWidth',2)
% hold on
% line1=plot(calibrationperiods1,params3(:,6),'r--')
% set(line1,'LineWidth',2)
% line1=plot(calibrationperiods1,params3(:,7),'r--')
% set(line1,'LineWidth',2)
% line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.3,'k:')
% set(line1,'LineWidth',2)
% 
% set(gca,'FontSize', 18);
% set(gcf,'color','white')
% axis([calibrationperiods1(1) calibrationperiods1(end) 0 4])


%\rho

subplot(3,3,4)
line1=plot(calibrationperiods1,params1(:,8),'bo-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params1(:,9),'b--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params1(:,10),'b--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.8,'k:')
set(line1,'LineWidth',2)

ylabel('\rho estimate')
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0 1])

subplot(3,3,5)
line1=plot(calibrationperiods1,params2(:,8),'go-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params2(:,9),'g--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params2(:,10),'g--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.8,'k:')
set(line1,'LineWidth',2)

set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0 1])

subplot(3,3,6)
line1=plot(calibrationperiods1,params3(:,8),'ro-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params3(:,9),'r--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params3(:,10),'r--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.8,'k:')
set(line1,'LineWidth',2)

set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0 1])

%\gamma


subplot(3,3,7)
line1=plot(calibrationperiods1,params1(:,11),'bo-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params1(:,12),'b--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params1(:,13),'b--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.2,'k:')
set(line1,'LineWidth',2)

ylabel('\gamma estimate')
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0 2])
xlabel('Time')

subplot(3,3,8)
line1=plot(calibrationperiods1,params2(:,11),'go-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params2(:,12),'g--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params2(:,13),'g--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.2,'k:')
set(line1,'LineWidth',2)

set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0 2])
xlabel('Time')

subplot(3,3,9)
line1=plot(calibrationperiods1,params3(:,11),'ro-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params3(:,12),'r--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params3(:,13),'r--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.2,'k:')
set(line1,'LineWidth',2)

set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0 2])
xlabel('Time')



%% GLM model

figure

%r
subplot(3,3,1)
line1=plot(calibrationperiods1,params4(:,2),'bo-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params4(:,3),'b--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params4(:,4),'b--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.5,'k:')
set(line1,'LineWidth',2)

ylabel('r estimate')
title(names(1))
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0 0.5])

subplot(3,3,2)
line1=plot(calibrationperiods1,params5(:,2),'go-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params5(:,3),'g--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params5(:,4),'g--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.5,'k:')
set(line1,'LineWidth',2)

title(names(2))
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0 0.5])

subplot(3,3,3)
line1=plot(calibrationperiods1,params6(:,2),'ro-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params6(:,3),'r--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params6(:,4),'r--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.5,'k:')
set(line1,'LineWidth',2)

title(names(3))
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0 0.5])


%p
subplot(3,3,4)
line1=plot(calibrationperiods1,params4(:,5),'bo-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params4(:,6),'b--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params4(:,7),'b--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.5,'k:')
set(line1,'LineWidth',2)

ylabel('p estimate')
title(names(1))
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0.5 1])

subplot(3,3,5)
line1=plot(calibrationperiods1,params5(:,5),'go-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params5(:,6),'g--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params5(:,7),'g--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.5,'k:')
set(line1,'LineWidth',2)

title(names(2))
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0.5 1])

subplot(3,3,6)
line1=plot(calibrationperiods1,params6(:,5),'ro-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params6(:,6),'r--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params6(:,7),'r--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.5,'k:')
set(line1,'LineWidth',2)

title(names(3))
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 0.5 1])



%K
subplot(3,3,7)
line1=plot(calibrationperiods1,params4(:,8),'bo-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params4(:,9),'b--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params4(:,10),'b--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.5,'k:')
set(line1,'LineWidth',2)

ylabel('K estimate')
title(names(1))
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 1000 1000000])

subplot(3,3,8)
line1=plot(calibrationperiods1,params5(:,8),'go-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params5(:,9),'g--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params5(:,10),'g--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.5,'k:')
set(line1,'LineWidth',2)

title(names(2))
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 1000 1000000])

subplot(3,3,9)
line1=plot(calibrationperiods1,params6(:,8),'ro-')
set(line1,'LineWidth',2)
hold on
line1=plot(calibrationperiods1,params6(:,9),'r--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,params6(:,10),'r--')
set(line1,'LineWidth',2)
line1=plot(calibrationperiods1,zeros(length(calibrationperiods1),1)+0.5,'k:')
set(line1,'LineWidth',2)

title(names(3))
set(gca,'FontSize', 18);
set(gcf,'color','white')
axis([calibrationperiods1(1) calibrationperiods1(end) 1000 1000000])

