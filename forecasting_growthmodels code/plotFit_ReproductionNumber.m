function param_Rt=plotFit_GrowthModels(tstart1_pass,tend1_pass,windowsize1_pass)

% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

% plot model fit to epidemic data with quantified uncertainty

close all

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 % Parameter estimation method


% <============================================================================>
% <=================== Load parameter values supplied by user =================>
% <============================================================================>

[cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP,flag1_INP,model_name1_INP,fixI0_INP,windowsize1_INP,tstart1_INP,tend1_INP]=options_fit;

[type_GId1_INP,mean_GI1_INP,var_GI1_INP]=options_Rt;


% <=======================================================================================>
% <========================== Reproduction number number parameters =======================>
% <========================================================================================>

type_GId1=type_GId1_INP; % type of GI distribution 1=Gamma, 2=Exponential, 3=Delta

mean_GI1=mean_GI1_INP;  % mean of the generation interval distribution

var_GI1=var_GI1_INP; % variance of the generation interval distribution


% <============================================================================>
% <================================ Datasets properties ==============================>
% <============================================================================>

cadfilename1=cadfilename1_INP;

caddisease=caddisease_INP;

datatype=datatype_INP;

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

%method1=0; % Type of estimation method: 0 = LSQ

d=1;

dist1=dist1_INP; %Define dist1 which is the type of error structure:

% LSQ=0,
% MLE Poisson=1,
% Pearson chi-squard=2,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

% % Define dist1 which is the type of error structure:
% switch method1
%
%     case 0
%
%         dist1=0; % Normnal distribution to model error structure
%
%         %dist1=2; % error structure type (Poisson=1; NB=2)
%
%         %factor1=1; % scaling factor for VAR=factor1*mean
%
%
%     case 3
%         dist1=3; % VAR=mean+alpha*mean;
%
%     case 4
%         dist1=4; % VAR=mean+alpha*mean^2;
%
%     case 5
%         dist1=5; % VAR=mean+alpha*mean^d;
%
% end

numstartpoints=numstartpoints_INP; % Number of initial guesses for optimization procedure using MultiStart

M=M_INP; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <============================== Growth model =====================================>
% <==============================================================================>

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

flag1=flag1_INP; % Sequence of subepidemic growth models considered in epidemic trajectory

model_name1=model_name1_INP;

% <==================================================================================>
% <=============================== other parameters=======================================>
% <==================================================================================>

fixI0=fixI0_INP; % 0=Estimate the initial number of cases; 1 = Fix the initial number of cases according to the first data point in the time series


% <==============================================================================>
% <======================== Load epiemic data ========================================>
% <==============================================================================>

data=load(strcat('./input/',cadfilename1,'.txt'));

if strcmp('CUMULATIVE',upper(cadfilename1(1:10)))==1

    data(:,2)=[data(1,2);diff(data(:,2))]; % Incidence curve

end

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

if exist('tstart1_pass','var')==1 & isempty(tstart1_pass)==0

    tstart1=tstart1_pass;

else
    tstart1=tstart1_INP;

end

if exist('tend1_pass','var')==1 & isempty(tend1_pass)==0

    tend1=tend1_pass;
else
    tend1=tend1_INP;

end

if exist('windowsize1_pass','var')==1 & isempty(windowsize1_pass)==0

    windowsize1=windowsize1_pass;
else
    windowsize1=windowsize1_INP;
end


% <=========================================================================================>
% <================================ Load short-term forecast results ==================================>
% <=========================================================================================>

factors=factor(length(tstart1:1:tend1));

if length(factors)==1
    rows=factors;
    cols=1;

elseif length(factors)==3
    rows=factors(1)*factors(2);
    cols=factors(3);
else
    rows=factors(1);
    cols=factors(2);
end

cc1=1;

param_Rt=[];

for i=tstart1:1:tend1 %rolling window analysis

    load(strcat('./output/Forecast-growthModel-',cadfilename1,'-flag1-',num2str(flag1(1)),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-0.mat'))


    % <======================================================================================>
    % <======================= Plot parameter distributions and model fit and forecast ========================>
    % <======================================================================================>

    LB1=quantile(forecast_model12',0.025)';
    LB1=(LB1>=0).*LB1;

    UB1=quantile(forecast_model12',0.975)';
    UB1=(UB1>=0).*UB1;

    median1=median(forecast_model12,2);

    if 1

        figure(101+i)
        subplot(2,4,1)
        hist(Phatss_model1(:,1))
        hold on

        line2=[param_r(1,2) 10;param_r(1,3) 10];
        line1=plot(line2(:,1),line2(:,2),'r--')
        set(line1,'LineWidth',2)

        xlabel('r')
        ylabel('Frequency')

        title(cad1)

        set(gca,'FontSize', 24);
        set(gcf,'color','white')

        subplot(2,4,2)
        hist(Phatss_model1(:,2))
        hold on

        line2=[param_p(1,2) 10;param_p(1,3) 10];
        line1=plot(line2(:,1),line2(:,2),'r--')
        set(line1,'LineWidth',2)

        xlabel('p')
        ylabel('Frequency')

        title(cad2)

        set(gca,'FontSize', 24);
        set(gcf,'color','white')

        subplot(2,4,3)

        hist(Phatss_model1(:,3))
        hold on

        line2=[param_a(1,2) 10;param_a(1,3) 10];
        line1=plot(line2(:,1),line2(:,2),'r--')
        set(line1,'LineWidth',2)

        xlabel('a')
        ylabel('Frequency')
        title(cad3)

        set(gca,'FontSize', 24);
        set(gcf,'color','white')


        subplot(2,4,4)
        hist(Phatss_model1(:,4))
        hold on

        line2=[param_K(1,2) 10;param_K(1,3) 10];
        line1=plot(line2(:,1),line2(:,2),'r--')
        set(line1,'LineWidth',2)

        xlabel('K')
        ylabel('Frequency')

        title(cad4)

        set(gca,'FontSize', 24);
        set(gcf,'color','white')


        % <========================================================================================>
        % <================================ Plot model fit and forecast ======================================>
        % <========================================================================================>

        subplot(2,4,[5 6 7 8])

        plot(timevect2,forecast_model12,'c')
        hold on

        % plot 95% PI

        line1=plot(timevect2,median1,'r-')
        set(line1,'LineWidth',2)
        hold on

        line1=plot(timevect2,LB1,'r--')
        set(line1,'LineWidth',2)
        hold on

        line1=plot(timevect2,UB1,'r--')
        set(line1,'LineWidth',2)

        % plot model fit

        color1=gray(8);
        line1=plot(timevect1,fit_model1,'color',color1(6,:))
        set(line1,'LineWidth',1)

        line1=plot(timevect2,median1,'r-')
        set(line1,'LineWidth',2)
        hold on


        % plot the data

        line1=plot(timevect_all,data_all,'bo')
        set(line1,'LineWidth',2)

        line2=[timevect1(end) 0;timevect1(end) max(quantile(forecast_model12',0.975))*1.5];

        if forecastingperiod>0
            line1=plot(line2(:,1),line2(:,2),'k--')
            set(line1,'LineWidth',2)
        end

        axis([timevect1(1) timevect2(end) 0 max(quantile(forecast_model12',0.975))*1.5])

        xlabel('Time')
        ylabel(strcat(caddisease,{' '},datatype))

        set(gca,'FontSize',24)
        set(gcf,'color','white')

        title(model_name1)

    end


    fitdata=[timevect1 data_all(1:length(timevect1)) median1 LB1 UB1];

    T = array2table(fitdata);
    T.Properties.VariableNames(1:5) = {'time','data','median','LB','UB'};
    writetable(T,strcat('./output/Fit-flag1-',num2str(flag1),'-tstart-',num2str(i),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))


    % Estimate effective reproduction number, R_t
    Rss=[];

    for realization1=1:M
        [Rs,ts]=get_Rt(timevect2,fit_model1(:,realization1),type_GId1,mean_GI1,var_GI1);

        Rss=[Rss Rs];
    end

    % <=========================================================================================>
    % <================================ Plot R_t ===============================================>
    % <=========================================================================================>

    Rtmedian=quantile(Rss',0.5);
    RtLB=quantile(Rss',0.025);
    RtUB=quantile(Rss',0.975);

    figure(500+cc1)

    subplot(2,1,1)

    plot(timevect2,forecast_model12,'c')
    hold on

    % plot 95% PI

    LB1=quantile(forecast_model12',0.025)';
    LB1=(LB1>=0).*LB1;

    UB1=quantile(forecast_model12',0.975)';
    UB1=(UB1>=0).*UB1;

    median1=median(forecast_model12,2);

    line1=plot(timevect2,median1,'r-')
    set(line1,'LineWidth',2)

    hold on
    line1=plot(timevect2,LB1,'r--')
    set(line1,'LineWidth',2)

    line1=plot(timevect2,UB1,'r--')
    set(line1,'LineWidth',2)

    % plot model fit

    color1=gray(8);
    line1=plot(timevect1,fit_model1,'color',color1(6,:))
    set(line1,'LineWidth',1)

    % plot the data

    line1=plot(timevect_all,data_all,'bo')
    set(line1,'LineWidth',2)

    line2=[timevect1(end) 0;timevect1(end) max(quantile(forecast_model12',0.975))*1.5];

    if forecastingperiod>0
        line1=plot(line2(:,1),line2(:,2),'k--')
        set(line1,'LineWidth',2)
    end

    axis([timevect1(1) timevect2(end) 0 max(quantile(forecast_model12',0.975))*1.5])

    %xlabel('Time (days)')
    %ylabel(strcat(caddisease,{' '},datatype))

    set(gca,'FontSize',24)
    set(gcf,'color','white')

    xlabel('Time')
    ylabel(strcat(caddisease,{' '},datatype))

    subplot(2,1,2)
    plot(ts(2:end),Rss,'c-')
    hold on

    line1=plot(ts(2:end),Rtmedian,'r-')
    set(line1,'LineWidth',2)
    hold on

    line1=plot(ts(2:end),RtLB,'r--')
    set(line1,'LineWidth',2)
    line1=plot(ts(2:end),RtUB,'r--')
    set(line1,'LineWidth',2)

    line2=[0 1;ts(end) 1];
    line1=plot(line2(:,1),line2(:,2),'k--')
    set(line1,'LineWidth',2)

    axis([timevect1(1) timevect2(end) 0 quantile(Rss(end,:),0.975)+4])
    line2=[timevect1(end) 0;timevect1(end) quantile(Rss(end,:),0.975)+4];
    line1=plot(line2(:,1),line2(:,2),'k--')
    set(line1,'LineWidth',2)

    xlabel('Time')

    ylabel('R_t')

    set(gca,'FontSize',24)
    set(gcf,'color','white')

    cc1=cc1+1;

    % <=============================================================================================>
    % <============================== Save file with R_t estimates ======================================>
    % <=============================================================================================>

    Rtdata=[ts(2:100:end)' Rtmedian(2:100:end)' RtLB(2:100:end)' RtUB(2:100:end)'];

    T = array2table(Rtdata);
    T.Properties.VariableNames(1:4) = {'time','Rt median','Rt 95%CI LB','Rt 95% CI UB'};
    writetable(T,strcat('./output/Rt-flag1-',num2str(flag1),'-tstart-',num2str(i),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))


    % <=========================================================================================>
    % <================================ Most recent R estimate =======================>
    % <=========================================================================================>

    param_Rt=[param_Rt;[ts(end) median(Rss(end,:)) quantile(Rss(end,:),0.025) quantile(Rss(end,:),0.975) std(Rss(end,:))]];

end

% <==================================================================================================>
% <====================== Plot temporal variation of parameters from rolling window analysis ============================>
% <==================================================================================================>

if tend1>tstart1
    figure

    subplot(2,4,[1 2 3 4])

    plot(data(:,1),data(:,2),'ro-')
    xlabel('Time')
    ylabel('Cases')
    set(gca,'FontSize',24)
    set(gcf,'color','white')


    subplot(2,4,5)

    plot(tstart1:1:tend1,param_rs(:,1),'ro-')
    hold on
    plot(tstart1:1:tend1,param_rs(:,2),'b--')
    plot(tstart1:1:tend1,param_rs(:,3),'b--')

    %line1=plot(tstart1:1:tend1,smooth(param_rs(:,1),5),'k--')
    %set(line1,'LineWidth',3)


    ylabel('r')
    set(gca,'FontSize',24)
    set(gcf,'color','white')
    xlabel('Time')

    subplot(2,4,6)
    plot(tstart1:1:tend1,param_as(:,1),'ro-')
    hold on
    plot(tstart1:1:tend1,param_as(:,2),'b--')
    plot(tstart1:1:tend1,param_as(:,3),'b--')

    %line1=plot(tstart1:1:tend1,smooth(param_as(:,1),5),'k--')
    %set(line1,'LineWidth',3)

    ylabel('a')
    set(gca,'FontSize',24)
    set(gcf,'color','white')
    xlabel('Time')

    subplot(2,4,7)
    plot(tstart1:1:tend1,param_ps(:,1),'ro-')
    hold on
    plot(tstart1:1:tend1,param_ps(:,2),'b--')
    plot(tstart1:1:tend1,param_ps(:,3),'b--')

    %line1=plot(tstart1:1:tend1,smooth(param_ps(:,1),5),'k--')
    %set(line1,'LineWidth',3)


    ylabel('p')
    set(gca,'FontSize',24)
    set(gcf,'color','white')
    xlabel('Time')

    subplot(2,4,8)
    plot(tstart1:1:tend1,param_Ks(:,1),'ro-')
    hold on
    plot(tstart1:1:tend1,param_Ks(:,2),'b--')
    plot(tstart1:1:tend1,param_Ks(:,3),'b--')


    %line1=plot(tstart1:1:tend1,smooth(param_Ks(:,1),5),'k--')
    %set(line1,'LineWidth',3)

    ylabel('K')
    set(gca,'FontSize',24)
    set(gcf,'color','white')
    xlabel('Time')

end

% <========================================================================================>
% <========================================================================================>
% <========================== Save file with calibration performance metrics ============================>
% <========================================================================================>
% <========================================================================================>

performanceC=[(tstart1:1:tend1)' zeros(length(MAECSS(:,1)),1)+windowsize1 MAECSS(:,end)  MSECSS(:,end) PICSS(:,end) WISCSS(:,end)];

T = array2table(performanceC);
T.Properties.VariableNames(1:6) = {'time','calibration_period','MAE','MSE','Coverage 95%PI','WIS'};
writetable(T,strcat('./output/performance-calibration-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))


% <=============================================================================================>
% <================= Save csv file with parameters from rolling window analysis ====================================>
% <=============================================================================================>

if method1==3 | method1==4  %save parameter alpha. VAR=mean+alpha*mean; VAR=mean+alpha*mean^2;

    rollparams=[(tstart1:1:tend1)' param_rs(:,1:end) param_ps(:,1:end) param_as(:,1:end) param_Ks(:,1:end) param_I0s(:,1:end) param_alphas(:,1:end)];
    T = array2table(rollparams);
    T.Properties.VariableNames(1:19) = {'time','r mean','r LB','r UB','p mean','p LB','p UB','a mean','a LB','a UB','K0 mean','K0 LB','K0 UB','I0 mean','I0 LB','I0 UB','alpha mean','alpha LB','alpha UB'};

elseif method1==5

    rollparams=[(tstart1:1:tend1)' param_rs(:,1:end) param_ps(:,1:end) param_as(:,1:end) param_Ks(:,1:end) param_I0s(:,1:end) param_alphas(:,1:end) param_ds(:,1:end)];
    T = array2table(rollparams);
    T.Properties.VariableNames(1:22) = {'time','r mean','r LB','r UB','p mean','p LB','p UB','a mean','a LB','a UB','K0 mean','K0 LB','K0 UB','I0 mean','I0 LB','I0 UB','alpha mean','alpha LB','alpha UB','d mean','d LB','d UB'};
else

    rollparams=[(tstart1:1:tend1)' param_rs(:,1:end) param_ps(:,1:end) param_as(:,1:end) param_Ks(:,1:end) param_I0s(:,1:end)];
    T = array2table(rollparams);
    T.Properties.VariableNames(1:16) = {'time','r mean','r LB','r UB','p mean','p LB','p UB','a mean','a LB','a UB','K0 mean','K0 LB','K0 UB','I0 mean','I0 LB','I0 UB'};
end

writetable(T,strcat('./output/parameters-rollingwindow-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

