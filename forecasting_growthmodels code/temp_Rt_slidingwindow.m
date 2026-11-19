
clear
close all

Rt_epiestim=csvread('influenza-like S-E-I-R model, R0=2, Binomial, EpiEstim Rt estimates with CrI (7-day and prior=5).csv',1,1);

%Rt_epiestim=csvread('influenza-like S-E-I-R model, R0=2, Binomial, EpiEstim Rt estimates with CrI (14-day and prior=5).csv',1,1);

Rt_true=csvread('Influenza-like SEIR R0=2 Binomial_Group 4 Gold standard Rt, E=5 (n=9767).csv',1,1);

incidence1=csvread('Influenza-like SEIR R0=2 Binomial_Group 4 daily incidence, E=5 (n=9767).csv',1,2);

% parameters
windowsize1=10;
tstart1=1;
tend1=110;

for i=1:1:10
    
    curve=[(0:1:length(incidence1(:,i))-1)' incidence1(:,i)];

    cadfilename1=strcat('./input/simulated.txts');

    save(cadfilename1,'curve','-ascii');

    Run_Fit_GrowthModels(tstart1,tend1,windowsize1)

    plotFit_ReproductionNumber(tstart1,tend1,windowsize1)

    [Rts,params,performanceC,AICcs]=plotRt(tstart1,tend1,windowsize1);

    Rt_true1=Rt_true(:,i);

    save(strcat('./output/Rts-R0-2-window-',num2str(windowsize1),'-tstart1-',num2str(tstart1),'-tend1-',num2str(tend1),'-normal-simulated-',num2str(i),'.mat'),'Rts','Rt_true1','params','performanceC','AICcs','curve','-mat')

end

%%

%Run_Forecasting_GrowthModels(1,80,10,4)

Run_Fit_GrowthModels(65,65,7)

plotFit_GrowthModels(65,65,7)
plotFit_ReproductionNumber(65,65,7)

%plotForecast_ReproductionNumber(1,80,10,4)

%plotFit_ReproductionNumber(68,69,10)

plotRt(65,65,7)

%%
clear

close all

% parameters
windowsize1=7;
tstart1=1;
tend1=110;

Rt_true=csvread('Influenza-like SEIR R0=2 Binomial_Group 4 Gold standard Rt, E=5 (n=9767).csv',1,1);

incidence1=csvread('Influenza-like SEIR R0=2 Binomial_Group 4 daily incidence, E=5 (n=9767).csv',1,2);

switch windowsize1
    case 7
      Rt_epiestim=csvread('influenza-like S-E-I-R model, R0=2, Binomial, EpiEstim Rt estimates with CrI (7-day and prior=5), I(0)=5.csv',1,1);
    case 10
      Rt_epiestim=csvread('influenza-like S-E-I-R model, R0=2, Binomial, EpiEstim Rt estimates with CrI (7-day and prior=5).csv',1,1);
    case 14
       Rt_epiestim=csvread('influenza-like S-E-I-R model, R0=2, Binomial, EpiEstim Rt estimates with CrI (14-day and prior=2), I(0)=5.csv',1,1);
end


cc1=1;
MAES=[]
MSES=[];

for i=1:10

    figure(100)

    subplot(2,2,cc1)

    %load(strcat('./output/Rts-R0-2-window-7-normal-simulated-',num2str(i),'.mat'))

    load(strcat('./output/Rts-R0-2-window-',num2str(windowsize1),'-tstart1-',num2str(tstart1),'-tend1-',num2str(tend1),'-normal-simulated-',num2str(i),'.mat'),'Rts','Rt_true1','params','performanceC','AICcs','curve','-mat')


    % Plot data on the left y-axis
    yyaxis left                   % Activate left y-axis
    plot(curve(:,1), curve(:,2), '-go', 'LineWidth', 1.5); % Blue line for the left y-axis
    ylabel('Cases');             % Label for the left y-axis
    grid on;                      % Enable grid


    % Plot data on the right y-axis
    yyaxis right                  % Activate right y-axis

    plot((0:1:length(Rt_true1)-1)', Rt_true1,'r','LineWidth', 1.5)
    hold on

    plot(fix(Rts(:,1)),Rts(:,2),'b-','LineWidth', 1.5)
    hold on
    %plot(fix(Rts(:,1)),Rts(:,3),'b--','LineWidth', 1.5)
    %plot(fix(Rts(:,1)),Rts(:,4),'b--','LineWidth', 1.5)


    %plot Rt EpiEstim
    %plot(Rts(:,1),Rt_epiestim(1:length(Rts(:,1)),(i-1)*3+1),'m--','LineWidth', 1.5)
    hold on
    plot(fix(Rts(:,1)),Rt_epiestim(1:length(Rts(:,1)),(i-1)*3+2),'m-','LineWidth', 1.5)
    hold on
    %plot(Rts(:,1),Rt_epiestim(1:length(Rts(:,1)),(i-1)*3+3),'m--','LineWidth', 1.5)

    % plot Rt=1
    Rtone=[curve(1,1) 1;curve(end,1) 1];

    plot(Rtone(:,1),Rtone(:,2),'k--','LineWidth', 1.5)

    xlabel('Time (days)')
    %legend('Cases','Rt (true)','Rt mean','Rt 95%CI LB','Rt 95%IC UB','Rt (EpiEstim)')
    legend('Cases','Rt (true)','Rt median (parametric)','Rt median (EpiEstim)')

    set(gca,'FontSize', GetAdjustedFontSize);
    set(gcf,'color','white')

    % estimate Rt errors

    Rvalues=[Rt_true1(fix(Rts(:,1))+1) Rts(:,2) Rt_epiestim(1:length(Rts(:,1)),(i-1)*3+2)];

    MAE1=mean(abs(Rvalues(:,2)-Rvalues(:,1)));
    MAE2=mean(abs(Rvalues(:,3)-Rvalues(:,1)));

    MSE1=mean((Rvalues(:,2)-Rvalues(:,1)).^2);
    MSE2=mean((Rvalues(:,3)-Rvalues(:,1)).^2);

    MAES=[MAES;[MAE1 MAE2]];

    MSES=[MSES;[MSE1 MSE2]];
    
    figure(101)

    subplot(2,2,cc1)
    plot(abs(Rvalues(:,1)-Rvalues(:,2)),'b-')
    hold on
    plot(abs(Rvalues(:,1)-Rvalues(:,3)),'r-')

    pause

    cc1=cc1+1;

end


%%
close all

incidence1=csvread('Influenza-like SEIR R0=2 Binomial_Group 4 daily incidence, E=5 (n=9767).csv',1,2);

windowsize1 = 7;

% Example: Compute CV over sliding windows and plot

for j=1:10
    % Generate example time series data (growth trend with noise)
    time_points = length(incidence1(:,j));
    time_series = incidence1(:,j);

    % Parameters for sliding window
    num_windows = time_points - windowsize1 + 1;

    % Initialize array to store CV values
    cv_values = zeros(1, num_windows);

    % Compute CV over sliding windows

    cv_time_index=[];
    for i = 1:num_windows
        window_data = time_series(i:i +  windowsize1 - 1);
        mean_val = mean(window_data);
        std_val = std(window_data, 1); % Use sample standard deviation
        if mean_val ~= 0
            cv_values(i) = (std_val / mean_val) * 100; % Compute CV
        else
            cv_values(i) = 0; % Handle zero mean case
        end

        %cv_time_index=[cv_time_index;i +  windowsize1 - 1];
        cv_time_index=[cv_time_index;i];

    end

    % Create time indices for CV values
    %cv_time_index = (ceil( windowsize1 / 2)):(time_points - floor( windowsize1 / 2));

    % Plot the results

    % Plot time series
    subplot(2, 1, 1);
    plot(1:time_points, time_series, 'b', 'LineWidth', 1.5);
    hold on
    title('Time Series Data');

    xlabel('Time');
    ylabel('Values');
    grid on;

    % Plot CV
    subplot(2, 1, 2);
    plot(cv_time_index, cv_values, 'g', 'LineWidth', 1.5);
    hold on
    title('Coefficient of Variation (CV) Over Sliding Windows');
    xlabel('Time');
    ylabel('CV (%)');
    grid on;

end
