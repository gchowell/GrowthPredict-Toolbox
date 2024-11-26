function [Rts,params,performanceC,AICcs]=plotRt(tstart1_pass,tend1_pass,windowsize1_pass)

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

% Check if fileName ends with '.txt'
if ~endsWith(cadfilename1, '.txt', 'IgnoreCase', true)
    % Append '.txt' extension if not present
    cadfilename1 = strcat(cadfilename1, '.txt');
end

% Create full file path
fullFilePath = fullfile('./input', cadfilename1);

% Check if the file exists before attempting to load
if exist(fullFilePath, 'file') == 2
    % Load the file
    data = load(fullFilePath);
else
    % Display an error message if the file is not found
    error('File "%s" not found in the specified directory.', fullFilePath);
end

if isempty(data)
    error('The dataset is empty')
end

if length(cadfilename1)>=10 & strcmp('CUMULATIVE',upper(cadfilename1(1:10)))==1

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

forecastingperiod=0;

Rts=zeros(length(tstart1:tend1),4);

cc1=1;

for i=tstart1:1:tend1 %rolling window analysis

    %strcat('./output/Rt-flag1-',num2str(flag1),'-tstart-',num2str(i),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv')

    %time	Rt median	Rt 95%CI LB	Rt 95% CI UB

    Rt=csvread(strcat('./output/Rt-flag1-',num2str(flag1),'-tstart-',num2str(i),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'),1,0);

    Rts(cc1,:)=Rt(end,1:4);

    cc1=cc1+1;

end


%time	r mean	r LB	r UB	p mean	p LB	p UB	a mean	a LB	a UB	K0 mean	K0 LB	K0 UB	I0 mean	I0 LB	I0 UB

%strcat('./output/parameters-rollingwindow-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv')

params=csvread(strcat('./output/parameters-rollingwindow-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'),1,0);

%%
% Define a separator for clarity
separator = '<=============================================================================>';

% Display options_fit.m header
disp(separator);
disp('**************************** options_fit.m ***********************************');
disp(separator);

% Display Datasets properties
disp('<================================ Datasets properties ========================>');
disp(separator);

fprintf('cadfilename1: %-30s\t(String) Name of the data file containing time-series data.\n', cadfilename1);
fprintf('caddisease:    %-30s\t(String) Name of the disease or subject related to the time-series data.\n', caddisease);
fprintf('datatype:      %-30s\t(String) Nature of the data (cases, deaths, hospitalizations, etc).\n', datatype);

disp(separator);

% Display Parameter estimation details
disp('<=========================== Parameter estimation ============================>');
disp(separator);

fprintf('method1:        %-20s\t(Method for optimization).\n', num2str(method1));
fprintf('dist1:          %-20s\t(Distribution used in the estimation).\n', num2str(dist1));
fprintf('numstartpoints: %-20s\t(Number of initial guesses for optimization).\n', num2str(numstartpoints));
fprintf('B:              %-20s\t(Number of bootstrap realizations for uncertainty).\n', num2str(M));

disp(separator);

% Display Growth model details
disp('<=========================== Growth model ====================================>');
disp(separator);

fprintf('flag1:       %-20s\t(Integer) Growth model to fit the time-series data.\n', num2str(flag1));
fprintf('model_name1: %-20s\t(String) Name of the model.\n', model_name1);
fprintf('fixI0:       %-20s\t(Boolean) Fix initial value in time-series (true/false).\n', num2str(fixI0));

disp(separator);

% Display Sliding window parameters
disp('<======================= Sliding window parameters ===========================>');
disp(separator);

fprintf('windowsize1: %-20s\t(Integer) Moving window size.\n', num2str(windowsize1));
fprintf('tstart11:    %-20s\t(Integer) Start time point for rolling window analysis.\n', num2str(tstart1));
fprintf('tend1:       %-20s\t(Integer) End time point for rolling window analysis.\n', num2str(tend1));

disp(separator);


%% plot epidemic curve, Reproduction number, and scaling of growth parameter p
 

subplot(2,1,1)

curve=data;

line1=plot(curve(tstart1:tend1+windowsize1-1,1),curve(tstart1:tend1+windowsize1-1,2),'bo-');
set(line1,'LineWidth',1.5)

set(gca,'FontSize',24)
set(gcf,'color','white')

xlabel('Time')

ylabel(strcat(caddisease,{' '},datatype))

set(gca,'FontSize',GetAdjustedFontSize)
set(gcf,'color','white')

subplot(2,1,2)

line1=plot(Rts(:,1),Rts(:,2),'r.-');
hold on
set(line1,'LineWidth',1.5)

line1=plot(Rts(:,1),Rts(:,3:4),'r--');
set(line1,'LineWidth',2)

line1=[tstart1 1;Rts(end,1) 1];

line1=plot(line1(:,1),line1(:,2),'k--');
set(line1,'LineWidth',1.5)

xlabel('Time')

axis([tstart1 Rts(end,1) min(Rts(:,3))-0.1 max(Rts(:,4))+0.1])

ylabel('Reproduction number, R_t')
set(gca,'FontSize',GetAdjustedFontSize)
set(gcf,'color','white')

% save the Rts values in a csv file
T = array2table(Rts);
T.Properties.VariableNames(1:4) = {'time','mean Rt','Rt LB','Rt UB'};
writetable(T,strcat('./output/Rts-flag1-',num2str(flag1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'));



% subplot(3,1,3)
% 
% line1=plot(tstart1+windowsize1-1:tend1+windowsize1-1,params(:,5),'g-.');
% set(line1,'LineWidth',2)
% hold on
% line1=plot(tstart1+windowsize1-1:tend1+windowsize1-1,params(:,6),'g--');
% set(line1,'LineWidth',2)
% line1=plot(tstart1+windowsize1-1:tend1+windowsize1-1,params(:,7),'g--');
% set(line1,'LineWidth',2)
% 
% ylabel('Scaling of growth parameter (p)')
% 
% set(gca,'FontSize',GetAdjustedFontSize)
% set(gcf,'color','white')
% 
% xlabel('Time')
% 
% axis([tstart1 Rts(end,1) 0 1])


%% plot parameter estimates

figure

subplot(2,2,1)

plot(tstart1:1:tend1,params(:,2),'ro-')
hold on
plot(tstart1:1:tend1,params(:,3),'b--')
plot(tstart1:1:tend1,params(:,4),'b--')

%line1=plot(tstart1:1:tend1,smooth(param_rs(:,1),5),'k--')
%set(line1,'LineWidth',3)


ylabel('r')
set(gca,'FontSize',GetAdjustedFontSize)
set(gcf,'color','white')
xlabel('Time')

subplot(2,2,2)
plot(tstart1:1:tend1,params(:,5),'ro-')
hold on
plot(tstart1:1:tend1,params(:,6),'b--')
plot(tstart1:1:tend1,params(:,7),'b--')

%line1=plot(tstart1:1:tend1,smooth(param_as(:,1),5),'k--')
%set(line1,'LineWidth',3)

ylabel('p')
set(gca,'FontSize',GetAdjustedFontSize)
set(gcf,'color','white')
xlabel('Time')

subplot(2,2,3)
plot(tstart1:1:tend1,params(:,8),'ro-')
hold on
plot(tstart1:1:tend1,params(:,9),'b--')
plot(tstart1:1:tend1,params(:,10),'b--')

%line1=plot(tstart1:1:tend1,smooth(param_ps(:,1),5),'k--')
%set(line1,'LineWidth',3)


ylabel('a')
set(gca,'FontSize',GetAdjustedFontSize)
set(gcf,'color','white')
xlabel('Time')

subplot(2,2,4)
plot(tstart1:1:tend1,params(:,11),'ro-')
hold on
plot(tstart1:1:tend1,params(:,12),'b--')
plot(tstart1:1:tend1,params(:,13),'b--')


ylabel('K')
set(gca,'FontSize',GetAdjustedFontSize)
set(gcf,'color','white')
xlabel('Time')

%% plot calibration performance metrics

%performance-calibration-flag1-1-fixI0-1-method-0-dist-0-tstart-1-tend-80-calibrationperiod-10-horizon-0-Cholera-cases.csv

performanceC=csvread(strcat('./output/performance-calibration-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'),1,0);

figure

subplot(2,2,1)
plot(tstart1+windowsize1-1:tend1+windowsize1-1,performanceC(:,3),'r-')
ylabel('MAE')
set(gca,'FontSize',GetAdjustedFontSize)
set(gcf,'color','white')
xlabel('Time')

disp(['mean MAE=',num2str(mean(performanceC(:,3)))])

%axis([tstart1 Rts(end,1) 0 max(performanceC(:,3))])

subplot(2,2,2)

plot(tstart1+windowsize1-1:tend1+windowsize1-1,performanceC(:,4),'r-')

ylabel('MSE')
set(gca,'FontSize',GetAdjustedFontSize)
set(gcf,'color','white')
xlabel('Time')

disp(['mean MSE=',num2str(mean(performanceC(:,4)))])

%axis([tstart1 Rts(end,1) 0 max(performanceC(:,4))])

subplot(2,2,3)

plot(tstart1+windowsize1-1:tend1+windowsize1-1,performanceC(:,5),'r-')

ylabel('Coverage 95% PI')
set(gca,'FontSize',GetAdjustedFontSize)
set(gcf,'color','white')
xlabel('Time')

disp(['mean coverage 95% PI=',num2str(mean(performanceC(:,5)))])

%axis([tstart1 Rts(end,1) 0 max(performanceC(:,5))])

subplot(2,2,4)

plot(tstart1+windowsize1-1:tend1+windowsize1-1,performanceC(:,6),'r-')

ylabel('WIS')
set(gca,'FontSize',GetAdjustedFontSize)
set(gcf,'color','white')
xlabel('Time')

disp(['mean WIS=',num2str(mean(performanceC(:,6)))])

%axis([tstart1 Rts(end,1) 0 max(performanceC(:,6))])

%% plot AICc

%performance-calibration-flag1-1-fixI0-1-method-0-dist-0-tstart-1-tend-80-calibrationperiod-10-horizon-0-Cholera-cases.csv

%time	AICc	AICc part1	AICc part2	numparams

AICcs=csvread(strcat('./output/AICcs-rollingwindow-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'),1,0);

figure

plot(tstart1+windowsize1-1:tend1+windowsize1-1,AICcs(:,2),'r-')

%axis([tstart1 Rts(end,1) 0 max(AICcs(:,2))])

ylabel('AICc')
set(gca,'FontSize',GetAdjustedFontSize)
set(gcf,'color','white')
xlabel('Time')

disp(['mean AICc=',num2str(mean(AICcs(:,2)))])
