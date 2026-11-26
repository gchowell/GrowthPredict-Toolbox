
clear
close all

%Epi-1	Epi-2	Epi-3	Epi-4	Epi-5	Epi-6	Epi-7	Epi-8	Epi-9	Epi-10	Epi-11	Epi-12	Epi-13	Epi-14	Epi-15	Epi-16	Epi-17	Epi-18	Epi-19	Epi-20	Epi-21	Epi-22	Epi-23	Epi-24	Epi-25	Epi-26	Epi-27	Epi-28	Epi-29	Epi-30	Epi-31	Epi-32	Epi-33	Epi-34	Epi-35	Epi-36	Epi-37	Epi-38	Epi-39	Epi-40	Epi-41	Epi-42	Epi-43	Epi-44	Epi-45	Epi-46	Epi-47	Epi-48	Epi-49	Epi-50	Epi-51	Epi-52	Epi-1	Epi-2	Epi-3	Epi-4	Epi-5	Epi-6	Epi-7	Epi-8	Epi-9	Epi-10	Epi-11	Epi-12	Epi-13	Epi-14	Epi-15	Epi-16	Epi-17	Epi-18	Epi-19	Epi-20	Epi-21	Epi-22	Epi-23	Epi-24	Epi-25	Epi-26	Epi-27	Epi-28	Epi-29	Epi-30	Epi-31	Epi-32	Epi-33	Epi-34	Epi-35	Epi-36	Epi-37	Epi-38	Epi-39	Epi-40	Epi-41	Epi-42	Epi-43	Epi-44	Epi-45	Epi-46	Epi-47	Epi-48	Epi-49	Epi-50	Epi-51	Epi-52	Epi-1	Epi-2	Epi-3
% District
% Dedza
% Dowa
% Kasungu
% Lilongwe
% Mchinji
% Nkhotakota
% Ntcheu
% Ntchisi
% Salima
% Chitipa
% Karonga
% Likoma
% Mzimba
% Nkhata Bay
% Rumphi
% Balaka
% Blantyre
% Chikwawa
% Chiradzulu
% Machinga
% Mangochi
% Mulanje
% Mwanza
% Nsanje
% Thyolo
% Phalombe
% Zomba
% Neno
baseDir = fullfile(pwd, 'input'); % Use current directory as the base

incidence1=load(fullfile(baseDir,'cholera-district-epidemic.txt'));

for i=1:28

    % Construct the curve variable
    curve = [(0:1:length(incidence1(:,i))-1)' incidence1(:,i)];

    % Define the minimum consecutive of non-decreasing part
    minWeeks = 3;

    % Search for monotonic (non-decreasing) part
    startIndex = -1; 
    for k = 1:(length(curve) - minWeeks + 1)
        % Check if the next 4 values are in a strictly increasing order
        if  all(diff(curve(k:k + minWeeks - 1, 2)) >=0) && curve(k, 2) > 0
            startIndex = k;
            break;
        end
    end

     if (startIndex == -1)
       continue
       i
     end

    % Trim the curve to start from the detected index
    trimmedCurve = curve(startIndex:end, :);

    % Ensure the time index starts from 0
    trimmedCurve(:, 1) = trimmedCurve(:, 1) - trimmedCurve(1, 1);

    % Save the trimmed data
    outputFileName = fullfile(baseDir, 'simulated_trimmed.txt');
    try
        % Ensure the directory exists
        outputDir = './input';
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end

        save(outputFileName, 'trimmedCurve', '-ascii');
        fprintf('Trimmed data successfully saved to %s\n', outputFileName);
    catch ME
        fprintf('Failed to save data: %s\n', ME.message);
    end

    tstart1=1;
    tend1=1;
    windowsize1=minWeeks;

    Run_Fit_GrowthModels(tstart1,tend1,windowsize1)

    plotFit_ReproductionNumber(tstart1,tend1,windowsize1)

    [Rts,params,performanceC,AICcs]=plotRt(tstart1,tend1,windowsize1);

    save(strcat('./output/Rts-window-',num2str(windowsize1),'-tstart1-',num2str(tstart1),'-tend1-',num2str(tend1),'-normal-cholera-malawi-district-',num2str(i),'.mat'),'Rts','params','performanceC','AICcs','trimmedCurve','startIndex','-mat')

    %save(strcat('./output/Rts-window-',num2str(windowsize1),'-tstart1-',num2str(tstart1),'-tend1-',num2str(tend1),'-normal-cholera-malawi-district-',num2str(i),'.mat'),'Rts','params','performanceC','AICcs','curve','-mat')

end



%%
clear
close all

windowsize1=6;
tstart1=1;
tend1=80;

cc1=1;

Rini=[];

for i=1:1:28

    figure

    %subplot(3,3,cc1)

    load(strcat('./output/Rts-window-',num2str(windowsize1),'-tstart1-',num2str(tstart1),'-tend1-',num2str(tend1),'-normal-cholera-malawi-district-',num2str(i),'.mat'),'Rts','params','performanceC','AICcs','curve','-mat')

          index1=find(curve(:,2)>0);

      index1=index1(1);
      
      if length(index1)>0

          % Plot data on the left y-axis
          yyaxis left                   % Activate left y-axis
          plot(curve(:,1), curve(:,2), '-go', 'LineWidth', 1.5); % Blue line for the left y-axis
          ylabel('Cases');             % Label for the left y-axis
          grid on;                      % Enable grid
          hold on

          yyaxis right
          plot(fix(Rts(:,1)),Rts(:,2),'b-','LineWidth', 1.5)
          hold on
          plot(fix(Rts(:,1)),Rts(:,3),'b--','LineWidth', 1.5)
          plot(fix(Rts(:,1)),Rts(:,4),'b--','LineWidth', 1.5)

          % plot Rt=1
          Rtone=[curve(1,1) 1;curve(end,1) 1];

          plot(Rtone(:,1),Rtone(:,2),'k--','LineWidth', 1.5)

          xlabel('Time (days)')

          if i==1
              legend('Cases','Rt mean','Rt 95%CI LB','Rt 95%IC UB')
              %legend('Cases','Rt (true)','Rt median (parametric)','Rt median (EpiEstim)')
          end

          set(gca,'FontSize', GetAdjustedFontSize);
          set(gcf,'color','white')

          cc1=cc1+1;


          Rini=[Rini;[index1+5 Rts(index1+6 ,1:4)]]
      end

end