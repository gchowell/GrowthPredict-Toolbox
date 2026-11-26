
close all

Rts=[];

for i=1:80
    i

    Rt=csvread(strcat('./output/Rt-flag1-1-tstart-',num2str(i),'-fixI0-1-method-0-dist-0-calibrationperiod-10-horizon-0-Cholera-cases.csv'),1,0);

    Rts=[Rts; Rt(end,1:4)];

end

%%


subplot(2,1,1)

curve=load('./input/cholera-epidemic.txt');

line1=plot(curve(:,1),curve(:,2),'bo-')

set(gca,'FontSize',24)
set(gcf,'color','white')

ylabel('Case incidence')

subplot(2,1,2)

plot(Rts(:,1),Rts(:,2),'r-')
hold on
plot(Rts(:,1),Rts(:,3:4),'r--')

line1=[Rts(1,1) 1;Rts(end,1) 1];

plot(line1(:,1),line1(:,2),'k--')
xlabel('Time (weeks)')

axis([Rts(1,1) Rts(end,1) min(Rts(:,3))-0.1 max(Rts(:,4))+0.1])

ylabel('Effective reproduction number, R_t')

set(gca,'FontSize',24)
set(gcf,'color','white')

