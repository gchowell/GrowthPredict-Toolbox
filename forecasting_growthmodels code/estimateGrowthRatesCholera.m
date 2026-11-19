close all
clear

data=load('cholera-data-Sally.txt');

for i=1:length(data(:,1))

    a=[(0:1:length(data(1,:))-1)' data(i,:)'];

    save('./input/cholera-epidemic.txt','a','-ascii')

    subplot(6,5,i)
    plot(a(:,1),a(:,2))

    hold on

    [max1,index1]=max(a(:,2))

    if index1<6
        continue
    end

    nd=5;
    np=1;

   %subplot(1,2,2)

    [index1-nd-np index1-np]

    line1=plot(a(index1-nd-np:index1-np,1),a(index1-nd-np:index1-np,2),'rx-')
    set(line1,'Linewidth',2,'markerSize',6)
    hold on

end


    %%
    Run_Fit_GrowthModels(index1-nd-np,index1-nd-np,nd)

%end








