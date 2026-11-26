

close all

figure

%function [curve,params]=plotGrowthModel(flag1_pass,r_pass,p_pass,a_pass,K_pass,I0_pass,windowsize1_pass)

legends=['A';'B';'C';'D'];

flag1s=[1 2 3 4];

for cc1=1:4

    subplot(2,2,cc1)

    curve=plotGrowthModel(flag1s(cc1),0.18,0.9,0.6,10000,1,200)

    text(5,max(curve)-20,strcat('\fontsize{24}',legends(cc1)))

end

