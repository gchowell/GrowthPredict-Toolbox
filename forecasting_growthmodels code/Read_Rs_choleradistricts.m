

Rs=[]

startIndexs=[];

for i=1:28

    if i==14
        Rs=[Rs;[NaN NaN NaN]];
        continue

    end

      load(strcat('./output/Rts-window-3-tstart1-1-tend1-1-normal-cholera-malawi-district-',num2str(i),'.mat'));

      Rs=[Rs;Rts(1,2:4)];
            
      startIndexs=[startIndexs;[startIndex]];


end


hist(startIndexs)

