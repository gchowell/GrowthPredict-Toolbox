function dataw=days2weeks(curve)

dataw=[];
for i=1:7:length(curve)
  if i+6<=length(curve)
        dataw=[dataw;sum(curve(i:i+6))];
    else
        dataw=[dataw;sum(curve(i:end))];
    end
 
end

