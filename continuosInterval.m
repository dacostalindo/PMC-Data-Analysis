function [updatedList] = continuosInterval( smdsmdCnt, time0,timeRaw,nameofprofile)
%Purpose: Find continuos time intervals
logicVector = ~isnan(smdsmdCnt(:,1))';
flag = 0;
Startcount = 1;
Stopcount = 1;
start = NaN;
stop = NaN;
for i = 1:length(logicVector)

    if flag == 0 && logicVector(i) == 1
        start(Startcount) = i;    
        flag = 1;   
        Startcount = Startcount+1;
    end
    
    
    if flag == 1 && logicVector(i) == 0
        stop(Stopcount) = i-1;    
        flag = 0;   
        Stopcount = Stopcount+1; 
    end
    
    if flag == 1 && i == length(logicVector)
        stop(Stopcount) = i;    
        flag = 0;   
        Stopcount = Stopcount+1; 
    end

end
updatedList = [];
if isnan(start(1))== 1
   return 
end

for i = 1:length(start)
  interval =   find(timeRaw >= time0(start(i)) & timeRaw <= time0(stop(i)));
  updatedList = [updatedList; nameofprofile(interval,:)];
   
end
end

