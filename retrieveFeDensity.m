function [count1,zuluday,UT_Hour]= retrieveFeDensity(filename)

fp = fopen(filename,'r');
fgetl(fp);
C = fscanf(fp,'%i %i %d %d %d %f %f %f %f',[9]);
day = C(4);
month = C(3);
UT_Hour = C(5);
startHour1 = C(6);
endHour1 = C(7);
startHour2 = C(8);
endHour2 = C(9);
count1 = endHour1 - startHour1;
%count2 = endHour2 - startHour2;
zulumonth = [31;28;31;30;31;30;31;31;30;31;30;31];		% day numbers in each month through whole year
zuluday = sum(zulumonth(1:month-1))+day+UT_Hour/24;	% convert to zulu day 
if zuluday > 182.5
zuluday = zuluday -365;    
end

if count1 < 0
   disp('error')
   fclose(fp);
   return
end


fclose(fp);

end

