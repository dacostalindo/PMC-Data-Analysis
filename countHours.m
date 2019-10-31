function [count1,count2,startHour1,endHour1,startHour2,endHour2]= countHours(filename)

fp = fopen(filename,'r');
fgetl(fp);
C = fscanf(fp,'%i %i %d %d %d %f %f %f %f',[9]);
startHour1 = C(6);
endHour1 = C(7);
startHour2 = C(8);
endHour2 = C(9);
count1 = endHour1 - startHour1;
count2 = endHour2 - startHour2;

if count1 < 0|| count2 < 0
   disp('error')
   fclose(fp);
   return
end

fclose(fp);

end

