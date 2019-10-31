function [year month day UT_Hour zuluday Rmax Bmax dBmax Btotal dBTotal ZcPMC RMSwidth FWHM ZBpk ZRpk factorPMC factorNoise duration]= retrievePMCData(filename)

fp = fopen(filename,'r');
fgetl(fp);
C = fscanf(fp,'%f %f %f %f %f %f %f %f %f',[9]);
year = C(2);
month = C(3);
day = C(4);
UT_Hour = C(5);
duration = C(7)-C(6);
zulumonth = [31;28;31;30;31;30;31;31;30;31;30;31];		% day numbers in each month through whole year
zuluday = sum(zulumonth(1:month-1))+day+UT_Hour/24;	% convert to zulu day 
if zuluday > 182.5
zuluday = zuluday -365;    
end
fgetl(fp);fgetl(fp);
D = fscanf(fp,'%f %f',[2]);
factorPMC = D(1);
factorNoise = D(2);
fgetl(fp);fgetl(fp);fgetl(fp);fgetl(fp);
E = fscanf(fp,'%f %f %f %f %f %f %f',[7]);
dBTotal = E(7);
Rmax = E(1);
Bmax = E(3);
dBmax = E(4);
Btotal = E(6);
ZBpk = E(5);
ZRpk = E(2);
fgetl(fp); fgetl(fp);
F = fscanf(fp,'%f %f %f',[3]);
ZcPMC = F(1);
RMSwidth = F(2);
FWHM = F(3);

fclose(fp);

end

