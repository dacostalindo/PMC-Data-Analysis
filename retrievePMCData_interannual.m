function [R_MAX,Beta_MAX,Beta_TOT,Z_c,Sigma_RMS]= retrievePMCData(filename)

fp = fopen(filename,'r');
fgetl(fp); fgetl(fp);fgetl(fp); fgetl(fp);fgetl(fp); fgetl(fp);fgetl(fp); 
C = fscanf(fp,'%f %f %f %f %f %f ',[6]);
fgetl(fp);fgetl(fp);
D = fscanf(fp,'%f %f %f',[3]);
%Assign data
R_MAX = C(1);
Beta_MAX =C(3);
Beta_TOT = C(6);
Z_c =D(1);
Sigma_RMS = D(2);

fclose(fp);

end

