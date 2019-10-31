function [ATPND] = msise00cal(year, month, day, time, latitude, longitude, fluxindex)

%% Function to calculate MSISE-00 model output using "atmosnrlmsise00" 
% for instantaneous MSISE-00 output at (year, DOY, UTseconds),
% Read in observed solar flux and geomagnetic index (F107A, F107, APH) as
% input to "atmosnrlmsise00" function
% Output includes temperature (T) and density (rho)
% Pressure (P) can be calculated from temperature and density
% P = rho(6)(kg/m^3) *287(J/K/kg) * T(K) = J/m^3 = Pascal
%% Calculate Day of Year, UT seconds, and UT hour, and assign altitude
zulumonth = [31;28;31;30;31;30;31;31;30;31;30;31];		% day numbers in each month through whole year
DOY = sum(zulumonth(1:month-1))+day;	% convert to zulu day (1-365)
UTseconds = 3600*time;   % UT 1 hour
UThour = floor(time);
altitude = (0:0.5:250).*1e3;  % unit: meters
%% Calculate MSISE-00 using default solar flux and geomagnetic index
% f107a = 150;      % default value of average F107 (81-day average centered on DOY)
% f107 = 150;       % default value of daily F107 (for previous day)
% aph = 4;          % default value of geomagnetic index
if (fluxindex == 1) % Use default values of solar flux and geomagnetic index
    [T rho] = atmosnrlmsise00(altitude, latitude, longitude, year, DOY, UTseconds);
                    % input "altitude" as a vector, but lat, long, year,
                    % DOY, and UTseconds as single point value
    rho = rho';     % for vector sum calculation of total number density
    % If A is a matrix, sum(A) treats the columns of A as vectors, 
    % returning a row vector of the sums of each column.
    Temperature = T(:,2)';
    NumDensity = sum(rho(1:5,:)) + sum(rho(7:9,:));   % total number density (unit: m^-3)
    Pressure = rho(6,:).*287.*T(:,2)';            % pressure (unit: Pascal = J/m^3) 
                    % R*=287 J/K/kg, total mass density rho(6) has unit of kg/m^3

%     figure; plot(Temperature, altitude/1000)
%     axis([160 280 0 120])
%     xlabel('Temperature (K)'); ylabel('Altitude (km)')
%     figure; semilogx(NumDensity, altitude/1000)
%     xlabel('Number Density'); ylabel('Altitude (km)')
%     figure; semilogx(Pressure/100, altitude/1000)
%     xlabel('Pressure (mbar)'); ylabel('Altitude (km)')
end
%% Calculate MSISE-00 using observed solar flux and geomagnetic index
if (fluxindex == 2) % Use observed solar flux F107 and geomagnetic index AP
    [f107a, f107, aph] = readindex(year, month, day, UThour);
%     flags = ones(1,23);
%     flags(9) = -1;
    [T rho] = atmosnrlmsise00(altitude, latitude, longitude, year, DOY, UTseconds, f107a, f107, aph');
                    % input "altitude" as a vector, but lat, long, year,
                    % DOY, and UTseconds as single point value
    rho = rho';     % for vector sum calculation of total number density
    % If A is a matrix, sum(A) treats the columns of A as vectors, 
    % returning a row vector of the sums of each column.
    Temperature = T(:,2)';
    NumDensity = sum(rho(1:5,:)) + sum(rho(7:9,:));   % total number density (unit: m^-3)
    Pressure = rho(6,:).*287.*T(:,2)';            % pressure (unit: Pascal = J/m^3) 
                    % R*=287 J/K/kg, total mass density rho(6) has unit of kg/m^3
                    
%     figure; plot(Temperature, altitude/1000)
%     axis([160 280 0 120])
%     xlabel('Temperature (K)'); ylabel('Altitude (km)')
%     figure; semilogx(NumDensity, altitude/1000)
%     xlabel('Number Density'); ylabel('Altitude (km)')
%     figure; semilogx(Pressure/100, altitude/1000)
%     xlabel('Pressure (mbar)'); ylabel('Altitude (km)')
end
%% Return altitude, temperature, pressure, and total number density
ATPND = [altitude; Temperature; Pressure; NumDensity];

