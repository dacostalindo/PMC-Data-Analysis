% February 13, 2011 @ Boulder
% Main code to process the raw data (photon counts) of Fe Boltzmann lidar 
% for wavelength scans at Boulder and in the future at McMurdo
% to determine the PAL laser lineshape and linewidth

% Utilize "global variables" to make the code universal and consistent with
% Na Doppler lidar code, and in the future Fe Doppler lidar code
% Save "step data" to hard drive to expedite future processing

% This new flowchart moves readrawdata and preprocess (saturation correction,
% chopper correction, and integrating range bins) from <setprocess.m> to this 
% <main.m> and adds temporal integration to <main.m> for versitle handling
% of different situations, while keeping <setprocess.m> and <profileprocess.m> 
% universal for all cases.

% The different situations include, e.g., 
% (1)     Good data profiles for every unit integration (e.g., 90 s), so processing 
% each data file without integration in range bin or time; 
% (2) Insufficient SNR for raw binwidth so have to integrate in range bins; 
% (3) Insufficient SNR for each saved profiles, so have to integrate in time and range bins 
% to form one set of useful data for derivation of T, W, and D.
function [] = overviewFunction(year,month,day) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Initialization of Fe Boltzmann Lidar Data Processing    %%%%%%%
%%%%%%% Claim global variables; Define constants and parameters %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Year Month Day (YMD) list for McMurdo lidar data
% YMDayList=[2010 12  16  17  18  NaN     % year, month, and multiple days
%     2010	12	24  25  26  NaN
%     2011    01  15  16  17  18
%     2011    01  23  NaN  NaN  NaN];
% year = YMDayList(2,1); month = YMDayList(2,2); day = YMDayList(2,3:end);
% year = 1999; month = 12; day = [02];

% % read in the date-list for Rothera lidar data
% load /Users/chu/Tiger/ScienceProjects/Rothera/DataProcess/RotheraPMC/RotheraPMCmdych.dat
% FePMCmdych=RotheraPMCmdych;         % [month, day, year, number-of-channels, 374-on/off, sys1-Fe/Rayleigh, 372-on/off, sys2-Fe/Rayleigh, bkgfit]

%% User-controllable data processing parameters
global PMTsaturation_correction_yn Chopper_correction_yn Background_fit_yn Smooth_PhotonCnt_yn RayleighNormChoice 
global Location BaseAlt latitude longitude fluxindex PlotRayFigure_yn
Location = 1;                       % 1 = McMurdo, 2 = Boulder, 3 = Rothera, 4 = South Pole
        % Location selection affects "Base Altitude", "Latitude", "Longitude" and "Local Time" in data processing
fluxindex = 1;                      % 1 = use default F107 and AP for MSISE00, 2 = use observed F107 and AP
PMTsaturation_correction_yn = 0;    % 0 = no correction, 1 = do correction
Chopper_correction_yn = 0;          % 0 = no correction, 1 = do correction
Background_fit_yn = 0;              % 0 = constant background, 1 = tilted background so fitting is needed
Smooth_PhotonCnt_yn = 1;            % 0 = not to smooth, 1 = to smooth photon counts
PlotRayFigure_yn = 0;
RayleighNormChoice = 2;             % 1 = Take one point via mean(log); 2 = Take sum of a range

OverviewRawData = 1;
ScreenRawData = 0;                  % 0 = no screening, 1 = do raw data screening
SaveStatData = 1;                   % 0 = not to save statistic data, 1 = save stat data to harddrive
SaveSmdData = 0;                    % 0 = not to save smoothed data, 1 = save smoothed data to harddrive
SaveStatFigure = 1;                 % 0 = not to save statistic figure, 1 = save stat figure to harddrive
SaveContFigure = 1;                 % 0 = not to save contour figure, 1 = save contour figure to harddrive
PlotStatFigure = 0;                 % 0 = not to plot stat figure, 1 = plot stat figure
PlotRawProfile = 0;                 % 0 = not to plot each raw profile, 1 = plot each raw profile
SmoothOrder = 1;                    % 1 = Smooth in time and then range, 2 = Smooth in range and then time
% SmoothWidthTime = 2;                % 1 = 0.25-h Full Width, 2 = 0.5-h FW, 3 = 1-h FW, 4 = 2-h FW, 5 = 3-h FW
% SmoothWidthAlt = 1;                 % 1 = 1-km Full Width, 2 = 2-km FW, 3 = 3-km FW, 4 = 4-km FW, 5 = 10-km FW
SmoothWidthTime = 3;                % 1 = 0.25-h Full Width, 2 = 0.5-h FW, 3 = 1-h FW, 4 = 2-h FW, 5 = 3-h FW
SmoothWidthAlt = 2;                 % 1 = 1-km Full Width, 2 = 2-km FW, 3 = 3-km FW, 4 = 4-km FW, 5 = 10-km FW
ResolutionTime = 2;                 % 1 = 3-min Resolution, 2 = 6-min res, 3 = 0.25-h res, 4 = 0.5-h res, 5 = 0.02-hr res
ResolutionAlt = 2;                  % 1 = 48-m Resolution, 2 = 0.1-km res, 3 = 0.25-km res, 4 = 0.5-km res, 5 = 1-km res
Process4temp = 0;                   % 0 = not for temperature (save statistic and high-resolution data),
                    % 1 = for temperature (not to save stat data but save low-res data to different names)
                    
Process4TempDen = 0;        % Formal processing routine #1 for temp and density
SaveProcessedData = 0;

Process4TempSmooth = 0;     % Formal processing routine #2 for temp and density
SavePreprocessedData = 0;
SaveSmoothedData = 0;

Process4DailyMean = 0;      % Formal processing routine #3 for daily temperature or long-integration temperature
SaveDailymeanData = 0;

%% Location Selection affects base altitude, latitude and longitude (local time vs. UT)
% Latitude and longitude will affect "atmosnrlmsise00" outputs of temperature and atmospheric density
if (Location == 1)          % McMurdo, Antarctica
    station = 'McMurdo';
    BaseAlt = 190;          % unit: m
    latitude = -77.85;      % 77.85 deg South
    longitude = 166.67;     % 166.67 deg East
elseif (Location == 2)      % Boulder, Colorado
    station = 'Boulder';
    BaseAlt = 1600;         % unit: meter
    latitude = 40.1278;     % 40 deg 07.665 min North
    longitude = -105.2437;  % 105 deg 14.620 min West
elseif (Location == 3)      % Rothera, Antarctica
    station = 'Rothera';
    BaseAlt = 30;           % unit: m
    latitude = -67.5;       % 67.5 deg South
    longitude = -68.0;      % 68.0 deg West
elseif (Location == 4)      % South Pole, Antarctica
    station = 'SouthPole';
    BaseAlt = 2840;         % unit: m
    latitude = -90;         % 90 deg South
    longitude = 0;          % All longitude lines cross at the South Pole
end

global Fe_start_alt Fe_end_alt vert_bin_res Rayleigh_sum_start_alt Rayleigh_sum_end_alt Rayleigh_norm_alt
global Rayleigh_fit_start_alt Rayleigh_fit_end_alt Rayleigh_fit_alt BG_start_alt BG_end_alt BG_start_alt2 BG_end_alt2
global temp_default wind_default density_default time_adjust max_temp min_temp max_wind minimum_signal Hamming_width
global RTUpperLim1 RTLowerLim1 RWUpperLim1 RWLowerLim1
% BG_start_alt = 130e3;             % unit: m for all altitudes
% BG_start_alt = 155e3;             % For 1-hr integration temperature profile
% BG_start_alt = 160e3;             % unit: m for all altitudes
BG_start_alt = 130e3;             % unit: m for all altitudes
BG_end_alt = 180e3;
BG_start_alt2 = 70e3; 
BG_end_alt2 = 75e3;
% Rayleigh_sum_start_alt = 45.0e3;
% Rayleigh_sum_end_alt = 55.0e3;
% Rayleigh_norm_alt = 50.0e3;
% Rayleigh_fit_start_alt = 45e3;
% Rayleigh_fit_end_alt = 55e3;
% Rayleigh_fit_alt = 50.0e3;
Rayleigh_sum_start_alt = 40.0e3;
Rayleigh_sum_end_alt = 50.0e3;
Rayleigh_norm_alt = 45.0e3;
Rayleigh_fit_start_alt = 40e3;
Rayleigh_fit_end_alt = 50e3;
Rayleigh_fit_alt = 45.0e3;
% Fe_start_alt = 75e3;
% Fe_end_alt = 115e3;
% Fe_start_alt = 70e3;
% Fe_end_alt = 155e3;
Fe_start_alt = 75e3;
Fe_end_alt = 115e3;
Hamming_width = 1920;       % Nominal FWHM = 1920 m
vert_bin_res = 1;
temp_default = 200;               % unit: K
wind_default = 0;                 % unit: m/s
density_default = 1e7;            % unit: ?
time_adjust = 0.;
max_temp = 500;
min_temp = 100;
max_wind = 200;
minimum_signal = 50;
RTUpperLim1 = 1.799294;         % limits of RT & RW when AOM shift = 630MHz
RTLowerLim1 = 0.059951;
RWUpperLim1 = 1.977745;
RWLowerLim1 = -1.063667;

% %% Bin Number Corresponding to Absolute Altitude for
% % RayStart, RayEnd, FeStart, FeEnd, BgStart, BgEnd and Base Altitude
% bino4RayStart = round(Rayleigh_fit_start_alt/Rawbinwid);
% bino4RayEnd = round(Rayleigh_fit_end_alt/Rawbinwid);
% bino4RayNorm = round(Rayleigh_norm_alt/Rawbinwid);
% bino4FeStart = round(Fe_start_alt/Rawbinwid);
% bino4FeEnd = round(Fe_end_alt/Rawbinwid);
% bino4BgStart = round(BG_start_alt/Rawbinwid);
% bino4BgEnd = round(BG_end_alt/Rawbinwid);
% bino4BaseAlt = round(BaseAlt/Rawbinwid);

% AMU = 1.6605402e-27;              % Atomic mass unit
% Radius = 6371e3;                  % radius of the Earth (unit: m)
% SigmaRay = 5.5504e-31;            % unit: m^2, including the factor of 4*pi
% MAX_Countrate = 55e6;             % maximum count rate

global smonth zulumonth
smonth = ['JA';'FB';'MR';'AR';'MY';'JN';'JL';'AG';'SP';'OT';'NV';'DC'];
zulumonth = [31;28;31;30;31;30;31;31;30;31;30;31];		% day numbers in each month through whole year
azimuth = 0;
elevation = 90;
zuluday = sum(zulumonth(1:month-1))+day;	% convert to zulu day (1-365)

%% Fundamental Constants and Atomic Fe Constants
global pi c h Me Qe E0 kB AMU NA mFeMean mFeIsotope AbundanceFe DeltaE_diff 
global Aki_Fe372 gk_Fe372 gi_Fe372 RB_Fe372 Lambda_Center_Fe372 fosc_Fe372 IsotopeShift_Fe372
global Aki_Fe374 gk_Fe374 gi_Fe374 RB_Fe374 Lambda_Center_Fe374 fosc_Fe374 
pi = 3.14159265;                    % Accurate value of pi
c = 2.99792458e8;                   % light speed (m/s)
h = 6.62606896e-34;                 % Planck constant (unit: Js)
Me = 9.10938215e-31;                % mass of electron (unit: kg)
Qe = 1.602176487e-19 ;              % charge of electron (unit: C)
E0 = 8.854187817e-12;               % (unit: F/m)
kB = 1.3806504e-23;                 % Boltzmann constant (Unit:  J/K)
AMU = 1.660538782e-27;              % Atomic mass constant (unit: kg)
NA = 6.02214179e23;                 % Avogadro constant (unit: mol^-1)

mFeMean = 55.845*1e-3/NA;           % mean mass of Fe atom (unit: kg)
mFeIsotope = [53.9396,55.9349,56.9354,57.9333]*1e-3/NA;
                                    % mass of Fe atom for each isotope
                                    % (unit: kg)
AbundanceFe = [0.0585,0.9175,0.0212,0.0028];
                                    % natural abundance for Fe 54, 56, 57, and 58 isotopes
DeltaE_diff = 415.932;              % energy level difference between ground state (a5D4) and (a5D3) (unit: cm^-1)
                                    
Aki_Fe372 = 0.163e8;                % Enstein A coefficient for Fe 372-nm line
gk_Fe372 = 2*5+1;                   % degeneracy factor of Fe first excited state for 372nm line
gi_Fe372 = 2*4+1;                   % degeneracy factor of Fe ground state for 372nm line
RB_Fe372 = 0.9959;                  % Branching ratio for Fe 372-nm absorption line (refer to spectroscopy HW#5)
Lambda_Center_Fe372 = 744.1992e-9 / 2;  % center wavelength of Fe 372 nm line in vacuum (m)
                                    % measured by Burleigh WA-1500 wavemeter at Doppler-broadened peak from sky return
fosc_Fe372 = E0*Me*c/(2*pi*Qe.^2)*gk_Fe372/gi_Fe372*Lambda_Center_Fe372.^2*Aki_Fe372;
                                    % oscillator strength of Fe 372nm line (dimensionless)
IsotopeShift_Fe372 = [-725e6,0,495e6,669e6];     
                                    % isotopic shift with respect to Fe56 (unit: Hz) 
              % isotopic shifts for Fe 54 and 57 are from hollow cathod discharge cell [Smeets et al., Appl. Phys. B, 2003]
              % isotopic shift for Fe 58 is from Kaletta thesis [1969]
              
Aki_Fe374 = 0.142e8;                % Enstein A coefficient for Fe 374-nm line
gk_Fe374 = 2*4+1;                   % degeneracy factor of Fe first excited state for 374nm line
gi_Fe374 = 2*3+1;                   % degeneracy factor of Fe ground state for 374nm line
RB_Fe374 = 0.9079;                  % Branching ratio for Fe 374-nm absorption line (refer to spectroscopy HW#5)
Lambda_Center_Fe374 = 747.6394e-9 / 2;  % center wavelength of Fe 374 nm line in vacuum (unit: m)
                                    % measured by Burleigh WA-1500 wavemeter at Doppler-broadened peak from sky return
fosc_Fe374 = E0*Me*c/(2*pi*Qe.^2)*gk_Fe374/gi_Fe374*Lambda_Center_Fe374.^2*Aki_Fe374;

%% Linear PAL laser information (central wavelength, peak freq offset, linewidth, energy portion/pedestal, ...)
global PAL374Eportion PAL372Eportion PAL374RMSwidth PAL372RMSwidth 
global PAL374Detune PAL372Detune PAL374WL PAL372WL dEoverKB
PAL374Eportion = [0.814,0.186];     % Energy portions of narrow peak and pedestal for 374-nm PAL (dimensionless)
PAL372Eportion = [0.896,0.104];     % Energy portions of narrow peak and pedestal for 372-nm PAL (dimensionless)
PAL374RMSwidth = [408.5,15e3]*1e6;  % RMS width of Gaussian linewidth for the peak and pedestal (unit: Hz)
PAL372RMSwidth = [572.1,15e3]*1e6;  % RMS width of Gaussian linewidth for the peak and pedestal (unit: Hz)
PAL374Detune = 51.38e6;             % PAL 374nm laser central frequency detuned related to 56Fe peak (unit: Hz)
        % The Fe peak freq is at -51.38 MHz relative to the laser wavelength 744.6394nm/2
        % If assuming the Fe peak is the 56Fe peak, the laser central freq detune is +51.38 MHz from the 56Fe peak
PAL372Detune = -23.8e6;             % PAL 372nm laser central frequency detuned related to 56Fe peak (unit: Hz)
PAL374WL = 747.6394e-9 / 2;         % unit: meter
PAL372WL = 744.1992e-9 / 2;         % unit: meter
dEoverKB = h*c*100*DeltaE_diff/kB;       % dEoverKB = 598.438 K;

%% Fe Boltzmann lidar receiver parameters
global Deadtime_PMT Deadtime_Disc MAX_Countrate
% MAX_Countrate = 55e6;           % maximum count rate
Deadtime_PMT = 20e-9;           % PMT dead time (s)
Deadtime_Disc = 10e-9;          % Discriminator dead time (s)

%% Overview the entire raw datasets to screen data, do statistics, and plot overview contours
if (OverviewRawData == 1)
    for system = 1:2
        clear timeraw nameofprofile processedCnt smdCnt smdsmdCnt time0 Bkg Rayleigh FeCntpershot NormFeCnt
        for dd = 1:length(day)
            %         if (dd == 2)    % special case for December 25, 2010 at McMurdo
            %             fp = fopen([sprintf('/Users/haoyuli/Desktop/ScienceProjects/SkyObservations/%04d%s%02dA/RX%d/filelist',year,smonth(month,:),day(dd),system)],'r');
            %         elseif (dd == 3)
            %             fp = fopen([sprintf('/Users/haoyuli/Desktop/ScienceProjects/SkyObservations/%04d%s%02dB/RX%d/filelist',year,smonth(month,:),day(dd),system)],'r');
            %         else
            %             fp = fopen([sprintf('/Users/haoyuli/Desktop/ScienceProjects/SkyObservations/%04d%s%02d/RX%d/filelist',year,smonth(month,:),day(dd),system)],'r');
            %         end
            if (ScreenRawData == 0)
                fp = fopen(sprintf('/Users/dacostalindo/Desktop/Research/SkyObservations/%04d%s%02d/RX%d/Goodfilelist.txt',year,smonth(month,:),day(dd),system),'r');
                Numofprofile(dd) = fscanf(fp,'%d',[1]);
            elseif (ScreenRawData == 1)
                fp = fopen(sprintf('/Users/dacostalindo/Desktop/Research/SkyObservations/%04d%s%02d/RX%d/filelist.txt',year,smonth(month,:),day(dd),system),'r');
                Numofprofile(dd) = fscanf(fp,'%d',[1]);
            end
            
            
            for ii=1:Numofprofile(dd)
                if (dd > 1)
                    jj = ii + sum(Numofprofile(1:dd-1));        % for handling multiple days of data
                elseif (dd == 1)
                    jj = ii;
                end
               
                nameofprofile(jj,:) = fscanf(fp,'%s',[1]);      % read profile name listed in the "filelist" or "Goodfilelist"
                 if length(nameofprofile)== 96
                    disp('here')
                end
            end
            fclose(fp);
            
            if (PlotRawProfile == 1)
                scrsz = get(0,'ScreenSize');
                figure('Position',[200 70 scrsz(3)*0.8 scrsz(4)*0.4])
            end
            
            for mm = 1:Numofprofile(dd)
                if (dd > 1)
                    nn = mm + sum(Numofprofile(1:dd-1));
                elseif (dd == 1)
                    nn = mm;
                end
                fpin = fopen([sprintf('/Users/dacostalindo/Desktop/Research/SkyObservations/%04d%s%02d/RX%d/%s',year,smonth(month,:),day(dd),system,nameofprofile(nn,:))],'r');
                fgetl(fpin); fgetl(fpin); fgetl(fpin); fgetl(fpin); fgetl(fpin);
                timeraw(nn) = fscanf(fpin,'%f', [1]) + (day(dd)-day(1))*24;
                fgetl(fpin); fgetl(fpin); fgetl(fpin); fgetl(fpin);
                bintime = fscanf(fpin,'%f', [1]);fgetl(fpin);
                binrange = c*bintime*1e-9/2;
                BinWidth = round(binrange);
                fgetl(fpin);
                ShotsNum = fscanf(fpin,'%d', [1]);fgetl(fpin);
                fgetl(fpin);
                BinNum = fscanf(fpin,'%d', [1]); fgetl(fpin);
                fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);
                fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);
                % rawphotoncnt = fscanf(fpin,'%f',[1 inf]);
                rawphotoncnt = fscanf(fpin,'%f',[1 BinNum]);    % BinNum = 4096
                fclose(fpin);
                
                DataInfo = [year month day(dd) timeraw(nn) azimuth elevation BaseAlt BinNum ShotsNum BinWidth binrange];
                [pureFe,NormdCnt,TotalCnt,FeCnt,BkgCnt] = profileprocessOverview(rawphotoncnt,DataInfo,bintime);
                %processedCnt(nn,:) = pureFe;
                 processedCnt(nn,:) = NormdCnt;
                %% raw data statistics for screening or viewing %%%
%                 AltRaw = (1:BinNum)*binrange/1000*sin(elevation*pi/180)+
%                 BaseAlt/1000; % unit: km   
                AltRaw = (1:BinNum)*binrange/1000*sin(elevation*pi/180);        % unit: km % h = R*cos(theta)
                BinRangeNum = round(([Rayleigh_sum_start_alt Rayleigh_sum_end_alt Fe_start_alt Fe_end_alt BG_start_alt BG_end_alt] - BaseAlt)./binrange);
                        % for May 2nd and 28th, 2011 (thermospheric Fe layer and fast gravity waves)

                Bkg(nn) = mean(rawphotoncnt(BinRangeNum(5):BinRangeNum(6)));
                Rayleigh(nn) = sum(rawphotoncnt(BinRangeNum(1):BinRangeNum(2))-Bkg(nn));
                FeCntSum = sum(rawphotoncnt(BinRangeNum(3):BinRangeNum(4))-Bkg(nn));
                FeCntpershot(nn) = FeCntSum/ShotsNum;
                NormFeCnt(nn) = FeCntSum/Rayleigh(nn);
                if (PlotRawProfile == 1)
%                     plot(rawphotoncnt/ShotsNum/(bintime/1000), AltRaw); title(sprintf('%s',nameofprofile(mm,:)))
%                     axis([0 50 0 120]); grid on; pause(0.1)
%                     plot(NormdCnt, AltRaw,[0 0],[0 120],'-y'); title(sprintf('%s',nameofprofile(nn,:)))
%                     axis([-5 2 0 120])
%                     pause(0.1)
                    subplot(1,3,1)
%                     plot(NormdCnt,AltRaw); title(sprintf('%s',nameofprofile(nn,:)))
                    plot(rawphotoncnt,AltRaw+BaseAlt/1000,[0 0],[30 195],'-y'); title(sprintf('%s',nameofprofile(nn,:)))
                    xlabel('Raw Photon Counts','fontsize',14)
                    if (system == 1)
                        axis([-0.5 20 0 180])
                        % axis([-200 2000 0 195])
                    elseif (system == 2)
                        axis([-0.5 60 0 180])
                        % axis([-0.5 2000 0 195])
                    end
                    subplot(1,3,2)
                    plot(NormdCnt,AltRaw,[0 0],[30 195],'-y'); title(sprintf('%s',nameofprofile(nn,:)))
                    xlabel('Normalized Counts','fontsize',14)
                    if (system == 1)
                        % axis([-0.5 60 0 180])
                        axis([-20 30 0 195])
                    elseif (system == 2)
                        % axis([-0.5 60 0 180])
                        axis([-20 30 0 195])
                    end
                    subplot(1,3,3)
                    plot(pureFe, AltRaw,[0 0],[30 180],'-y'); title(sprintf('%s',nameofprofile(nn,:)))
                    xlabel('Pure Fe Counts','fontsize',14)
                    if (system == 1)
                        % axis([-0.5 5 30 180])
                        axis([-8 8 30 195])
                    elseif (system == 2)
                        % axis([-0.5 15 30 180])
                        axis([-8 8 30 195])
                    end
                    pause(0.1)
                end
            end
        end
        if (PlotStatFigure == 1)
            scrsz = get(0,'ScreenSize');
            figure('Position',[200 70 scrsz(3)*0.3 scrsz(4)*0.9])
            orient portrait; set(gcf,'PaperPositionMode','auto')
            subplot(4,1,1)
            plot(timeraw,Rayleigh); ylabel('Rayleigh Sum')
            if (length(day) == 1)
                if (system == 2)
                    title(sprintf('372-nm Channel Statistics on %d %s %04d (UT) ',day,smonth(month,:),year),'fontsize',18)
                elseif (system == 1)
                    title(sprintf('374-nm Channel Statistics on %d %s %04d (UT) ',day,smonth(month,:),year),'fontsize',18)
                end
            elseif (length(day) > 1)
                if (system == 2)
                    title(sprintf('372-nm Channel Statistics on %d-%d %s %04d (UT) ',day(1), day(end),smonth(month,:),year),'fontsize',18)
                elseif (system == 1)
                    title(sprintf('374-nm Channel Statistics on %d-%d %s %04d (UT) ',day(1), day(end),smonth(month,:),year),'fontsize',18)
                end
            end
            subplot(4,1,2)
            plot(timeraw,Bkg); ylabel('Background Mean')
            subplot(4,1,3)
            plot(timeraw,FeCntpershot); ylabel('Fe Count/shot');
            subplot(4,1,4)
            plot(timeraw,NormFeCnt); ylabel('Normalized Fe Count'); xlabel('Time (UT)')
            if (SaveStatFigure == 1)
                if (ScreenRawData == 1)
                    if (system == 1)
                        saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/FIG/%s%04d%s%02dstat374nmNew.fig',station,station,year,smonth(month,:),day(1)));
                        print('-depsc2',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/EPS/%s%04d%s%02dstat374nmNew',station,station,year,smonth(month,:),day(1)));
                        print('-djpeg',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/JPG/%s%04d%s%02dstat374nmNew',station,station,year,smonth(month,:),day(1)));
                    elseif (system == 2)
                        saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/FIG/%s%04d%s%02dstat372nmNew.fig',station,station,year,smonth(month,:),day(1)));
                        print('-depsc2',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/EPS/%s%04d%s%02dstat372nmNew',station,station,year,smonth(month,:),day(1)));
                        print('-djpeg',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/JPG/%s%04d%s%02dstat372nmNew',station,station,year,smonth(month,:),day(1)));
                    end
                elseif (ScreenRawData == 0)
                    if (system == 1)
                        saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/FIG/%s%04d%s%02dstat374nmGoodNew.fig',station,station,year,smonth(month,:),day(1)));
                        print('-depsc2',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/EPS/%s%04d%s%02dstat374nmGoodNew',station,station,year,smonth(month,:),day(1)));
                        print('-djpeg',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/JPG/%s%04d%s%02dstat374nmGoodNew',station,station,year,smonth(month,:),day(1)));
                    elseif (system == 2)
                        saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/FIG/%s%04d%s%02dstat372nmGoodNew.fig',station,station,year,smonth(month,:),day(1)));
                        print('-depsc2',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/EPS/%s%04d%s%02dstat372nmGoodNew',station,station,year,smonth(month,:),day(1)));
                        print('-djpeg',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/JPG/%s%04d%s%02dstat372nmGoodNew',station,station,year,smonth(month,:),day(1)));
                    end
                end
            end
            pause(0.1)
        end
        
        %%% Criterion #1 is one third of maximum(Rayleigh) %%%
        if (ScreenRawData == 1)
%             RayleighThreshold = max(Rayleigh)*0.5;  % for daytime configuration data
%             RayleighThreshold = max(Rayleigh)/3;  % for daytime configuration data
%             RayleighThreshold = max(Rayleigh)/4;  % for daytime configuration data
            RayleighThreshold = max(Rayleigh)/5;  % for daytime configuration data
%             RayleighThreshold = max(Rayleigh)/10;   % for nightime and daytime switching data
                        % daytime signal levels are much lower than nighttime configuration, 
                        % so 1/3 is too high standard for daytime data when compared to
                        % nighttime data
            %         goodfilenames = nameofprofile(find(Rayleigh > RayleighThreshold),:);
            %         badfilenames = nameofprofile(find(Rayleigh <= RayleighThreshold),:);
            indexgood = find(Rayleigh > RayleighThreshold);
            indexbad = find(Rayleigh <= RayleighThreshold);
            for dd = 1:length(day)
                clear goodfilenames badfilenames
                if (dd == 1)
                    goodfilenames = nameofprofile(indexgood(find(indexgood<=Numofprofile(dd))),:);
                    badfilenames = nameofprofile(indexbad(find(indexbad<=Numofprofile(dd))),:);
                elseif (dd > 1)
                    goodfilenames = nameofprofile(indexgood(find(indexgood>sum(Numofprofile(1:dd-1)) & indexgood<=sum(Numofprofile(1:dd)))),:);
                    badfilenames = nameofprofile(indexbad(find(indexbad>sum(Numofprofile(1:dd-1)) & indexbad<=sum(Numofprofile(1:dd)))),:);
                end
                if (SaveStatData == 1)
                    nametobesaved1 = sprintf('/Users/dacostalindo/Desktop/Research/SkyObservations/%04d%s%02d/RX%d/GoodfilelistNew',year,smonth(month,:),day(dd),system);
                    dlmwrite(nametobesaved1,goodfilenames,'')
                    nametobesaved2 = sprintf('/Users/dacostalindo/Desktop/Research/SkyObservations/%04d%s%02d/RX%d/BadfilelistNew',year,smonth(month,:),day(dd),system);
                    dlmwrite(nametobesaved2,badfilenames,'')
                end
            end
            
            scrsz = get(0,'ScreenSize');
            figure('Position',[200 70 scrsz(3)*0.3 scrsz(4)*0.9])
            orient portrait; set(gcf,'PaperPositionMode','auto')
            subplot(4,1,1)
            plot(timeraw(indexgood),Rayleigh(indexgood)); ylabel('Rayleigh Sum')
            if (length(day) == 1)
                if (system == 2)
                    title(sprintf('372-nm Channel Statistics on %d %s %04d (UT) ',day,smonth(month,:),year),'fontsize',18)
                elseif (system == 1)
                    title(sprintf('374-nm Channel Statistics on %d %s %04d (UT) ',day,smonth(month,:),year),'fontsize',18)
                end
            elseif (length(day) > 1)
                if (system == 2)
                    title(sprintf('372-nm Channel Statistics on %d-%d %s %04d (UT) ',day(1), day(end),smonth(month,:),year),'fontsize',18)
                elseif (system == 1)
                    title(sprintf('374-nm Channel Statistics on %d-%d %s %04d (UT) ',day(1), day(end),smonth(month,:),year),'fontsize',18)
                end
            end
            subplot(4,1,2)
            plot(timeraw(indexgood),Bkg(indexgood)); ylabel('Background Mean')
            subplot(4,1,3)
            plot(timeraw(indexgood),FeCntpershot(indexgood)); ylabel('Fe Count/shot');
            subplot(4,1,4)
            plot(timeraw(indexgood),NormFeCnt(indexgood)); ylabel('Normalized Fe Count'); xlabel('Time (UT)')
            if (SaveStatFigure == 1)
                if (system == 1)
                    saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/FIG/%s%04d%s%02dstat374nmGoodNew.fig',station,station,year,smonth(month,:),day(1)));
                    print('-depsc2',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/EPS/%s%04d%s%02dstat374nmGoodNew',station,station,year,smonth(month,:),day(1)));
                    print('-djpeg',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/JPG/%s%04d%s%02dstat374nmGoodNew',station,station,year,smonth(month,:),day(1)));
                elseif (system == 2)
                    saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/FIG/%s%04d%s%02dstat372nmGoodNew.fig',station,station,year,smonth(month,:),day(1)));
                    print('-depsc2',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/EPS/%s%04d%s%02dstat372nmGoodNew',station,station,year,smonth(month,:),day(1)));
                    print('-djpeg',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/StatisticsPlots/JPG/%s%04d%s%02dstat372nmGoodNew',station,station,year,smonth(month,:),day(1)));
                end
            end
            pause(0.1)
        end
        
        if sum(Numofprofile(1:end)) ~= 0
        % Plot processed counts at raw data resolutions
        Tmin = floor(min(timeraw));
        Tmax = ceil(max(timeraw));
        figure
        colormap('jet')
        pcolor(timeraw, AltRaw,processedCnt')
        shading interp; colorbar
        axis([Tmin Tmax 70 120])
        set(gca,'xtick',[Tmin:2:Tmax],'fontsize',18)
        set(gca,'tickdir','out')
        if (system == 2)
            caxis([-1 10])
        elseif (system == 1)
            caxis([-0.05 0.5])
        end
        xlabel('Time UT (hr)','fontsize',18)
        ylabel('Altitude (km)','fontsize',18)
        if (system == 2)
            title(sprintf('372-nm Fe Layers on %02d %s %04d (UT) ',day(1),smonth(month,:),year),'fontsize',18)
        elseif (system == 1)
            title(sprintf('374-nm Fe Layers on %02d %s %04d (UT) ',day(1),smonth(month,:),year),'fontsize',18)
        end
        pause(0.1)
        
        if (SmoothOrder == 1)
            % Hamming smooth the processed photon counts in time domain
            if (SmoothWidthTime == 1)
                FWtime = 0.25;  % unit: hour (Full Width in time is 0.5 hour)
            elseif (SmoothWidthTime == 2)
                FWtime = 0.5;  % unit: hour (Full Width in time is 0.5 hour)
            elseif (SmoothWidthTime == 3)
                FWtime = 1;  % unit: hour (Full Width in time is 0.5 hour)
            elseif (SmoothWidthTime == 4)
                FWtime = 2;  % unit: hour (Full Width in time is 0.5 hour)
            elseif (SmoothWidthTime == 5)
                FWtime = 3;  % unit: hour (Full Width in time is 0.5 hour)
            end
            Tmin = floor(min(timeraw));
            Tmax = ceil(max(timeraw));
            if (ResolutionTime == 1)
                time0 = [Tmin:0.05:Tmax];  % unit: hour (Resolution in time is 0.25 hour)
            elseif (ResolutionTime == 2)
                time0 = [Tmin:0.1:Tmax];    % unit: hour (Resolution in time is 0.25 hour)
            elseif (ResolutionTime == 3)
                time0 = [Tmin:0.25:Tmax];    % unit: hour (Resolution in time is 0.25 hour)
            elseif (ResolutionTime == 4)
                time0 = [Tmin:0.5:Tmax];    % unit: hour (Resolution in time is 0.25 hour)
            elseif (ResolutionTime == 5)
                time0 = [Tmin:0.02:Tmax];    % unit: hour (Resolution in time is 0.25 hour)
            end
            for kk = 1:length(AltRaw)
                [smdCnt(:,kk)] = HammingSmth_time(FWtime,time0,timeraw,processedCnt(:,kk));
            end
            
            figure
            pcolor(time0, AltRaw,smdCnt')
            shading interp; colorbar
            axis([Tmin Tmax 70 120])
            set(gca,'xtick',[Tmin:2:Tmax],'fontsize',18)
            set(gca,'tickdir','out')
            if (system == 2)
                caxis([-1 10])
            elseif (system == 1)
                caxis([-0.05 0.5])
            end
            xlabel('Time UT (hr)','fontsize',18)
            ylabel('Altitude (km)','fontsize',18)
            if (system == 2)
                title(sprintf('372-nm Fe Layers on %02d %s %04d (UT) ',day(1),smonth(month,:),year),'fontsize',18)
            elseif (system == 1)
                title(sprintf('374-nm Fe Layers on %02d %s %04d (UT) ',day(1),smonth(month,:),year),'fontsize',18)
            end
            pause(0.1)
            
            % Hamming smooth the photon counts in altitude range
            if (SmoothWidthAlt == 1)
                FWalt = 1;  % unit: km (Full Width in altitude is 1 km)
            elseif (SmoothWidthAlt == 2)
                FWalt = 2;  % unit: km (Full Width in altitude is 1 km)
            elseif (SmoothWidthAlt == 3)
                FWalt = 3;  % unit: km (Full Width in altitude is 1 km)
            elseif (SmoothWidthAlt == 4)
                FWalt = 4;  % unit: km (Full Width in altitude is 1 km)
            elseif (SmoothWidthAlt == 5)
                FWalt = 10; % unit: km (Full Width in altitude is 1 km)
            end
            Altmin = floor(min(AltRaw));
            Altmax = ceil(max(AltRaw));
            if (ResolutionAlt == 1)
                alt0 = AltRaw;  % unit: km
            elseif (ResolutionAlt == 2)
                alt0 = [Altmin:0.1:Altmax];    % unit: km
            elseif (ResolutionAlt == 3)
                alt0 = [Altmin:0.25:Altmax];    % unit: km
            elseif (ResolutionAlt == 4)
                alt0 = [Altmin:0.5:Altmax];    % unit: km
            elseif (ResolutionAlt == 5)
                alt0 = [Altmin:1:Altmax];    % unit: km
            end
            for tt = 1:length(time0)
                [smdsmdCnt(tt,:)] = HammingSmth(FWalt,alt0,AltRaw,smdCnt(tt,:));
            end
            figure
            pcolor(time0, alt0,smdsmdCnt')      % correct the original AltRaw typo
            shading interp; colorbar
            axis([Tmin Tmax 70 120])
            set(gca,'xtick',[Tmin:2:Tmax],'fontsize',18)
            set(gca,'tickdir','out')
            if (system == 2)
                caxis([-1 10])
            elseif (system == 1)
                caxis([-0.05 0.5])
            end
            xlabel('Time UT (hr)','fontsize',18)
            ylabel('Altitude (km)','fontsize',18)
            if (system == 2)
                title(sprintf('372-nm Fe Layers on %02d %s %04d (UT) ',day(1),smonth(month,:),year),'fontsize',18)
            elseif (system == 1)
                title(sprintf('374-nm Fe Layers on %02d %s %04d (UT) ',day(1),smonth(month,:),year),'fontsize',18)
            end
            pause(0.1)
        elseif (SmoothOrder == 2)
            % Hamming smooth the photon counts in altitude range (first)
            if (SmoothWidthAlt == 1)
                FWalt = 1;  % unit: km (Full Width in altitude is 1 km)
            elseif (SmoothWidthAlt == 2)
                FWalt = 2;  % unit: km (Full Width in altitude is 1 km)
            elseif (SmoothWidthAlt == 3)
                FWalt = 3;  % unit: km (Full Width in altitude is 1 km)
            elseif (SmoothWidthAlt == 4)
                FWalt = 4;  % unit: km (Full Width in altitude is 1 km)
            elseif (SmoothWidthAlt == 5)
                FWalt = 10;  % unit: km (Full Width in altitude is 1 km)
            end
            Altmin = floor(min(AltRaw));
            Altmax = ceil(max(AltRaw));
            if (ResolutionAlt == 1)
                alt0 = AltRaw;                  % unit: km
            elseif (ResolutionAlt == 2)
                alt0 = [Altmin:0.1:Altmax];     % unit: km
            elseif (ResolutionAlt == 3)
                alt0 = [Altmin:0.25:Altmax];    % unit: km
            elseif (ResolutionAlt == 4)
                alt0 = [Altmin:0.5:Altmax];     % unit: km
            elseif (ResolutionAlt == 5)
                alt0 = [Altmin:1:Altmax];       % unit: km
            end
            for tt = 1:length(timeraw)
                [smdCnt(tt,:)] = HammingSmth(FWalt,alt0,AltRaw,processedCnt(tt,:));
            end
            figure
            pcolor(timeraw, alt0,smdCnt')
            shading interp; colorbar
            axis([Tmin Tmax 70 120])
            set(gca,'xtick',[Tmin:2:Tmax],'fontsize',18)
            set(gca,'tickdir','out')
            if (system == 2)
                caxis([-1 10])
            elseif (system == 1)
                caxis([-0.05 0.5])
            end
            xlabel('Time UT (hr)','fontsize',18)
            ylabel('Altitude (km)','fontsize',18)
            if (system == 2)
                title(sprintf('372-nm Fe Layers on %02d %s %04d (UT) ',day(1),smonth(month,:),year),'fontsize',18)
            elseif (system == 1)
                title(sprintf('374-nm Fe Layers on %02d %s %04d (UT) ',day(1),smonth(month,:),year),'fontsize',18)
            end
            pause(0.1)
            
            % Hamming smooth the processed photon counts in time domain (second)
            if (SmoothWidthTime == 1)
                FWtime = 0.25;  % unit: hour (Full Width in time is 0.5 hour)
            elseif (SmoothWidthTime == 2)
                FWtime = 0.5;  % unit: hour (Full Width in time is 0.5 hour)
            elseif (SmoothWidthTime == 3)
                FWtime = 1;  % unit: hour (Full Width in time is 0.5 hour)
            elseif (SmoothWidthTime == 4)
                FWtime = 2;  % unit: hour (Full Width in time is 0.5 hour)
            elseif (SmoothWidthTime == 5)
                FWtime = 3;  % unit: hour (Full Width in time is 0.5 hour)
            end
            Tmin = floor(min(timeraw));
            Tmax = ceil(max(timeraw));
            if (ResolutionTime == 1)
                time0 = [Tmin:0.05:Tmax];  % unit: hour (Resolution in time is 0.25 hour)
            elseif (ResolutionTime == 2)
                time0 = [Tmin:0.1:Tmax];    % unit: hour (Resolution in time is 0.25 hour)
            elseif (ResolutionTime == 3)
                time0 = [Tmin:0.25:Tmax];    % unit: hour (Resolution in time is 0.25 hour)
            elseif (ResolutionTime == 4)
                time0 = [Tmin:0.5:Tmax];    % unit: hour (Resolution in time is 0.25 hour)
            elseif (ResolutionTime == 5)
                time0 = [Tmin:0.02:Tmax];    % unit: hour (Resolution in time is 0.25 hour)
            end
            for kk = 1:length(alt0)
                [smdsmdCnt(:,kk)] = HammingSmth_time(FWtime,time0,timeraw,smdCnt(:,kk));
            end
            figure
            pcolor(time0, alt0,smdsmdCnt')
            shading interp; colorbar
            axis([Tmin Tmax 70 120])
            set(gca,'xtick',[Tmin:2:Tmax],'fontsize',18)
            set(gca,'tickdir','out')
            if (system == 2)
                caxis([-1 10])
            elseif (system == 1)
                caxis([-0.05 0.5])
            end
            xlabel('Time UT (hr)','fontsize',18)
            ylabel('Altitude (km)','fontsize',18)
            if (system == 2)
                title(sprintf('372-nm Fe Layers on %02d %s %04d (UT) ',day(1),smonth(month,:),year),'fontsize',18)
            elseif (system == 1)
                title(sprintf('374-nm Fe Layers on %02d %s %04d (UT) ',day(1),smonth(month,:),year),'fontsize',18)
            end
            pause(0.1)
        end
        
        %%% Plot 2-D smoothed data and save data and figure %%%
        scrsz = get(0,'ScreenSize');
        figure('Position',[1 70 scrsz(3)*0.5 scrsz(4)*0.3])
        orient portrait
        set(gcf,'PaperPositionMode','auto')
        pcolor(time0, alt0,smdsmdCnt')
        colormap('jet')
        shading interp; colorbar
        axis([Tmin Tmax 75 115])
        set(gca,'xtick',[Tmin:2:Tmax],'fontsize',16)
        set(gca,'tickdir','out')
        xlabel('Time UT (hr)','fontsize',18)
        ylabel('Altitude (km)','fontsize',18)
        if (system == 2)
            caxis([-1 10])
        elseif (system == 1)
            caxis([-0.05 0.5])
        end
        if (length(day) == 1)
            if (system == 2)
                title(sprintf('372-nm Fe Layers on %d %s %04d (UT) ',day,smonth(month,:),year),'fontsize',18)
            elseif (system == 1)
                title(sprintf('374-nm Fe Layers on %d %s %04d (UT) ',day,smonth(month,:),year),'fontsize',18)
            end
        elseif (length(day) > 1)
            if (system == 2)
                title(sprintf('372-nm Fe Layers on %d-%d %s %04d (UT) ',day(1), day(end),smonth(month,:),year),'fontsize',18)
            elseif (system == 1)
                title(sprintf('374-nm Fe Layers on %d-%d %s %04d (UT) ',day(1), day(end),smonth(month,:),year),'fontsize',18)
            end
        end
        
        if (SaveContFigure == 1)
            if (Process4temp == 0)
                if (ScreenRawData == 1)
                    if (system == 1)
                        saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/FIG/%sFe374nm%04d%s%02dAllNew.fig',station,station,year,smonth(month,:),day(1)));
                        print('-depsc2',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/EPS/%sFe374nm%04d%s%02dAllNew',station,station,year,smonth(month,:),day(1)));
                        print('-djpeg',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/JPG/%sFe374nm%04d%s%02dAllNew',station,station,year,smonth(month,:),day(1)));
                    elseif (system == 2)
                        saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/FIG/%sFe372nm%04d%s%02dAllNew.fig',station,station,year,smonth(month,:),day(1)));
                        print('-depsc2',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/EPS/%sFe372nm%04d%s%02dAllNew',station,station,year,smonth(month,:),day(1)));
                        print('-djpeg',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/JPG/%sFe372nm%04d%s%02dAllNew',station,station,year,smonth(month,:),day(1)));
                    end
                elseif (ScreenRawData == 0)
                    if (system == 1)
                        saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/FIG/%sFe374nm%04d%s%02dGoodNew.fig',station,station,year,smonth(month,:),day(1)));
                        print('-depsc2',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/EPS/%sFe374nm%04d%s%02dGoodNew',station,station,year,smonth(month,:),day(1)));
                        print('-djpeg',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/JPG/%sFe374nm%04d%s%02dGoodNew',station,station,year,smonth(month,:),day(1)));
                    elseif (system == 2)
                        saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/FIG/%sFe372nm%04d%s%02dGoodNew.fig',station,station,year,smonth(month,:),day(1)));
                        print('-depsc2',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/EPS/%sFe372nm%04d%s%02dGoodNew',station,station,year,smonth(month,:),day(1)));
                        print('-djpeg',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/JPG/%sFe372nm%04d%s%02dGoodNew',station,station,year,smonth(month,:),day(1)));
                    end
                end
            elseif (Process4temp == 1)
                if (system == 1)
                    saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/FIG/%sFe374nm%04d%s%02dTempNew.fig',station,station,year,smonth(month,:),day(1)));
                    %                 print('-djpeg',sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/OverviewPlots/JPG/%sFe374nm%04d%s%02dTemp',station,station,year,smonth(month,:),day(dd)));
                elseif (system == 2)
                    saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewPlots/FIG/%sFe372nm%04d%s%02dTempNew.fig',station,station,year,smonth(month,:),day(1)));
                    %                 print('-djpeg',sprintf('/Users/Chu/Tiger/ScienceProjects/%s/Results/McMurdo/OverviewStatistics/OverviewPlots/JPG/%sFe372nm%04d%s%02dTemp',station,station,year,smonth(month,:),day(dd)));
                end
            end
        end
        if (SaveSmdData == 1)
            clear datatobesaved
%             Fe_start_bin = round(Fe_start_alt/binrange);
%             Fe_end_bin = round(Fe_end_alt/binrange);
            [NUM,Fe_start_bin] = min(abs(alt0 - Fe_start_alt/1000));
            [NUM,Fe_end_bin] = min(abs(alt0 - Fe_end_alt/1000));
            datatobesaved = [355 time0;alt0(Fe_start_bin:Fe_end_bin)' smdsmdCnt(:,Fe_start_bin:Fe_end_bin)'];
            if (Process4temp == 0)
                if (ScreenRawData == 1)
                    if (system == 1)
                        filename=sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewData/%sFe374nm%04d%s%02dAllNew.txt',station,station,year,smonth(month,:),day(1));
                    elseif (system == 2)
                        filename=sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewData/%sFe372nm%04d%s%02dAllNew.txt',station,station,year,smonth(month,:),day(1));
                    end
                elseif (ScreenRawData == 0)
                    if (system == 1)
                        filename=sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewData/%sFe374nm%04d%s%02dGoodNew.txt',station,station,year,smonth(month,:),day(1));
                    elseif (system == 2)
                        filename=sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewData/%sFe372nm%04d%s%02dGoodNew.txt',station,station,year,smonth(month,:),day(1));
                    end
                end
            elseif (Process4temp == 1)
                if (system == 1)
                    filename=sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewData/%sFe374nm%04d%s%02dTempNew.txt',station,station,year,smonth(month,:),day(1));
                elseif (system == 2)
                    filename=sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/OverviewData/%sFe372nm%04d%s%02dTempNew.txt',station,station,year,smonth(month,:),day(1));
                end
            end
            save(filename,'datatobesaved','-ASCII')
        end
    
        
        if ScreenRawData == 0
            
            updatedFileList = continuosInterval( smdsmdCnt, time0,timeraw,nameofprofile);
            fpout = fopen(sprintf('/Users/dacostalindo/Desktop/Research/SkyObservations/%04d%s%02d/RX%d/GoodfilelistUpdated.txt',year,smonth(month,:),day(dd),system),'w');
            fprintf(fpout,'%d \n',size(updatedFileList,1));
            
             for i = 1:size(updatedFileList,1)
                    fprintf(fpout,[updatedFileList(i,:) '\n']);
             end
            
            fclose(fpout);
  
 
       end 
        else
           
        fpout = fopen(sprintf('/Users/dacostalindo/Desktop/Research/SkyObservations/%04d%s%02d/RX%d/GoodfilelistUpdated.txt',year,smonth(month,:),day(dd),system),'w');
        fprintf(fpout,'%d \n',0);     
            
            
            
        end
    %%%%end
   end 
    
 
    
    
    clear DataInfo processedCnt smdCnt smdsmdCnt nameofprofile AltRaw time0 timeraw
end

%% Process raw photon counts for Fe temp and data with standard Rayleigh normalization (at one particular altitude)
if (Process4TempDen == 1)
    SmoothOrder = 1;                    % 1 = Smooth in time and then range, 2 = Smooth in range and then time
    SmoothWidthTime = 4;                % 1 = 0.25-h Full Width, 2 = 0.5-h FW, 3 = 1-h FW, 4 = 2-h FW, 5 = 3-h FW
    SmoothWidthAlt = 3;                 % 1 = 1-km Full Width, 2 = 1.5km FW, 3 = 2km FW, 4 = 2.5km FW, 5 = 3km FW
    ResolutionTime = 3;                 % 1 = 1.5-min Resolution, 2 = 6-min res, 3 = 0.25-h res, 4 = 0.5-h res, 5 = 1-hr res
    ResolutionAlt = 1;                  % 1 = 48-m Resolution, 2 = 0.1-km res, 3 = 0.5-km res, 4 = 1-km res, 5 = 2-km res
    Process4temp = 1;                   % 0 = not for temperature (save statistic and high-resolution data),

    %% Read in raw data profiles for an entire dataset, form DataInfo and RawData matrixes
    % DataInfo matrix and RawData matrix are used in the PMC 1 or 2 channel
    % setprocess and unitprocess subroutines
    for system = 1:2
        for dd = 1:length(day)
            goodfilelist = textread(sprintf('/Users/haoyuli/Desktop/ScienceProjects/SkyObservations/%04d%s%02d/RX%d/Goodfilelist',year,smonth(month,:),day(dd),system),'%s');
            Numofprofile(dd) = size(goodfilelist,1); clear goodfilelist
            fp = fopen(sprintf('/Users/haoyuli/Desktop/ScienceProjects/SkyObservations/%04d%s%02d/RX%d/Goodfilelist',year,smonth(month,:),day(dd),system),'r');
            for ii=1:Numofprofile(dd)
                if (dd > 1)
                    jj = ii + sum(Numofprofile(1:dd-1));        % for handling multiple days of data
                elseif (dd == 1)
                    jj = ii;
                end
                nameofprofile(jj,:) = fscanf(fp,'%s',[1]);      % read profile name listed in the "filelist" or "Goodfilelist"
            end
            fclose(fp);

            for mm = 1:Numofprofile(dd)
                if (dd > 1)
                    nn = mm + sum(Numofprofile(1:dd-1));
                elseif (dd == 1)
                    nn = mm;
                end
                fpin = fopen([sprintf('/Users/haoyuli/Desktop/ScienceProjects/SkyObservations/%04d%s%02d/RX%d/%s',year,smonth(month,:),day(dd),system,nameofprofile(nn,:))],'r');
                fgetl(fpin); fgetl(fpin); fgetl(fpin); fgetl(fpin); fgetl(fpin);
                time(nn) = fscanf(fpin,'%f', [1]) + (day(dd)-day(1))*24; 
                fgetl(fpin); fgetl(fpin); fgetl(fpin); fgetl(fpin); 
                BinTime(nn) = fscanf(fpin,'%f', [1]);fgetl(fpin);
                binrange = c*BinTime(nn)*1e-9/2;
                BinWidth = round(binrange);
                fgetl(fpin);
                ShotsNum(nn) = fscanf(fpin,'%f', [1]);fgetl(fpin);
                fgetl(fpin);
                BinNum = fscanf(fpin,'%f', [1]); BinNum = round(BinNum);fgetl(fpin);
                fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);
                fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);
                RawPhotonCnt(nn,:) = fscanf(fpin,'%f',[1 BinNum]);
                fclose(fpin);
                DataInfo(nn,:) = [year month day(dd) time(nn) ShotsNum(nn) BinTime(nn)];
            end
        end
        if (system == 1)
            DataInfo1 = DataInfo;
            RawData1 = RawPhotonCnt;
        elseif (system == 2)
            DataInfo2 = DataInfo;
            RawData2 = RawPhotonCnt;
        end
        clear DataInfo RawPhotonCnt time ShotsNum BinTime nameofprofile Numofprofile
    end

    %% Load atmosphere number density, temperature, pressure vs altitude database
%     fmsis=fopen([sprintf('/Users/chu/Tiger/ScienceProjects/Rothera/DataProcess/MSIS00/MSISE00zTPND103%03d',zuluday(end))],'r'); 		
%     fscanf(fmsis,'%s %s	%s	%s',[4]);		% to skip
%     ATPND=fscanf(fmsis,'%f	%f	%f	%f',[4 inf]);	% Altitude (km), temperature (K), pressure (mbar) and number density (cm^-3)
%     ATPND=ATPND';
%     fclose(fmsis);
%     %% Compute MSISE00 data using MatLab built-in "atmosnrlmsise00.m" model
%     ATPND = msise00cal(year, month, day, time, latitude, longitude, fluxindex);
%     ATPND = ATPND';
% 
%     %% Find the temperature, pressure and number density at 50km
%     ATPND_R = ATPND(find(ATPND(:,1)==round(Rayleigh_norm_alt/1000)),1:4);
%     AR = ATPND_R(1);
%     TR = ATPND_R(2);		% Temperature at Rayleigh Normalization altitude
%     PR = ATPND_R(3);		% Pressure at Rayleigh Normalization altitude
%     NDR = ATPND_R(4);		% Number density at Rayleigh Normalization altitude
%     
%     %% Use interpolation to obtain the atmosphere number density at each bin altitude
%     ND = interp1(ATPND(:,1),ATPND(:,4),(1:BinNum)*binrange/1000);
%     TP = interp1(ATPND(:,1),ATPND(:,2),(1:BinNum)*binrange/1000);

    %% Divide data into unit time interval
    time1 = DataInfo1(:,4);
    time2 = DataInfo2(:,4);

    tmin = min(min(floor(time1)),min(floor(time2)));    % minimum UT hour
    tmax = max(max(ceil(time1)),max(ceil(time2)));      % maximum UT hour
    if (year == 2011 && month == 05 && day == 2)   % special for May 2, 2011 to make 13.5 to 14.5UT integration, centered on 14UT
        tmin = min(min(floor(time1)),min(floor(time2))) + 0.5;    % minimum UT hour
        tmax = max(max(ceil(time1)),max(ceil(time2)));      % maximum UT hour
    end
    if (year == 2011 && month == 05 && day == 28)   % special for May 28, 2011 to make 14.5 to 15.5UT integration, centered on 15UT
        tmin = min(min(floor(time1)),min(floor(time2))) + 0.5;    % minimum UT hour
        tmax = max(max(ceil(time1)),max(ceil(time2)));      % maximum UT hour
    end
    dt = 1;         % the time interval for each integration (normal PMC and Fe computation)
    t0 = [tmin:dt:tmax];      % set the time bin/interval for unit integration

%     HammingFWHM = 20;       % raw bin width = 48 m, so 48 x 20 = 960 m
%     HammingFWHM = 40;       % raw bin width = 48 m, so 48 x 20 = 1920 m
    HammingFWHM = 100;       % raw bin width = 48 m, so 48 x 100 = 4800 m
%     HammingFWHM = 150;       % raw bin width = 48 m, so 48 x 150 = 7200 m

    %% Process for entire datasets
    if (PlotRayFigure_yn == 1), figure; end
    for tt = 1:length(t0)-1
        clear sp1 sp2 RawCnt1 RawCnt2 Info1 Info2
        NumofProf1(tt) = length(find(time1>=t0(tt) & time1<t0(tt+1)));
        NumofProf2(tt) = length(find(time2>=t0(tt) & time2<t0(tt+1)));
        %% Compute MSISE00 data using MatLab built-in "atmosnrlmsise00.m" model
        timecurrent = t0(tt) + dt/2;
        [Y,I] = min(abs(timecurrent - time2)); daycurrent = DataInfo2(I,3);
        ATPND = msise00cal(year, month, daycurrent, timecurrent, latitude, longitude, fluxindex);
        ATPND = ATPND';
        
        %% Find the temperature, pressure and number density at 50km
        ATPND_R = ATPND(find(ATPND(:,1)==round(Rayleigh_norm_alt)),1:4);
        AR = ATPND_R(1);
        TR = ATPND_R(2);		% Temperature at Rayleigh Normalization altitude
        PR = ATPND_R(3);		% Pressure at Rayleigh Normalization altitude
        NDR = ATPND_R(4);		% Number density at Rayleigh Normalization altitude
        
        %% Use interpolation to obtain the atmosphere number density at each bin altitude
        ND = interp1(ATPND(:,1),ATPND(:,4),(1:BinNum)*binrange);
        TP = interp1(ATPND(:,1),ATPND(:,2),(1:BinNum)*binrange);

        Fe_start_bin = round(Fe_start_alt/binrange);
        Fe_end_bin = round(Fe_end_alt/binrange);

        if (NumofProf1(tt)>0 & NumofProf2(tt)>0)
            Info1 = DataInfo1(find(time1>=t0(tt) & time1<t0(tt+1)),:);
            RawCnt1 = RawData1(find(time1>=t0(tt) & time1<t0(tt+1)),:);
            Info2 = DataInfo2(find(time2>=t0(tt) & time2<t0(tt+1)),:);
            RawCnt2 = RawData2(find(time2>=t0(tt) & time2<t0(tt+1)),:);
            %% Standard pre-process of raw photon counts for both channels
            system = 1;
            [Normdsumsig374(tt,:),sumsig374(tt,:),RayNormSig374(tt,:),M(tt,:),sumsmdM(tt,:),sumdM(tt,:),B(tt,:),RayRemoval374(tt,:)] = PMCunitprocess(Info1,RawCnt1,HammingFWHM,system);
            if (PlotRayFigure_yn == 1), pause(0.5); end
            PuresigFe374(tt,:) = Normdsumsig374(tt,:) - ND/NDR;
            system = 2;
            [Normdsumsig372(tt,:),sumsig372(tt,:),RayNormSig372(tt,:),MM(tt,:),sumsmdMM(tt,:),sumdMM(tt,:),BB(tt,:),RayRemoval372(tt,:)] = PMCunitprocess(Info2,RawCnt2,HammingFWHM,system);
            if (PlotRayFigure_yn == 1), pause(0.5); end
            PuresigFe372(tt,:) = Normdsumsig372(tt,:) - ND/NDR;
        elseif (NumofProf1(tt)==0 & NumofProf2(tt)>0)
            Normdsumsig374(tt,1:BinNum) = NaN;
            PuresigFe374(tt,1:BinNum) = NaN;
            M(tt,1:BinNum) = NaN;
            sumsmdM(tt,1:BinNum) = NaN;
            sumdM(tt,1:BinNum) = NaN;
            B(tt) = NaN;
            RayRemoval374(tt,1:(Fe_end_bin-Fe_start_bin+1)) = NaN;
            Info2 = DataInfo2(find(time2>=t0(tt) & time2<t0(tt+1)),:);
            RawCnt2 = RawData2(find(time2>=t0(tt) & time2<t0(tt+1)),:);
            %% Pre-process of raw photon counts for 372-nm channels
            system = 2;
            [Normdsumsig372(tt,:),sumsig372(tt,:),RayNormSig372(tt,:),MM(tt,:),sumsmdMM(tt,:),sumdMM(tt,:),BB(tt,:),RayRemoval372(tt,:)] = PMCunitprocess(Info2,RawCnt2,HammingFWHM,system);
            if (PlotRayFigure_yn == 1), pause(0.5); end
            PuresigFe372(tt,:) = Normdsumsig372(tt,:) - ND/NDR;
        elseif (NumofProf1(tt)>0 & NumofProf2(tt)==0)
            Info1 = DataInfo1(find(time1>=t0(tt) & time1<t0(tt+1)),:);
            RawCnt1 = RawData1(find(time1>=t0(tt) & time1<t0(tt+1)),:);
            %% Pre-process of raw photon counts for 374-nm channels
            system = 1;
            [Normdsumsig374(tt,:),sumsig374(tt,:),RayNormSig374(tt,:),M(tt,:),sumsmdM(tt,:),sumdM(tt,:),B(tt,:),RayRemoval374(tt,:)] = PMCunitprocess(Info1,RawCnt1,HammingFWHM,system);
            if (PlotRayFigure_yn == 1), pause(0.5); end
            PuresigFe374(tt,:) = Normdsumsig374(tt,:) - ND/NDR;
            Normdsumsig372(tt,1:BinNum) = NaN;
            PuresigFe372(tt,1:BinNum) = NaN;
            MM(tt,1:BinNum) = NaN;
            sumsmdMM(tt,1:BinNum) = NaN;
            sumdMM(tt,1:BinNum) = NaN;
            BB(tt) = NaN;
            RayRemoval372(tt,1:(Fe_end_bin-Fe_start_bin+1)) = NaN;
        elseif (NumofProf1(tt)==0 & NumofProf2(tt)==0)
            Normdsumsig374(tt,1:BinNum) = NaN;
            PuresigFe374(tt,1:BinNum) = NaN;
            M(tt,1:BinNum) = NaN;
            sumsmdM(tt,1:BinNum) = NaN;
            sumdM(tt,1:BinNum) = NaN;
            B(tt) = NaN;
            RayRemoval374(tt,1:(Fe_end_bin-Fe_start_bin+1)) = NaN;
            Normdsumsig372(tt,1:BinNum) = NaN;
            PuresigFe372(tt,1:BinNum) = NaN;
            MM(tt,1:BinNum) = NaN;
            sumsmdMM(tt,1:BinNum) = NaN;
            sumdMM(tt,1:BinNum) = NaN;
            BB(tt) = NaN;
            RayRemoval372(tt,1:(Fe_end_bin-Fe_start_bin+1)) = NaN;
        end
    end
    timebins = t0(1:end-1) + dt/2;  % there was a typo before 5/18/2011 (-)
    altitude = (1:BinNum)*binrange/1000;
    %% Plot processed counts at raw data resolutions
    figure
    pcolor(timebins,(1:BinNum)*binrange/1000,PuresigFe374')
    shading interp; colorbar
    axis([tmin tmax 70 120])
    set(gca,'ytick',[70:5:120],'fontsize',18)
    set(gca,'tickdir','out')
    caxis([-0.05 1])
    xlabel('Time UT (hr)','fontsize',18)
    ylabel('Altitude (km)','fontsize',18)
    title(sprintf('374-nm Fe Layers on %02d %s %04d (UT) ',day(1),smonth(month,:),year),'fontsize',18)
    pause(0.1)

    figure
    pcolor(timebins,(1:BinNum)*binrange/1000,PuresigFe372')
    shading interp; colorbar
    axis([tmin tmax 70 120])
    set(gca,'ytick',[70:5:120],'fontsize',18)
    set(gca,'tickdir','out')
    caxis([-1 25])
    xlabel('Time UT (hr)','fontsize',18)
    ylabel('Altitude (km)','fontsize',18)
    title(sprintf('372-nm Fe Layers on %02d %s %04d (UT) ',day(1),smonth(month,:),year),'fontsize',18)
    pause(0.1)
    
    if (SaveProcessedData == 1)
        clear datatobesaved
        Fe_start_bin = round(Fe_start_alt/binrange);
        Fe_end_bin = round(Fe_end_alt/binrange);
        datatobesaved = [355 timebins;altitude(Fe_start_bin:Fe_end_bin)' PuresigFe374(:,Fe_start_bin:Fe_end_bin)'];
        filename=sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/ProcessedData/%sFeCnt374nm%04d%s%02dTDnew.txt',station,station,year,smonth(month,:),day(1));          
        save(filename,'datatobesaved','-ASCII')
        datatobesaved = [355 timebins;altitude(Fe_start_bin:Fe_end_bin)' PuresigFe372(:,Fe_start_bin:Fe_end_bin)'];
        filename=sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/ProcessedData/%sFeCnt372nm%04d%s%02dTDnew.txt',station,station,year,smonth(month,:),day(1));
        save(filename,'datatobesaved','-ASCII')
        datatobesaved = [355 timebins;altitude(Fe_start_bin:Fe_end_bin)' Normdsumsig374(:,Fe_start_bin:Fe_end_bin)'];
        filename=sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/ProcessedData/%sNormCnt374nm%04d%s%02dnew.txt',station,station,year,smonth(month,:),day(1));          
        save(filename,'datatobesaved','-ASCII')
        datatobesaved = [355 timebins;altitude(Fe_start_bin:Fe_end_bin)' Normdsumsig372(:,Fe_start_bin:Fe_end_bin)'];
        filename=sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/ProcessedData/%sNormCnt372nm%04d%s%02dnew.txt',station,station,year,smonth(month,:),day(1));
        save(filename,'datatobesaved','-ASCII')
        datatobesaved = [355 timebins;altitude' M'];
        filename=sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/ProcessedData/%sM374nm%04d%s%02dnew.txt',station,station,year,smonth(month,:),day(1));          
        save(filename,'datatobesaved','-ASCII')
        datatobesaved = [355 timebins;altitude' MM'];
        filename=sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/ProcessedData/%sM372nm%04d%s%02dnew.txt',station,station,year,smonth(month,:),day(1));
        save(filename,'datatobesaved','-ASCII')
        datatobesaved = [355 timebins;altitude' sumsmdM'];
        filename=sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/ProcessedData/%ssumsmdM374nm%04d%s%02dnew.txt',station,station,year,smonth(month,:),day(1));          
        save(filename,'datatobesaved','-ASCII')
        datatobesaved = [355 timebins;altitude' sumsmdMM'];
        filename=sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/ProcessedData/%ssumsmdM372nm%04d%s%02dnew.txt',station,station,year,smonth(month,:),day(1));
        save(filename,'datatobesaved','-ASCII')
        datatobesaved = [timebins; B'];
        filename=sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/ProcessedData/%sBkg374nm%04d%s%02dnew.txt',station,station,year,smonth(month,:),day(1));          
        save(filename,'datatobesaved','-ASCII')
        datatobesaved = [timebins; BB'];
        filename=sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/ProcessedData/%sBkg372nm%04d%s%02dnew.txt',station,station,year,smonth(month,:),day(1));
        save(filename,'datatobesaved','-ASCII')
        datatobesaved = [355 timebins;altitude(Fe_start_bin:Fe_end_bin)' RayRemoval374'];
        filename=sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/ProcessedData/%sRayRemoval374nm%04d%s%02dnew.txt',station,station,year,smonth(month,:),day(1));          
        save(filename,'datatobesaved','-ASCII')
        datatobesaved = [355 timebins;altitude(Fe_start_bin:Fe_end_bin)' RayRemoval372'];
        filename=sprintf('/Users/haoyuli/Desktop/PMC_processing/Results/%s/OverviewStatistics/ProcessedData/%sRayRemoval372nm%04d%s%02dnew.txt',station,station,year,smonth(month,:),day(1));
        save(filename,'datatobesaved','-ASCII')
    end

end

%% Process raw data using Hamming smooth method for temperature
if (Process4TempSmooth == 1)
	SmoothOrder = 1;                    % 1 = Smooth in time and then range, 2 = Smooth in range and then time
%     SmoothWidthTime = 3;                % 1 = 0.25-h Full Width, 2 = 0.5-h FW, 3 = 1-h FW, 4 = 2-h FW, 5 = 3-h FW
%     SmoothWidthAlt = 2;                 % 1 = 1-km Full Width, 2 = 2km FW, 3 = 4km FW, 4 = 6km FW, 5 = 10km FW
    SmoothWidthTime = 2;                % 1 = 0.25-h Full Width, 2 = 0.5-h FW, 3 = 1-h FW, 4 = 2-h FW, 5 = 3-h FW
    SmoothWidthAlt = 1;                 % 1 = 1-km Full Width, 2 = 2km FW, 3 = 4km FW, 4 = 6km FW, 5 = 10km FW
    ResolutionTime = 2;                 % 1 = 0.05-h Resolution, 2 = 0.1-h res, 3 = 0.2-h res, 4 = 0.5-h res, 5 = 1-hr res
    ResolutionAlt = 2;                  % 1 = 48-m Resolution, 2 = 0.1-km res, 3 = 0.25-km res, 4 = 0.5-km res, 5 = 1-km res
    % 2-h and 2-km FW were the original selection
    % Hamming smooth the processed photon counts in time domain
    if (SmoothWidthTime == 1)
        FWtime = 0.25;  % unit: hour (Full Width in time is 0.5 hour)
    elseif (SmoothWidthTime == 2)
        FWtime = 0.5;  % unit: hour (Full Width in time is 0.5 hour)
    elseif (SmoothWidthTime == 3)
        FWtime = 1;  % unit: hour (Full Width in time is 0.5 hour)
    elseif (SmoothWidthTime == 4)
        FWtime = 2;  % unit: hour (Full Width in time is 0.5 hour)
    elseif (SmoothWidthTime == 5)
        FWtime = 3;  % unit: hour (Full Width in time is 0.5 hour)
    end
    % Hamming smooth the photon counts in altitude range
    if (SmoothWidthAlt == 1)
        FWalt = 1;  % unit: km (Full Width in altitude is 1 km)
    elseif (SmoothWidthAlt == 2)
        FWalt = 2;  % unit: km (Full Width in altitude is 1 km)
    elseif (SmoothWidthAlt == 3)
        FWalt = 4;  % unit: km (Full Width in altitude is 1 km)
    elseif (SmoothWidthAlt == 4)
        FWalt = 6;  % unit: km (Full Width in altitude is 1 km)
    elseif (SmoothWidthAlt == 5)
        FWalt = 10;  % unit: km (Full Width in altitude is 1 km)
    end
    %% Read in raw data profiles for an entire dataset, form DataInfo and RawData matrixes
    % DataInfo matrix and RawData matrix are used in the PMC 1 or 2 channel
    % setprocess and unitprocess subroutines
    if (PlotRayFigure_yn == 1), figure; end
	HammingFWHM = 20;
    for system = 1:2
        for dd = 1:length(day)
            goodfilelist = textread(sprintf('/Users/dacostalindo/Desktop/Research/SkyObservations/%04d%s%02d/RX%d/GoodfilelistNew',year,smonth(month,:),day(dd),system),'%s');
            Numofprofile(dd) = size(goodfilelist,1); clear goodfilelist
            fp = fopen(sprintf('/Users/dacostalindo/Desktop/Research/SkyObservations/%04d%s%02d/RX%d/GoodfilelistNew',year,smonth(month,:),day(dd),system),'r');
            for ii=1:Numofprofile(dd)
                if (dd > 1)
                    jj = ii + sum(Numofprofile(1:dd-1));        % for handling multiple days of data
                elseif (dd == 1)
                    jj = ii;
                end
                nameofprofile(jj,:) = fscanf(fp,'%s',[1]);      % read profile name listed in the "filelist" or "Goodfilelist"
            end
            fclose(fp);

            for mm = 1:Numofprofile(dd)
                if (dd > 1)
                    nn = mm + sum(Numofprofile(1:dd-1));
                elseif (dd == 1)
                    nn = mm;
                end
                fpin = fopen([sprintf('/Users/dacostalindo/Desktop/Research/SkyObservations/%04d%s%02d/RX%d/%s',year,smonth(month,:),day(dd),system,nameofprofile(nn,:))],'r');
                fgetl(fpin); fgetl(fpin); fgetl(fpin); fgetl(fpin); fgetl(fpin);
                time(nn) = fscanf(fpin,'%f', [1]) + (day(dd)-day(1))*24; 
                fgetl(fpin); fgetl(fpin); fgetl(fpin); fgetl(fpin); 
                BinTime(nn) = fscanf(fpin,'%f', [1]);fgetl(fpin);
                binrange = c*BinTime(nn)*1e-9/2;
                BinWidth = round(binrange);
                fgetl(fpin);
                ShotsNum(nn) = fscanf(fpin,'%f', [1]);fgetl(fpin);
                fgetl(fpin);
                BinNum = fscanf(fpin,'%f', [1]); BinNum = round(BinNum);fgetl(fpin);
                fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);
                fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);fgetl(fpin);
                RawPhotonCnt(nn,:) = fscanf(fpin,'%f',[1 BinNum]);
                fclose(fpin);
                DataInfo(nn,:) = [year month day(dd) time(nn) ShotsNum(nn) BinTime(nn)];
                [Normdsumsig(nn,:),sumsig(nn,:),RayNormSig(nn,:),M(nn,:),sumsmdM(nn,:),sumdM(nn,:),B(nn,:),RayRemoval(nn,:)] = Feunitprocessupdated(DataInfo(nn,:),RawPhotonCnt(nn,:),HammingFWHM,system);
                if (PlotRayFigure_yn == 1), pause(0.1); end
                %% Compute MSISE00 data using MatLab built-in "atmosnrlmsise00.m" model
                timecurrent = time(nn); daycurrent = DataInfo(nn,3);
                ATPND = msise00cal(year, month, daycurrent, timecurrent, latitude, longitude, fluxindex);
                ATPND = ATPND';
                
                %% Find the temperature, pressure and number density at 50km
                ATPND_R = ATPND(find(ATPND(:,1)==round(Rayleigh_norm_alt)),1:4);
                AR = ATPND_R(1);
                TR = ATPND_R(2);		% Temperature at Rayleigh Normalization altitude
                PR = ATPND_R(3);		% Pressure at Rayleigh Normalization altitude
                NDR = ATPND_R(4);		% Number density at Rayleigh Normalization altitude
                
                %% Use interpolation to obtain the atmosphere number density at each bin altitude
                ND = interp1(ATPND(:,1),ATPND(:,4),(1:BinNum)*binrange);
                TP = interp1(ATPND(:,1),ATPND(:,2),(1:BinNum)*binrange);
                PuresigFe(nn,:) = Normdsumsig(nn,:) - ND/NDR;
                end
        end
        if (system == 1)
            DataInfo1 = DataInfo;
            RawData1 = RawPhotonCnt;
            PuresigFe374 = PuresigFe;
            Normdsumsig374 = Normdsumsig;
        elseif (system == 2)
            DataInfo2 = DataInfo;
            RawData2 = RawPhotonCnt;
            PuresigFe372 = PuresigFe;
            Normdsumsig372 = Normdsumsig;
        end
        clear DataInfo RawPhotonCnt time ShotsNum BinTime nameofprofile Numofprofile
        clear PuresigFe Normdsumsig sumsig RayNormSig M sumsmdM sumdM B RayRemoval
    end
	altitude = (1:BinNum)*binrange/1000;
    if (SavePreprocessedData == 1)
        clear datatobesaved
        Fe_start_bin = round(Fe_start_alt/binrange);
        Fe_end_bin = round(Fe_end_alt/binrange);
        datatobesaved = [355 DataInfo1(:,4)';altitude(Fe_start_bin:Fe_end_bin)' PuresigFe374(:,Fe_start_bin:Fe_end_bin)'];
        filename=sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/ProcessedData/%sPreproFeCnt374nm%04d%s%02dNew.txt',station,station,year,smonth(month,:),day(1));          
        save(filename,'datatobesaved','-ASCII')
        clear datatobesaved
        datatobesaved = [355 DataInfo2(:,4)';altitude(Fe_start_bin:Fe_end_bin)' PuresigFe372(:,Fe_start_bin:Fe_end_bin)'];
        filename=sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/ProcessedData/%sPreproFeCnt372nm%04d%s%02dNew.txt',station,station,year,smonth(month,:),day(1));
        save(filename,'datatobesaved','-ASCII')
        clear datatobesaved
        datatobesaved = [355 DataInfo1(:,4)';altitude(Fe_start_bin:Fe_end_bin)' Normdsumsig374(:,Fe_start_bin:Fe_end_bin)'];
        filename=sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/ProcessedData/%sPreproNormCnt374nm%04d%s%02dNew.txt',station,station,year,smonth(month,:),day(1));          
        save(filename,'datatobesaved','-ASCII')
        clear datatobesaved
        datatobesaved = [355 DataInfo2(:,4)';altitude(Fe_start_bin:Fe_end_bin)' Normdsumsig372(:,Fe_start_bin:Fe_end_bin)'];
        filename=sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/ProcessedData/%sPreproNormCnt372nm%04d%s%02dNew.txt',station,station,year,smonth(month,:),day(1));
        save(filename,'datatobesaved','-ASCII')
        clear datatobesaved
    end

    %% Divide data into unit time interval
    time1 = DataInfo1(:,4);
    time2 = DataInfo2(:,4);

    tmin = min(min(floor(time1)),min(floor(time2)));    % minimum UT hour
    tmax = max(max(ceil(time1)),max(ceil(time2)));      % maximum UT hour
    tmin = min(floor(time2));    % minimum UT hour
    tmax = max(ceil(time2));      % maximum UT hour
    if (ResolutionTime == 1)
        timebinsmooth = [tmin:0.05:tmax];  % unit: hour (Resolution in time is 0.25 hour)
    elseif (ResolutionTime == 2)
        timebinsmooth = [tmin:0.1:tmax];    % unit: hour (Resolution in time is 0.25 hour)
    elseif (ResolutionTime == 3)
        timebinsmooth = [tmin:0.2:tmax];    % unit: hour (Resolution in time is 0.25 hour)
    elseif (ResolutionTime == 4)
        timebinsmooth = [tmin:0.5:tmax];    % unit: hour (Resolution in time is 0.25 hour)
    elseif (ResolutionTime == 5)
        timebinsmooth = [tmin:1:tmax];    % unit: hour (Resolution in time is 0.25 hour)
    end

    for kk = 1:length(altitude)
        [smdPureFe374(:,kk)] = HammingSmth(FWtime,timebinsmooth,time1,PuresigFe374(:,kk));
        [smdPureFe372(:,kk)] = HammingSmth(FWtime,timebinsmooth,time2,PuresigFe372(:,kk));
    end
    
    altmin = floor(min(altitude));
    altmax = ceil(max(altitude));
    if (ResolutionAlt == 1)
        altbinsmooth = altitude;  % unit: km
    elseif (ResolutionAlt == 2)
        altbinsmooth = [altmin:0.1:altmax];    % unit: km
    elseif (ResolutionAlt == 3)
        altbinsmooth = [altmin:0.25:altmax];    % unit: km
    elseif (ResolutionAlt == 4)
        altbinsmooth = [altmin:0.5:altmax];    % unit: km
    elseif (ResolutionAlt == 5)
        altbinsmooth = [altmin:1:altmax];    % unit: km
    end

    for tt = 1:length(timebinsmooth)
        [smdsmdPureFe374(tt,:)] = HammingSmth(FWalt,altbinsmooth,altitude,smdPureFe374(tt,:));
        [smdsmdPureFe372(tt,:)] = HammingSmth(FWalt,altbinsmooth,altitude,smdPureFe372(tt,:));
    end
    
    if (SaveSmoothedData == 1)
        clear datatobesaved
        Fe_start_bin = round(Fe_start_alt/binrange);
        Fe_end_bin = round(Fe_end_alt/binrange);
        [NUM, Fe_start_bin] = min(abs(altbinsmooth - Fe_start_alt/1000));
        [NUM, Fe_end_bin] = min(abs(altbinsmooth - Fe_end_alt/1000));
        datatobesaved = [355 timebinsmooth;altbinsmooth(Fe_start_bin:Fe_end_bin)' smdsmdPureFe374(:,Fe_start_bin:Fe_end_bin)'];
        filename=sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/ProcessedData/%sSmoothFeCnt374nm%04d%s%02dNew.txt',station,station,year,smonth(month,:),day(1));          
        save(filename,'datatobesaved','-ASCII')
        clear datatobesaved
        datatobesaved = [355 timebinsmooth;altbinsmooth(Fe_start_bin:Fe_end_bin)' smdsmdPureFe372(:,Fe_start_bin:Fe_end_bin)'];
        filename=sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/ProcessedData/%sSmoothFeCnt372nm%04d%s%02dNew.txt',station,station,year,smonth(month,:),day(1));
        save(filename,'datatobesaved','-ASCII')
        clear datatobesaved
    end
    hold off
    %% Plot 2-D smoothed data and save data and figure %%%
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 70 scrsz(3)*0.5 scrsz(4)*0.3])
    orient portrait
    set(gcf,'PaperPositionMode','auto')
    pcolor(timebinsmooth, altbinsmooth,smdsmdPureFe372')
    shading interp; colorbar
    axis([tmin tmax 75 115])
    set(gca,'xtick',[tmin:2:tmax],'fontsize',16)
    set(gca,'tickdir','out')
    % caxis([-1000 15000])
    xlabel('Time UT (hr)','fontsize',18)
    ylabel('Altitude (km)','fontsize',18)
    % title('372-nm Fe Layers on 15-18 January 2011 (UT) ','fontsize',18)
    caxis([-1 10])
    if (length(day) == 1)
        title(sprintf('372-nm Fe Layers on %d %s %04d (UT) ',day,smonth(month,:),year),'fontsize',18)
    elseif (length(day) > 1)
        title(sprintf('372-nm Fe Layers on %d-%d %s %04d (UT) ',day(1), day(end),smonth(month,:),year),'fontsize',18)
    end
    saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/GoodfilesFePlots/FIG/%sFe372nm%04d%s%02dGoodNew.fig',station,station,year,smonth(month,:),day(1)));
    print('-depsc2',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/GoodfilesFePlots/EPS/%sFe372nm%04d%s%02dGoodNew',station,station,year,smonth(month,:),day(1)));
    print('-djpeg',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/GoodfilesFePlots/JPG/%sFe372nm%04d%s%02dGoodNew',station,station,year,smonth(month,:),day(1)));
    
    
    
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 70 scrsz(3)*0.5 scrsz(4)*0.3])
    orient portrait
    set(gcf,'PaperPositionMode','auto')
    pcolor(timebinsmooth, altbinsmooth,smdsmdPureFe374')
    shading interp; colorbar
    axis([tmin tmax 75 115])
    set(gca,'xtick',[tmin:2:tmax],'fontsize',16)
    set(gca,'tickdir','out')
    % caxis([-1000 15000])
    xlabel('Time UT (hr)','fontsize',18)
    ylabel('Altitude (km)','fontsize',18)
    % title('372-nm Fe Layers on 15-18 January 2011 (UT) ','fontsize',18)
    caxis([-0.05 0.5])
    if (length(day) == 1)
        title(sprintf('374-nm Fe Layers on %d %s %04d (UT) ',day,smonth(month,:),year),'fontsize',18)
    elseif (length(day) > 1)
        title(sprintf('374-nm Fe Layers on %d-%d %s %04d (UT) ',day(1), day(end),smonth(month,:),year),'fontsize',18)
    end
    saveas(gcf,sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/GoodfilesFePlots/FIG/%sFe374nm%04d%s%02dGoodNew.fig',station,station,year,smonth(month,:),day(1)));
    print('-depsc2',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/GoodfilesFePlots/EPS/%sFe374nm%04d%s%02dGoodNew',station,station,year,smonth(month,:),day(1)));
    print('-djpeg',sprintf('/Users/dacostalindo/Desktop/Research/Results/%s/OverviewStatistics/GoodfilesFePlots/JPG/%sFe374nm%04d%s%02dGoodNew',station,station,year,smonth(month,:),day(1)));

end


end


%% Process raw data for daily mean temperature or long-integration temperature (pushing the temperature boundaries)