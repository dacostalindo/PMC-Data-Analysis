function [processedPhotonCnt,znormphotoncnt,TotalPhotonCnt,FePhotonCnt,BkgPhotonCnt] = profileprocessOverview(photoncnt,DataInfo,timebin)

%% This subroutine is to process each frequency data profile for Fe Boltzmann Lidar.
% The process includes (1) PMT/discriminator saturation correction, (2) chopper correction, 
% (3) tilted background subtraction, (4) range-dependence removal, (5) base-altitude adjustment, 
% (6) Rayleigh normalization, and (7) Rayleigh signal removal.

%% constants and parameters (declare global variables)
global pi c Me Qe E0 Oss kB NA AMU Radius freqi strength
global FreqNaD2a FreqNaCrossOver FreqAOshift Lambda_center SigmaRay
global Deadtime_Disc No_of_Freqs PDA_Freq_Offset SigmaL MAX_Countrate
global Fe_start_alt Fe_end_alt vert_bin_res Rayleigh_sum_start_alt Rayleigh_sum_end_alt Rayleigh_norm_alt
global Rayleigh_fit_start_alt Rayleigh_fit_end_alt Rayleigh_fit_alt BG_start_alt BG_end_alt BG_start_alt2 BG_end_alt2
global temp_default wind_default density_default time_adjust max_temp min_temp max_wind minimum_signal Hamming_width
global smonth zulumonth
global PMTsaturation_correction_yn Chopper_correction_yn RayleighNormChoice
global latitude longitude fluxindex

%% profile information transferred from setprocess.m
year = DataInfo(1);
month = DataInfo(2);
day = DataInfo(3);
time = DataInfo(4);
azimuth = DataInfo(5);
elevation = DataInfo(6);
baseAlt = DataInfo(7);
binnum = DataInfo(8);
shotsnum = DataInfo(9);
binwid = DataInfo(10);
binrange = DataInfo(11);

year2 = year-2000;
zuluday = sum(zulumonth(1:month-1))+day;	% convert to zulu day (1-365)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Start to Process Data %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bino4baseAlt=round(baseAlt/sin(elevation*pi/180)/binrange);        % bin number corresponding to the base altitude
TotalPhotonCnt=photoncnt;       % total photon count before background subtraction for error analysis
TotalPhotonCnt(1:bino4baseAlt)=0;
TotalPhotonCnt(bino4baseAlt+1:binnum)=photoncnt(1:binnum-bino4baseAlt);     % adjustment to add base altitude
    
% % Subtract background
%     % To remove the tilted background, a linear fit is performed in the range of 5700-7200 bins 
%     % where chopper is fully open while no Na and Rayleigh signals present anymore,
%     % then extrapolate the fitting into the Na range and the Rayleigh range to estimate the background in these ranges
%     % and then subtract the estimated background from the Na and Rayleigh counts
%     % 3 frequency profiles have different tilted background due to different laser power causing different saturation of PMT,
%     % so the linear fit needs to perform for each frequency profiles.
% % linear fit to estimate the tilted background
% %xdata=[5700:7200];
% xdata=[5700/vert_bin_res:7200/vert_bin_res];
% %ydata=photoncnt(5700/vert_bin_res:7200/vert_bin_res);         % linear decay?
% %ydata=log(photoncnt(5700:7200));       % exponential decay?
% ydata=log(photoncnt(5700/vert_bin_res:7200/vert_bin_res));       % exponential decay?
% f=inline('x(1)+x(2)*xdata','x','xdata');
% x0(1)=mean(ydata);
% x0(2)=0;
% x=lsqcurvefit(f,x0,xdata/1000,ydata);         % Should we fit to the linear scale photon counts or the log scale photon counts?
%             % Here, dividing xdata by 1000 is to improve the precision of
%             % linear fit to the log-scale data to avoid -0.00000 slope
% sprintf('%.9f     %.9f',x)
% %bkg=f(x,[1:binnum]/1000);          % linear decay?
% bkg=exp(f(x,[1:binnum]/1000));          % exponential decay?
% %bkg(1:binnum)=mean(ydata);
%     %hold on
%     %plot([1:binnum],log(bkg),'g')
%     %pause(0.5)
%     %axis([2000 8000 0 100])
%% subtract a constant background
BG_start_bin = round(BG_start_alt/binrange);
BG_end_bin = round(BG_end_alt/binrange);
bkg(1:binnum) = mean(photoncnt(BG_start_bin:BG_end_bin));
photoncnt = photoncnt - bkg;

% subtract a constant background
% figure
% plot([1:binnum],photoncnt)
% axis([2000 8000 -10 100])
% axis([0 400 -100 1000])
% hold on
% plot([1,binnum],[0,0],'r')
% title('Photon count after background subtraction'); pause(0.2)

% save the Fe raw photon counts and background counts for error calculation (before removing range-dependence)
%bino4baseAlt=round(baseAlt*1000/cos(zenith*pi/180)/binrange);        % bin number corresponding to the base altitude
FePhotonCnt = photoncnt;
FePhotonCnt(1:bino4baseAlt) = 0;
FePhotonCnt(bino4baseAlt+1:binnum) = photoncnt(1:binnum-bino4baseAlt);
BkgPhotonCnt = bkg;
BkgPhotonCnt(1:bino4baseAlt) = 0;
BkgPhotonCnt(bino4baseAlt+1:binnum) = bkg(1:binnum-bino4baseAlt);

%% Remove the range-dependence
photoncnt = photoncnt.*((((1:binnum)-1/2).*binwid).^2);
        % use the center of each bin to denote its altitude

%% Base-altitude Adjustment (add base altitude)
bino4baseAlt = round(baseAlt/sin(elevation*pi/180)/binrange);
        % bin number corresponding to the base altitude
zphotoncnt = photoncnt;
zphotoncnt(1:bino4baseAlt) = 0;
zphotoncnt(bino4baseAlt+1:binnum) = photoncnt(1:binnum-bino4baseAlt);
% figure
% plot([1:binnum]*binrange/1000*cos(zenith*pi/180),log(max(zphotoncnt,1)))
% plot([1:binnum]*binrange/1000*cos(zenith*pi/180),zphotoncnt)
% title('Photon count after removal of range-dependence and adding base altitude')

%% Take Rayleigh normalization signal
bino4RayStart = round(Rayleigh_sum_start_alt/sin(elevation*pi/180)/binrange);
bino4RayEnd = round(Rayleigh_sum_end_alt/sin(elevation*pi/180)/binrange);
bino4RayNorm = round(Rayleigh_norm_alt/sin(elevation*pi/180)/binrange);
if (RayleighNormChoice == 1)
    PhotonCountRayNorm = exp(mean(log(zphotoncnt(bino4RayStart:bino4RayEnd))));
elseif (RayleighNormChoice == 2)
    PhotonCountRayNorm = mean(zphotoncnt(bino4RayStart:bino4RayEnd));
   % Noise = sqrt(abs(zphotoncnt(bino4RayStart:bino4RayEnd)));
end

% Take Rayleigh fit to get Rayleigh normalization signal
%xxdata=[bino4RayStart:bino4RayEnd];
%yydata=log10(zphotoncnt(bino4RayStart:bino4RayEnd));       % exponential decay?
%ff=inline('x(1)+x(2)*xxdata','xx','xxdata');
%xx0(1)=mean(yydata);
%xx0(2)=0;
%xx=lsqcurvefit(f,xx0,xxdata,yydata);         % Should we fit to the linear scale photon counts or the log scale photon counts?
%sprintf('%.9f     %.9f',xx)
%PCRayNorm=10^(f(xx,bino4RayNorm));          % exponential decay?
    %hold on
    %plot(Rayleigh_norm_alt/1000,log(PhotonCountRayNorm),'or',Rayleigh_norm_alt/1000,log(PCRayNorm),'sg')
    %axis([0 150 10 18])
    %plot(Rayleigh_norm_alt/1000,PhotonCountRayNorm,'or',Rayleigh_norm_alt/1000,PCRayNorm,'sg')
    %axis([20 80 1e5 3e5])

%% Rayleigh normalization
znormphotoncnt = zphotoncnt/PhotonCountRayNorm;       % using Rayleigh sum
%znormphotoncnt=zphotoncnt/PCRayNorm;                   % using Rayleigh fit

%% Read in MSIS-00 number density at deisred location
%Error
% fmsis = fopen([sprintf('/Users/chu/Tiger/ScienceProjects/Rothera/DataProcess/MSIS00/MSISE00zTPND103%03d',zuluday)],'r'); 		
% fscanf(fmsis,'%s %s	%s	%s',[4]);		% to skip
% ATPND = fscanf(fmsis,'%f	%f	%f	%f',[4 inf]);	% Altitude (km), temperature (K), pressure (mbar) and number density (cm^-3)
% ATPND = ATPND';
%fclose(fmsis);
ATPND = msise00cal(year, month, day, time, latitude, longitude, fluxindex);
ATPND = ATPND';
% ATPND = msise00cal(year, month, day, time, latitude, longitude, fluxindex);
% ATPND = ATPND';
MSIS00NumDen = ATPND(:,4);       % unit: cm^-3
xi=(1:binnum)*binrange/1000*sin(elevation*pi/180);        % the off-zenith angle must be considered here
% xi=(1:binnum)*binrange*sin(elevation*pi/180);        % the off-zenith angle must be considered here
yi=interp1(ATPND(:,1),MSIS00NumDen,xi,'spline');
% Find the number density at Rayleigh normalization altitude
NumDenRayNorm=interp1(ATPND(:,1),MSIS00NumDen,Rayleigh_norm_alt/1000,'spline');
% NumDenRayNorm=interp1(ATPND(:,1),MSIS00NumDen,Rayleigh_norm_alt,'spline');
    %figure
    %plot(log10(MSIS00NumDen),[0:1:250])
    %hold on
    %plot(log10(yi),xi,'r',log10(NumDenRayNorm),Rayleigh_norm_alt/1000,'og')

% Polynomial fit to Rayleigh + background range to estimate Rayleigh signals in the Na layer range
%bino4RayleighStart=round(35e3/cos(zenith*pi/180)/binrange);
%bino4RayleighEnd=round(75e3/cos(zenith*pi/180)/binrange);
%bino4BkgStart=round(130e3/cos(zenith*pi/180)/binrange);
%bino4BkgEnd=round(150e3/cos(zenith*pi/180)/binrange);
%RayBkgNormCnt(1:bino4RayleighEnd-bino4RayleighStart+1)=znormphotoncnt(bino4RayleighStart:bino4RayleighEnd);
%RayBkgNormCnt(bino4RayleighEnd-bino4RayleighStart+2:(bino4BkgEnd-bino4BkgStart+1)+(bino4RayleighEnd-bino4RayleighStart+1))=znormphotoncnt(bino4BkgStart:bino4BkgEnd);
%size(RayBkgNormCnt)
%RayBkg3Bino=[bino4RayleighStart:bino4RayleighEnd bino4BkgStart:bino4BkgEnd];
%p=polyfit(RayBkgBino,RayBkgNormCnt,4);
%RayBkgFit=polyval(p,[1:binnum]);
%processedPhotonCnt=znormphotoncnt-RayBkgFit;

% Subtract normalized Rayleigh signal
processedPhotonCnt = znormphotoncnt - yi/NumDenRayNorm;


% N_real = sum(processedPhotonCnt);
% N_noise = sum(sqrt(abs(photoncnt)));
% 
% SNR = N_real/N_noise;





% figure
% plot(xi,znormphotoncnt,'b',xi,yi/NumDenRayNorm,'r')
% plot(xi,znormphotoncnt,'b',xi,RayBkgFit,'r')
% axis([10 200 -1 12])
% title('Normalized photon count before subtracting Rayleigh')
% figure
% plot(xi,processedPhotonCnt,'c')
% hold on
% plot([0,250],[0,0],'-m')
% axis([10 200 -1 12])
% title('Normalized photon count after subtracting Rayleigh')