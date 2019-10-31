function varargout = interactive(varargin)
%Varargin: 
%1-year, 2- month,3-day, 4-handle)

%% Check for empty arguments


% INTERACTIVE MATLAB code for interactive.fig
%      INTERACTIVE, by itself, creates a new INTERACTIVE or raises the existing
%      singleton*.
%
%      H = INTERACTIVE returns the handle to a new INTERACTIVE or the handle to
%      the existing singleton*.
%
%      INTERACTIVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERACTIVE.M with the given input arguments.
%
%      INTERACTIVE('Property','Value',...) creates a new INTERACTIVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before interactive_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to interactive_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help interactive

% Last Modified by GUIDE v2.5 14-Feb-2018 11:59:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @interactive_OpeningFcn, ...
                   'gui_OutputFcn',  @interactive_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before interactive is made visible.
function interactive_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to interactive (see VARARGIN)

% Choose default command line output for interactive
handles.output = zeros(1,2);
handles.fig1 = 0;
handles.fig2 = 0;
handles.PMCflag =0;
handles.contour=0;

%set static text


if nargin == 0 || isempty(varargin)
  errorMessage = sprintf('Error: you did not enter any input(s)!');
  uiwait(warndlg(errorMessage));
  return
end

%% Search PMC from two lidar channels on-resonance and Compute PMC parameters

%% Declare global variables
%% Fundamental Constants and Atomic Fe Constants
global pi c h Me Qe E0 kB AMU NA mFeMean mFeIsotope AbundanceFe DeltaE_diff 
global Aki_Fe372 gk_Fe372 gi_Fe372 RB_Fe372 Lambda_Center_Fe372 fosc_Fe372 IsotopeShift_Fe372
global Aki_Fe374 gk_Fe374 gi_Fe374 RB_Fe374 Lambda_Center_Fe374 fosc_Fe374 

%% Linear PAL laser information (central wavelength, peak freq offset, linewidth, energy portion/pedestal, ...)
global PAL374Eportion PAL372Eportion PAL374RMSwidth PAL372RMSwidth 
global PAL374Detune PAL372Detune PAL374WL PAL372WL dEoverKB

%% Declare Global Variables for OK button, Calculate Button
global PMCl PMCh Rawbinwid AltCorrection R dR Beta_PMC dBeta year smonth month day filenum timebin time374 time372 factorPMC factorNoise PMTsaturation_correction_yn
global Chopper_correction_yn Rmax ZRpk Beta_max ZBpk Beta_Total dBeta_Total SigmaRMS FWHM PMCZcl PMCZch Zpkbin Zc Altl Alth
%% Constants and parameters (declare global variables)
global Fe_start_alt Fe_end_alt vert_bin_res Rayleigh_sum_start_alt Rayleigh_sum_end_alt Rayleigh_norm_alt
global Rayleigh_fit_start_alt Rayleigh_fit_end_alt Rayleigh_fit_alt BG_start_alt BG_end_alt BG_start_alt2 BG_end_alt2 
global temp_default wind_default density_default time_adjust max_temp min_temp max_wind minimum_signal Hamming_width
global  zulumonth
global RayleighNormChoice
global SavePMCData SaveFeData PlotFigure SaveFigure PlotRayFigure PlotRayFigure_yn
global BaseAlt latitude longitude fluxindex checkE checkEE

%% assign name to inputs varargin
Info1 = varargin{1};
RawCnt1 = varargin{2};
Info2 = varargin{3};
RawCnt2 = varargin{4};
timebin = varargin{5};
HammingFWHM = varargin{6};
filenum = varargin{7};


% UIWAIT makes interactive wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%% Data information transferred from mainprocess.m and setprocess.m
year = Info1(1,1);
month = Info1(1,2);
day = Info1(1,3);
time374 = Info1(:,4);
time372 = Info2(:,4);
ShotsNum = Info1(:,5);
BinTime = Info1(1,6);
Rawbinwid = BinTime/2*1e-9*c;       % Bin width (unit: m), considering round-trip
NumofProfile = size(RawCnt1,1);     % number of profiles in the unit time interval
BinNum = size(RawCnt1,2);           % total bin number in each profile
zuluday = sum(zulumonth(1:month-1)) + day;	% convert to zulu day (1-365)

showUTimeProcessing = timebin;

%% Bin Number Corresponding to Absolute Altitude for
% RayStart, RayEnd, FeStart, FeEnd, BgStart, BgEnd and Base Altitude
bino4RayStart = round(Rayleigh_fit_start_alt/Rawbinwid);
bino4RayEnd = round(Rayleigh_fit_end_alt/Rawbinwid);
bino4RayNorm = round(Rayleigh_norm_alt/Rawbinwid);
bino4FeStart = round(Fe_start_alt/Rawbinwid);
bino4FeEnd = round(Fe_end_alt/Rawbinwid);
bino4BgStart = round(BG_start_alt/Rawbinwid);
bino4BgEnd = round(BG_end_alt/Rawbinwid);
bino4BaseAlt = round(BaseAlt/Rawbinwid);

%% Standard pre-process of raw photon counts for both channels
system = 1;  %figure;
[Normdsumsig374,sumsig374,RayNormSig374,M,sumsmdM,sumdM,B,RayRemoval374] = PMCunitprocess(Info1,RawCnt1,HammingFWHM,system);
system = 2; 
             %figure;
[Normdsumsig372,sumsig372,RayNormSig372,MM,sumsmdMM,sumdMM,BB,RayRemoval372] = PMCunitprocess(Info2,RawCnt2,HammingFWHM,system);
% figure; plot((1:BinNum)*Rawbinwid/1000,Normdsumsig374)

%% Load atmosphere number density, temperature, pressure vs altitude database
% fmsis = fopen([sprintf('/Users/chu/Tiger/ScienceProjects/Rothera/DataProcess/MSIS00/MSISE00zTPND103%03d',zuluday)],'r'); 		
% fscanf(fmsis,'%s %s	%s	%s',[4]);		% to skip
% ATPND = fscanf(fmsis,'%f	%f	%f	%f',[4 inf]);	% Altitude (km), temperature (K), pressure (mbar) and number density (cm^-3)
% ATPND = ATPND';
% fclose(fmsis);
% %% Find the temperature, pressure and number density at Rayleigh normalization altitude 50km
% ATPND_R = ATPND(find(ATPND(:,1) == round(Rayleigh_norm_alt/1000)),1:4);
% AR = ATPND_R(1);
% TR = ATPND_R(2);		% Temperature at Rayleigh Normalization altitude
% PR = ATPND_R(3);		% Pressure at Rayleigh Normalization altitude
% NDR = ATPND_R(4);		% Number density at Rayleigh Normalization altitude
% %% Use interpolation to obtain the atmosphere number density at each bin altitude
% NDmsis = interp1(ATPND(:,1),ATPND(:,4),(1:BinNum)*Rawbinwid/1000);
% Tmsis = interp1(ATPND(:,1),ATPND(:,2),(1:BinNum)*Rawbinwid/1000);
% figure; plot(log(NDmsis),(1:BinNum)*Rawbinwid/1000,log(NDR),Rayleigh_norm_alt/1000,'o')

%% Compute the Beta_R(z_R) at z_R=50km, i.e., Rayleigh backscatter coefficient at normalization altitude 50 km
%% Compute MSISE00 data using MatLab built-in "atmosnrlmsise00.m" model
ATPND = msise00cal(year, month, day, timebin, latitude, longitude, fluxindex);
ATPND = ATPND';

%% Find the temperature, pressure and number density at Z_R (50 km)
ATPND_R = ATPND(find(ATPND(:,1) == round(Rayleigh_norm_alt)),1:4);
AR = ATPND_R(1);
TR = ATPND_R(2);		% Temperature at Rayleigh Normalization altitude
PR = ATPND_R(3)/100;	% Pressure at Rayleigh Normalization altitude (unit: mbar)
NDR = ATPND_R(4);		% Number density at Rayleigh Normalization altitude
%% Use interpolation to obtain the atmosphere number density at each bin altitude
NDmsis = interp1(ATPND(:,1),ATPND(:,4),(1:BinNum)*Rawbinwid);
Tmsis = interp1(ATPND(:,1),ATPND(:,2),(1:BinNum)*Rawbinwid);

%% Compute the Beta_R(z_R) using the model temperature, pressure, and number density
% Beta_RNm=(0.7629*(1+0.9324)/(4*pi)*9.807*1e-23*273*PR/(TR*1013*(lambda(system)*100)^4.0117));     % for wavelength in cm
Beta_RNm(1) = (0.7629*(1+0.9324)/(4*pi)*9.2925*1e-31*273*PR/(TR*1013*Lambda_Center_Fe374^4.0117));	% for wavelength in meter (374nm channel)
Beta_RNm(2) = (0.7629*(1+0.9324)/(4*pi)*9.2925*1e-31*273*PR/(TR*1013*Lambda_Center_Fe372^4.0117));	% for wavelength in meter (372nm channel)

%% subtract Rayleigh signal from normalized photon count
Normdpuresig374 = Normdsumsig374 - NDmsis/NDR;
 %figure; plot((1:BinNum)*Rawbinwid/1000,Normdpuresig374); pause(0.1)
% title('Normalized photo count for 374 nm wavelength')
 %xlabel('Altitude [km]')
% ylabel('photo count')
%% Old parameters for PAL laser information (central wavelength, peak freq offset, linewidth, energy portion/pedestal, ...)
% lambda=[373.8195e-9,372.0995e-9];	% laser wavelength in meter unit
% dEoverKB = h*c*100*DeltaE/kB;   % Delta_E/kB for the energy level splitting
% % dEoverKB=598.438;
% 
% T0=[176,168];						% Assumed temperature for 374 nm scan on 30 Aug 2002, 372 nm scan on 14 Aug 2002
% SigmaD0=[433e6,425e6];				% Doppler Broadening for temperature T0
% SigmaL_1=[437.0e6,435.3e6];			% Laser Linewidth (narrow peak)
% SigmaL_2=[14212e6,14212e6];			% Laser Linewidth (wide peak)
% %A=[0.448,0.552];					% energy portion of narrow peak and wide peak for 372nm and old 374nm laser
% %AA=[0.75,0.25];						% energy portion of narrow peak and wide peak for new 374nm laser
% % to read in the AA and A factors for each date at Rothera
% fin=fopen([sprintf('/Users/chu/Tiger/ScienceProjects/Rothera/DataProcess/RotheraPMC/RotheraAAfactor.dat')],'r');
% fscanf(fin,'%s	%s	%s	%s	%s	%s	%s',[7]);
% RotheraAA=fscanf(fin, '%d	%d	%d	%f	%f	%f	%f', [7 inf]);
% RotheraAA=RotheraAA';
% fclose(fin);
% AAmdy=RotheraAA(:,1:3);
% AA=RotheraAA(find(AAmdy(:,1)==month & AAmdy(:,2)==day & AAmdy(:,3)==year),4:5)  % for 374nm channel 
% A=RotheraAA(find(AAmdy(:,1)==month & AAmdy(:,2)==day & AAmdy(:,3)==year),6:7)   % for 372nm channel
% 
% sigmaDsimga0=[461.8e6*0.876e-16,463.8e6*0.943e-16];		
% 									% sigmaD(T=200K)*sigma0(T=200K)
% lambda_ref=[373.8195e-9,372.0995e-9];	% reference Wavelength [747.6390, 744.1990]/2*1e-9 (m)
% %lambda=[373.8196e-9,372.0996e-9];	% Laser Wavelength [747.6392, 744.1992]/2*1e-9 (m) used during taking data
% lambda=[373.81965e-9,372.0996e-9];	% Laser Wavelength [747.6393, 744.1992]/2*1e-9 (m) used during taking data
% %lambda=[373.8196e-9,372.09955e-9];	% Laser Wavelength [747.6392, 744.1991]/2*1e-9 (m) used during taking data
% 
% % the following dnu and dnu_real are obtained after wavemeter calibration at BAS, and should be used for data at BAS and at Rothera
% dnu_0=[-314.6e6,-176.5e6];			% Frequency Offset (Hz) for Wavelength Scan at BAS (obtained from Gaussian fit): 
% 									% i.e., Fe spectrum peak frequency relative to wavelegnth [747.6390, 744.1990]/2
% dnu_lambda=(c./lambda-c./lambda_ref);	
% 									% Frequency offset (Hz) of laser wavelength relative to reference wavelength
% dnu_real=dnu_lambda-dnu_0;			% Frequency offset (Hz) of laser wavelength relative to Fe spectrum peak freuqency

%% PAL laser central wavelength and frequency tuning
LaserWL372 = 744.1992e-9 / 2;
LaserFreqOffset372 = -c/Lambda_Center_Fe372^2 * (LaserWL372 - Lambda_Center_Fe372) + PAL372Detune;
LaserWL374 = 747.6394e-9 / 2;
LaserFreqOffset374 = -c/Lambda_Center_Fe374^2 * (LaserWL374 - Lambda_Center_Fe374) + PAL374Detune;

%% To compute Beta_Fe(z), Density_Fe(z)
Altl = round(75e3/Rawbinwid);         % lower limit bin for PMC
Alth = round(115e3/Rawbinwid);		% upper limit bin for PMC
% Altl = round(70e3/Rawbinwid);         % lower limit bin for PMC
% Alth = round(120e3/Rawbinwid);		% upper limit bin for PMC
R = zeros(1,BinNum);				% backscatter ratio
dR = zeros(1,BinNum);              % error of backscatter ratio R
Beta_PMC = zeros(1,BinNum);		% Volume backscatter coefficient of pure PMC
dBeta = zeros(1,BinNum);           % error of volume backscatter coefficient

E = 1;
EE = 1;
for kk = Altl:Alth
%     % to calculate SigmaEff(1) and SigmaEff(2) separately to keep the precision of digits
% 	SigmaD(1)=SigmaD0(1)*sqrt(TP(kk)/T0(1));	
% 	SigmaD(2)=SigmaD0(2)*sqrt(TP(kk)/T0(2));	
% 	SigmaRMS_1(1)=sqrt(SigmaL_1(1)^2+SigmaD(1)^2);
% 	SigmaRMS_1(2)=sqrt(SigmaL_1(2)^2+SigmaD(2)^2);
% 	SigmaRMS_2(1)=sqrt(SigmaL_2(1)^2+SigmaD(1)^2);
% 	SigmaRMS_2(2)=sqrt(SigmaL_2(2)^2+SigmaD(2)^2);
% 	SigmaEff(1)=A(1)*sigmaDsimga0(1)/SigmaRMS_1(1)*exp(-dnu_real(1)^2/(2*SigmaRMS_1(1)^2))+A(2)*sigmaDsimga0(1)/SigmaRMS_2(1)*exp(-dnu_real(1)^2/(2*SigmaRMS_2(1)^2));
% 	SigmaEffnew(1)=AA(1)*sigmaDsimga0(1)/SigmaRMS_1(1)*exp(-dnu_real(1)^2/(2*SigmaRMS_1(1)^2))+AA(2)*sigmaDsimga0(1)/SigmaRMS_2(1)*exp(-dnu_real(1)^2/(2*SigmaRMS_2(1)^2));
% 							% 374nm effective cross-section with new fraction of pedestal in laser spectrum
% 	SigmaEff(2)=A(1)*sigmaDsimga0(2)/SigmaRMS_1(2)*exp(-dnu_real(2)^2/(2*SigmaRMS_1(2)^2))+A(2)*sigmaDsimga0(2)/SigmaRMS_2(2)*exp(-dnu_real(2)^2/(2*SigmaRMS_2(2)^2));
%     % to compute the ratio of N372 over N374
% 	K(kk)=SigmaEff(2)/SigmaEffnew(1)/0.9114*9/7*exp(598.438/TP(kk));

    %% Calculate effective cross sections
    effCrossecion374(kk) = EffCrossectionFe374(LaserFreqOffset374,Tmsis(kk),0);
    effCrossecion372(kk) = EffCrossectionFe372(LaserFreqOffset372,Tmsis(kk),0);
    %% to compute the ratio of N372 over N374% pg274
	K(kk) = effCrossecion372(kk)/effCrossecion374(kk)*RB_Fe372/RB_Fe374*gi_Fe372/gi_Fe374*exp(dEoverKB/Tmsis(kk));
	%% Calculate Beta, R with the consideration of subtracting 374nm Fe signal from the 374nm photon counts
	Beta372(kk) = (Normdsumsig372(kk)-NDmsis(kk)/NDR)*Beta_RNm(2);
	Beta374(kk) = (Normdsumsig374(kk)-NDmsis(kk)/NDR)*Beta_RNm(1);
   
	Beta_Fe(kk) = (Normdsumsig372(kk)-Normdsumsig374(kk))*Beta_RNm(2)/(1-(Lambda_Center_Fe374/Lambda_Center_Fe372)^4/K(kk))/(E^2);
    
	CrtFe(kk) = 1/(K(kk)-1)*(Normdsumsig372(kk)-Normdsumsig374(kk));
	R(kk) = (Normdsumsig374(kk)-CrtFe(kk))*(NDR/NDmsis(kk));	% Backscattering ratio after subtracting Fe 374 signal
	Beta_PMC(kk) = (Normdsumsig374(kk)-NDmsis(kk)/NDR-CrtFe(kk))*Beta_RNm(1);
                % volume backscattering coefficient after subtracting Fe 374 signal
	%% Compute the pure photon counts (i.e., remove *R^2) for error analysis
	NFe374(kk-Altl+1) = CrtFe(kk)*RayNormSig374/(((kk-bino4BaseAlt)*Rawbinwid)^2);
	NR374(kk-Altl+1) = RayRemoval374(kk-Altl+1)/(((kk-bino4BaseAlt)*Rawbinwid)^2);
	%% Compute the equivalent Fe density for the corresponding Beta
	Density_Fe(kk) = 4*pi*Beta_Fe(kk)/RB_Fe372/effCrossecion372(kk)/(E^2);		% Unit: m^-3  for pure Fe
	DenFe372(kk) = 4*pi*Beta372(kk)/RB_Fe372/effCrossecion372(kk)/(E^2);		% Unit: m^-3  for 372 channel
	DenFe374(kk) = 4*pi*Beta374(kk)/1/effCrossecion374(kk)/(EE^2);	% Unit: m^-3  for 374 channel
			% the signals around 85 km in 374 channel in summer are mainly PMC (not Fe), it doesn't obey the branching ratio
            % so the branching ratio is set to 1
    %% Update the extinction coefficient
	E = E*exp(-effCrossecion372(kk)*Rawbinwid*Density_Fe(kk));		% extinction coefficient for 372nm due to Fe absorption
	EE = EE*exp(-effCrossecion374(kk)*Rawbinwid*DenFe374(kk));		% extinction coefficient for 374nm due to Fe absorption
end
checkE = E
checkEE = EE
% figure; plot((Altl:Alth),K(Altl:Alth))
Density_Fe = Density_Fe*1e-6;
DenFe372 = DenFe372*1e-6;
DenFe374 = DenFe374*1e-6;
FeInfo = 2;       % Fe density is derived from two channels

PMCl = Altl;
PMCh = Alth;

%% Calculate Backscatter ratio error and volume backscattering coefficient error from 374nm channel
Mp = M(PMCl:PMCh);				% total photon counts in PMC range
dMp = sumsmdM(PMCl:PMCh);       % photon noise in PMC range
Mr = M(bino4RayStart:bino4RayEnd);				% photon counts for Rayleigh normalization
%dMr=sumsmdM(Rayl:Rayh);		% smoothed photon noise in Rayleigh range (lower bound)
%dB=sumsmdM(Bgl:Bgh);			% smoothed photon noise in Background range (lower bound)
dMr = sumdM(bino4RayStart:bino4RayEnd);         % un-smoothed photon noise in Rayleigh range (upper-bound)
dB = sumdM(bino4BgStart:bino4BgEnd);			% un-smoothed photon noise in Background range (upper-bound)

%varP=(dMp).^2./((Mp-B).^2);
varP = (dMp).^2./((Mp-B-NR374).^2);
                                % photon noise due to PMC photon counts
varR = sum((dMr).^2./((Mr-B).^2)/((bino4RayEnd-bino4RayStart+1)^2));
                                % photon noise due to Rayleigh normalization photon counts
%varB=sum(dB.^2/((Bgh-Bgl+1)^2))*(sum(1./(Mr-B)/(Rayh-Rayl+1))-1./(Mp-B)).^2;
varB = sum(dB.^2/((bino4BgEnd-bino4BgStart+1)^2))*(sum(1./(Mr-B)/(bino4RayEnd-bino4RayStart+1))-1./(Mp-B-NR374)).^2;
                                % photon noise due to subtracted background
MpB = Mp-B;

RR = R(PMCl:PMCh);
dR(PMCl:PMCh) = RR.*sqrt(varP+varR+varB);		% absolute error of backscattering ratio
BetaB = Beta_PMC(PMCl:PMCh);
dBeta(PMCl:PMCh) = BetaB.*sqrt(varP+varR+varB);	% absolute error of volume backscattering coefficient

% Hamming Smooth
HammingSmthyn = 0;
if (HammingSmthyn == 1)
    [smdBeta372]=HammingSmooth(Beta372,HammingFWHM);
    [smdBeta374]=HammingSmooth(Beta374,HammingFWHM);
    [smdBeta_Fe]=HammingSmooth(Beta_Fe,HammingFWHM);
    [smdBeta_PMC]=HammingSmooth(Beta_PMC,HammingFWHM);
    [smdR]=HammingSmooth(R,HammingFWHM);
    %[smdDensity_Fe]=HammingSmooth(Density_Fe,HammingFWHM);
   % [smdDenFe372]=HammingSmooth(DenFe372,HammingFWHM);
   % [smdDenFe374]=HammingSmooth(DenFe374,HammingFWHM);
end

%% Altitude correction for base altitude non-integer number of Rawbinwid and center of each bin:
% (1) When convert altitude to bin-number, non-integer part of base altitude will be ignored or over-count;
% (2) For each bin, the altitude Z should be the center of the bin: Z=(binnumber-0.5)*Rawbinwid, instead of Z=binnumber*Rawbinwid.
AltCorrection = BaseAlt-bino4BaseAlt*Rawbinwid-0.5*Rawbinwid;
%AltCorrection = BaseAlt-(bino4BaseAlt*Rawbinwid-0.5*Rawbinwid);
				% Only absolute altitude needs correction (Z=Z+AltCorrection)
				% All rms width, FHWM and their errors don't need correction because of (Zi-Zc)
				% Error of Zc also doesn't need correction because of (Zi-Zc)
				% Modification made on September 22, 2000 around 23:00 LST

%% Use several critiria to judge whether the processed profile is qualified to be PMC
%% Search PMC peak in the range of 75-90km, and Beta_max has to be larger than zero
PMCZl = round(75e3/Rawbinwid);
% PMCZl = round(80e3/Rawbinwid);      % for 5UT on 2011/Feb/05
PMCZh = round(90e3/Rawbinwid);
% PMCZh = round(86e3/Rawbinwid);          % for 2 UT on 2011/Feb/12
%% Find the maximum peak in 75-90km range and its corresponding index
[Beta_max,ZpkIndex] = max(Beta_PMC(PMCZl:PMCZh));

Zpkbin = PMCZl-1+ZpkIndex;                  % bin number for the peak Beta
ZBpk = Rawbinwid*Zpkbin + AltCorrection;    % convert the index to altitude
%% !!!!!!!!!!!!!
dBeta(Zpkbin)                               % display for monitoring purpose
%% Judge whether Bmax is larger than 2*dBeta
if ((Beta_max - 2*abs(dBeta(Zpkbin))) > 0)
	BmaxLG2Err = 1;	% Bmax is larger than 2 times of error
else
	BmaxLG2Err = 0;
end

if (Beta_max > 0)
	BmaxLG0 = 1;
	FINDHfl = 0;
	m = 1;
	while (FINDHfl == 0 & m < (Zpkbin-PMCZl))		% to find the FWHM point at lower end
		if (Beta_PMC(Zpkbin-m) <= Beta_max/2 & Beta_PMC(Zpkbin-m+1) > Beta_max/2)
			ZHMlbin = Zpkbin-m;
			FINDHfl = 1;
		end
		m = m+1;
	end
	FINDHfh = 0;
	n = 1;	
	while (FINDHfh == 0 & n< (PMCZh-Zpkbin))		% to find the FWHM point at higher end
		if (Beta_PMC(Zpkbin+n-1) > Beta_max/2 & Beta_PMC(Zpkbin+n) <= Beta_max/2)
			ZHMhbin = Zpkbin+n;
			FINDHfh = 1;
		end
		n = n+1;
	end
	if (FINDHfl == 1 & FINDHfh == 1)
		FWHM = (ZHMhbin-ZHMlbin+1)*Rawbinwid;
		if (min(Beta_PMC(ZHMlbin:ZHMhbin)) > 0)
			BFWHMLG0 = 1;
			[Diff1,L1index] = min(Beta_PMC(ZHMlbin:ZHMhbin) - abs(dBeta(ZHMlbin:ZHMhbin)));
			if (Diff1>0)
				BLG1Err = 1;	% All Beta within FWHM range are larger than 1 time of error
			else
				BLG1Err = 0;
			end
            factorPMC = 1.5;  % change to 1.5 times of error (for most days)
%             factorPMC = 1.2;  % change to 1.2 times of error (for a few hours in 12/27/03, 2 hours in 12/12/03, 1 hour on 2010/12/21)
%             factorPMC = 1.05;  % change to 1.1 times of error (for a few hours in 12/10/04, 1 hour on 2010/12/26)
			%[Diff2,L2index]=min(Beta_PMC(ZHMlbin:ZHMhbin)-2*dBeta(ZHMlbin:ZHMhbin));   
			%[Diff2,L2index]=min(Beta_PMC(ZHMlbin:ZHMhbin)-1.5*dBeta(ZHMlbin:ZHMhbin));  % change to 1.5 times of error (for most days)
			[Diff2,L2index] = min(Beta_PMC(ZHMlbin:ZHMhbin)-factorPMC*dBeta(ZHMlbin:ZHMhbin));  % change to 1.2 times of error (for a few hours in 12/27/03, 2 hours in 12/12/03)
			%[Diff2,L2index]=min(Beta_PMC(ZHMlbin:ZHMhbin)-1.1*dBeta(ZHMlbin:ZHMhbin)); 
			if (Diff2 > 0)	% This standard is too high!!!
				BLG2Err = 1;	% All Beta within FWHM range are larger than 2 times of error
            else
				BLG2Err = 0;
            end
			
		else
			BFWHMLG0 = 0;		% all bins within FWHM range should be larger than zero
			sprintf('there are negative bins within FWHM,although Bmax > 0')
			[Diff1,L1index] = min(Beta_PMC(ZHMlbin:ZHMhbin)-abs(dBeta(ZHMlbin:ZHMhbin)));
			if (Diff1 > 0)
				BLG1Err = 1;	% All Beta within FWHM range are larger than 1 time of error
			else
				BLG1Err = 0;
			end
            factorPMC = 1.5;  % change to 1.5 times of error (for most days)
			[Diff2,L2index] = min(Beta_PMC(ZHMlbin:ZHMhbin)-factorPMC*dBeta(ZHMlbin:ZHMhbin));  % change to 1.2 times of error (for a few hours in 12/27/03, 2 hours in 12/12/03)
			if (Diff2 > 0)	% This standard is too high!!!
				BLG2Err = 1;	% All Beta within FWHM range are larger than 2 times of error
            else
				BLG2Err = 0;
            end
		end
	else
		sprintf('Cannot find FWHM points: FINDHfl, FINDHfh')
		BFWHMLG0 = NaN;
		BLG1Err = NaN;
        factorPMC = NaN;
	end
else
	BmaxLG0 = 0;		% Bmax should be larger than zero
	sprintf('All are negative photon counts (Bmax < 0), check background subtraction')
    PMCyn = NaN;
    FeInfo = NaN;
    return
end

%% Compute the standard deviation of noise (take range of 95-105km for now)
% change the range for background noise computation to 95-105 km
noiseZl = round(95e3/Rawbinwid);	% for most data
noiseZh = round(105e3/Rawbinwid);	% for most data
%noiseZl=round(90e3/Rawbinwid);	% for 29 Dec 2002 data
%noiseZh=round(100e3/Rawbinwid);	% for 29 Dec 2002 data
%noiseZl=round(90e3/Rawbinwid);	% for 19 Jan 2002 59UT data
%noiseZh=round(100e3/Rawbinwid);	% for 19 Jan 2002 59UT data
noisemean = mean(Beta_PMC(noiseZl:noiseZh));
stdnoise = std(Beta_PMC(noiseZl:noiseZh))		% divided by (N-1)
%factorNoise=2;
factorNoise = 1.5;
% factorNoise=1.2;        % for 2UT on 2011/02/12
if (Beta_max > factorNoise*stdnoise)		% define as 2 times of std deviation
	BmaxLGnoise=1;
else
	BmaxLGnoise=0;
end

%% Judge whether the peak is qualified to be a PMC
%if (BmaxLG0==1 & BFWHMLG0==1 & BLG2Err==1 & BmaxLGnoise==1)
%if (BmaxLG0==1 & BFWHMLG0==1 & BLG1Err==1 & BmaxLG2Err==1 & BmaxLGnoise==1)
%if (BmaxLG0==1 & BFWHMLG0==1 & BLG1Err==1 & BLG2Err==1 & BmaxLG2Err==1)    % remove stdnoise condition for 17UT on 12/27/2003, 15UT on 12/12/03
if (BmaxLG0==1 & BFWHMLG0==1 & BLG1Err==1 & BLG2Err==1 & BmaxLG2Err==1 & BmaxLGnoise==1)    % BLG2Err is actually for 1.5 times of error
	PMCyn = 1;
    handles.PMCflag = 1;
else
	PMCyn = 0;
end






%% Compute the PMC parameters after qualifying a peak as PMC
if (PMCyn==1)
	FIND0l=0;
	m=1;
	while (FIND0l==0)		% to find the zero point at lower end
		if (Beta_PMC(Zpkbin-m)<=0 & Beta_PMC(Zpkbin-m+1)>0)
			PMCZcl=Zpkbin-m;		% PMC layer start bin
			FIND0l=1
		end
		m=m+1;
	end
	FIND0h=0;
	n=1;	
	while (FIND0h==0)		% to find the zero point at higher end
		if (Beta_PMC(Zpkbin+n-1)>0 & Beta_PMC(Zpkbin+n)<=0)
			PMCZch=Zpkbin+n;		% PMC layer stop bin
			FIND0h=1
		end
		n=n+1;
    end
    
	if (FIND0l==1 & FIND0h==1)
		if ((PMCZch-Zpkbin)*Rawbinwid>3.5e3)	% if half side of PMC layer is larger than 3.5 km
			rangeHfPMC = round(3e3/Rawbinwid);			% 3.0km corresponding bin number		
			[Bminh,minhIndex] = min(Beta_PMC(Zpkbin:Zpkbin+rangeHfPMC));
			if (abs(Bminh) < Beta_max/5)		% if the minimumm less than 1/10th of Bmax, take this minimum as zero point for PMC stop
				PMCZch = Zpkbin-1+minhIndex;
				refinedPMCZch = PMCZch*Rawbinwid
            else
                PMCZch=Zpkbin+rangeHfPMC;
			end
        end
        
		if ((Zpkbin-PMCZcl)*Rawbinwid>4e3)
			rangeHfPMC=round(3.5e3/Rawbinwid)			% 3.5km corresponding bin number		
			[Bminl,minlIndex]=min(Beta_PMC(Zpkbin-rangeHfPMC:Zpkbin))
			if (abs(Bminl)<Beta_max/5)		% if the minimumm less than 1/10th of Bmax, take this minimum as zero point for PMC start
				PMCZcl=Zpkbin-rangeHfPMC-1+minlIndex;
            else
                PMCZcl=Zpkbin-rangeHfPMC;
			end
		end
	else
		sprintf('Cannot find zero point for PMC layer start or stop')
	end
	% after find zero points for both start and stop, compute PMC parameters
	if (FIND0l==1 & FIND0h==1)
		% centroid altitude (km)
		Zc = sum(Beta_PMC(PMCZcl:PMCZch).*(PMCZcl:PMCZch)*Rawbinwid)/sum(Beta_PMC(PMCZcl:PMCZch));
		% rms width (km)
		SigmaRMS = sqrt(sum(((PMCZcl:PMCZch)*Rawbinwid-Zc).^2.*Beta_PMC(PMCZcl:PMCZch))/sum(Beta_PMC(PMCZcl:PMCZch)));
		% total backscatter coefficient
		Beta_Total = sum(Beta_PMC(PMCZcl:PMCZch).*Rawbinwid);
        dBeta_Total = sum(dBeta(PMCZcl:PMCZch))*1e6;
		% maximum backscatter ratio of PMC layer
		[Rmax,RmaxIndex] = max(R(PMCZcl:PMCZch));
		ZRpkbin = PMCZcl-1+RmaxIndex;
		ZRpk = ZRpkbin*Rawbinwid + AltCorrection;
        Zc = Zc + AltCorrection;
	end
end

%% Plot PMC and photon count signals
%%
time374;
DenFe374;
(Altl:Alth)*Rawbinwid/1000;



% figure;
% hold on;
% colormap('jet');
% contourf(time374,(Altl:Alth)*Rawbinwid/1000,DenFe374(Altl:Alth))
% colorbar;
% title('Fe Density 374 nm')
% xlabel('time[s]');
% ylabel('Altitude [km]');
% hold off;
% 
% figure;
% hold on;
% colormap('jet');
% contourf(time372,(Altl:Alth)*Rawbinwid/1000,DenFe372(Altl:Alth))
% colorbar;
% xlabel('time [s]');
% ylabel('Altitude [km]');
% title('Fe Density 372 nm')
% hold off;
%% Save beta(z) data no matter there are PMC or not (for future contour purpose)

    AltRBeta = [((PMCl:PMCh)*Rawbinwid+AltCorrection)/1000;R(PMCl:PMCh);dR(PMCl:PMCh);1e9*Beta_PMC(PMCl:PMCh);1e9*dBeta(PMCl:PMCh);1e9*Beta374(PMCl:PMCh);1e9*Beta372(PMCl:PMCh)];	% Altitude, R, dR, BetaPMC, dBeta, Beta374, Beta372
	fpout=fopen(sprintf('/Users/dacostalindo/Desktop/Research/Results/McMurdo/ProcessedData/betaData/Beta%04d%s%02d_%03d.dat',year,smonth(month,:),day,filenum),'w');
	fprintf(fpout,'NumChannel	Year	Month	Date	UTHour	StartTime374(UT)	StopTime374(UT)    factorPMC   factorNoise\n');
	fprintf(fpout,'%d	%d	%02d	%02d	%02d	%4.2f	%4.2f   %3.2f   %3.2f\n',2,year,month,day,timebin,time374(1),time374(end),factorPMC,factorNoise);
	fprintf(fpout,'Altitude(km)	BSRatio	BSRErr	BetaPMC(10^-9 m-1sr-1)	BetaErr Beta374 Beta372\n');
	fprintf(fpout,'%7.3f	%7.3f	%7.3f	%7.3f	%7.3f   %7.3f   %7.3f\n',AltRBeta);
	fclose(fpout);

%% Save Fe density data (pure Fe, 372nm Fe, and 374nm Fe) for future contour purpose

    FeDensityAll = [((PMCl:PMCh)*Rawbinwid+AltCorrection)/1000;Density_Fe(PMCl:PMCh);DenFe374(PMCl:PMCh);DenFe372(PMCl:PMCh)];	% Altitude, Pure Fe density, 372nm Fe, 374nm Fe
	fpout = fopen(sprintf('/Users/dacostalindo/Desktop/Research/Results/McMurdo/ProcessedData//FeDensityData/FeDensity%04d%s%02d_%03d.dat',year,smonth(month,:),day,filenum),'w');
	fprintf(fpout,'NumChannel	Year	Month	Date	UTHour    StartTime374    StopTime374	StartTime372(UT)	StopTime372(UT)\n');
	fprintf(fpout,'%d	%d	%02d	%02d	%02d	%4.2f   %4.2f   %4.2f   %4.2f\n',2,year,month,day,timebin,time374(1),time374(end),time372(1),time372(end));
	fprintf(fpout,'Altitude(km)	PureFeDensity	374nmFe 372nmFe\n');
	fprintf(fpout,'%7.3f	%7.3f	%7.3f	%7.3f\n',FeDensityAll);
	fclose(fpout);

%% Plot section

    
%     subplot('Position',[0.15,0.65,0.75,0.3])
%     plot((1:BinNum)*Rawbinwid/1000,Normdpuresig374)
%     axis([20 150 -1 3])
%     hold on
%     plot([0 150],[0 0],'r')
%     xlabel('Altitude (km)')
%     ylabel('Normalized pure photon counts (Rayleigh removed)')
%     title(sprintf('374nm UT %4.2f-%4.2f [%02d/%02d/%04d]',time374(1),time374(end),month,day,year))
%     subplot('Position',[0.15,0.06,0.75,0.47])
   
    axes(handles.axes3)
    plot((1:BinNum)*Rawbinwid/1000,Normdpuresig374)
    axis([20 150 -1 3])
    hold on
    plot([0 150],[0 0],'r')
    xlabel('Altitude (km)')
    ylabel('Normalized pure photon counts (Rayleigh removed)')
    title(sprintf('374nm UT %4.2f-%4.2f [%02d/%02d/%04d]',time374(1),time374(end),month,day,year))
    hold off
    
%%    
    axes(handles.axes2)
    I =imread(sprintf('/Users/dacostalindo/Desktop/Research/Results/McMurdo/OverviewStatistics/OverviewPlots/JPG/McMurdoFe374nm%04d%s%02dGoodNew.jpg',year,smonth(month,:),day(end)) );
    imshow(I);     




%%

%       axes(handles.axes4)
%     I =imread(sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/OverviewStatistics/OverviewPlots/JPG/McMurdoFe372nm%04d%s%02dGoodNew.jpg',year,smonth(month,:),day(end)) );
%     imshow(I);     



    axes(handles.axes1)
   
    plot(([Altl:Alth-100]*Rawbinwid+AltCorrection)/1000,1e9*Beta_PMC([Altl:Alth-100]),'b',([Altl:Alth-100]*Rawbinwid+AltCorrection)/1000,1e9*(Beta_PMC([Altl:Alth-100])-dBeta([Altl:Alth-100])),'r')
     hold on
    xlabel('Altitude [km]')
    ylabel('Pure PMC signal Beta (x 10^-9 m^-1 sr^-1)')
    axis([75 110 -0.5 12 ])
    
    title(sprintf('374nm UT %4.2f-%4.2f [%02d/%02d/%04d]',time374(1),time374(end),month,day,year))
    if (FINDHfl==1 & FINDHfh==1)
	    title(sprintf('PMCyn= %d, BmaxLG0= %d, BFWHMLG0= %d, BLG1Err= %d, BLG2Err= %d, BmaxLG2Err= %d\n BmaxLGnoise=%d, ZBpk=%4.1fkm, FWHM=%3.1fkm, Bmax=%4.2f, ZHMl= %4.1fkm, ZHMh= %4.1fkm',PMCyn,BmaxLG0,BFWHMLG0,BLG1Err,BLG2Err,BmaxLG2Err,BmaxLGnoise,ZBpk/1000,FWHM/1000,1e9*Beta_max,ZHMlbin*Rawbinwid/1000,ZHMhbin*Rawbinwid/1000))
    else
	    title(sprintf('FINDHfl= %d, FINDHfh= %d, PMCyn= %d, BmaxLG0= %d, BmaxLG2Err= %d\n BmaxLGnoise=%d, ZBpk=%4.1fkm, Bmax=%4.2f',FINDHfl,FINDHfh,PMCyn,BmaxLG0,BmaxLG2Err,BmaxLGnoise,ZBpk/1000,1e9*Beta_max))
    end
    
    if (PMCyn==1)
	    text(97,1.5,sprintf('Zc=%4.2fkm\n Sigma=%3.2fkm\n Btotal=%5.2f\n Rmax=%5.2f\n PMCZcl=%4.2fkm\n PMCZch=%4.2fkm',Zc/1000,SigmaRMS/1000,1e6*Beta_Total,Rmax,PMCZcl*Rawbinwid/1000,PMCZch*Rawbinwid/1000))
    end
    pause(0.5)		% in order to show previous figures while processing the next hour
    hold off

    %Set the output values to the handles so that it can be output when output
    %function is called
    handles.output = [PMCyn;FeInfo];
%% TO PLOT IRON LAYERS
%     scrsz = get(0,'ScreenSize');
%     hndl = figure('position',[750 100 scrsz(3)*0.3 scrsz(4)*0.8])
%     orient portrait
%     set(gcf,'PaperPositionMode','auto')
%     subplot(3,1,1)
%     plot((Altl:Alth)*Rawbinwid/1000,DenFe374(Altl:Alth))
%     axis([70 120 -100 1000])
%     grid
%     xlabel('Altitude (km)')
%     ylabel('374nm Fe Density')
%     title(sprintf('372nm UT %4.2f-%4.2f [%02d/%02d/%04d]',time372(1),time372(end),month,day,year))
%     subplot(3,1,2)
%     plot((Altl:Alth)*Rawbinwid/1000,DenFe372(Altl:Alth))
%     grid
%     xlabel('Altitude [km]')
%     ylabel('372nm Fe Density')
%     axis([70 120 -100 2e4])
%     subplot(3,1,3)
%     plot((Altl:Alth)*Rawbinwid/1000,Density_Fe(Altl:Alth))
%     grid
%     xlabel('Altitude [km]')
%     ylabel('Pure Fe Density')
%     axis([70 120 -100 2e4])
%     pause(0.5)
%    saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/plots/formatFig/McMurdoFeDen%04d%s%02d_%03d.png',year,smonth(month,:),day(end),filenum));          
%    print('-djpeg',sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/plots/McMurdoFeDen%04d%s%02d_%03d',year,smonth(month,:),day(end),filenum)); 
%     close hndl


%%  Update handles structure and wait for OK button to be pressed
    guidata(hObject, handles);
    uiwait
	



% --- Outputs from this function are returned to the command line.
function varargout = interactive_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.output;



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Read upper and lower bound and display warnings in case the input is out
%of range or not numerical
upperbound = get(handles.upperbound, 'string');
lowerbound = get(handles.lowerbound, 'string')

if isempty(str2num(upperbound)) | ~isnumeric(str2num(upperbound)) | str2num(upperbound) > 100
    set(handles.upperbound,'string','insert upper bound');
    warndlg('Input must be numerical and less than 100 km');
    return
end


if isempty(str2num(lowerbound)) | ~isnumeric(str2num(lowerbound)) | str2num(lowerbound) < 50 | str2num(lowerbound) > str2num(upperbound)
    set(handles.lowerbound,'string','insert lower bound');
    warndlg('Input must be numerical and bigger than 50 km');
    return
end
axes(handles.axes1)

if handles.fig1 ~= 0 & handles.fig2 ~= 0
children = get(gca, 'children');
delete(findobj(children,'Tag','Line1'))  
delete(findobj(children,'Tag','Line2'))
    
end

hold on
h1 = plot([str2num(upperbound) str2num(upperbound)], [-0.5 12],'k','LineWidth',1)
set(h1,'Tag','Line1')
handles.fig1 = h1;
h2= plot([str2num(lowerbound) str2num(lowerbound)], [-0.5 12],'k','LineWidth',1)
set(h2,'Tag','Line2')
handles.fig2 = h2;
hold off

guidata(hObject,handles);





function lowerbound_Callback(hObject, eventdata, handles)
% hObject    handle to lowerbound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowerbound as text
%        str2double(get(hObject,'String')) returns contents of lowerbound as a double
%lb = str2double(get(hObject, 'String'));
%line([lb lb],get(gca,'YLim'))


% --- Executes during object creation, after setting all properties.
function lowerbound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowerbound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upperbound_Callback(hObject, eventdata, handles)
% hObject    handle to upperbound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upperbound as text
%        str2double(get(hObject,'String')) returns contents of upperbound as a double


% --- Executes during object creation, after setting all properties.
function upperbound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upperbound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)
%% saves PMC Data if PMC found, otherwise goes to next dataset
global PMCl PMCh Rawbinwid AltCorrection R dR Beta_PMC dBeta year smonth month day filenum timebin time374 time372 factorPMC factorNoise PMTsaturation_correction_yn
global Chopper_correction_yn Rmax ZRpk Beta_max ZBpk Beta_Total dBeta_Total SigmaRMS FWHM PMCZcl PMCZch Zpkbin Zc
%reset factor PMC and Factor Noise
factorPMC = 1.5;
factorNoise = 1.5;
set(handles.FactorNoise,'string','Factor Noise');
set(handles.FactorPMC,'string','Factor PMC');

if handles.PMCflag == 1

		AltRdRBdB = [((PMCl:PMCh)*Rawbinwid+AltCorrection)/1000;R(PMCl:PMCh);dR(PMCl:PMCh);1e9*Beta_PMC(PMCl:PMCh);1e9*dBeta(PMCl:PMCh)];	% Altitude, R, dR, Beta, dBeta
	    fpout = fopen(sprintf('/Users/dacostalindo/Desktop/Research/Results/McMurdo/ProcessedData/pmcData/PMCAWRB%04d%s%02d_%03d.dat',year,smonth(month,:),day,filenum),'w');
	    fprintf(fpout,'NumChannel	Year	Month	Date	UTHour	StartTime374(UT)	StopTime374(UT)    StartTime372(UT)	StopTime372(UT)\n');
        fprintf(fpout,'%d	%04d	%02d	%02d	%02d	%4.2f	%4.2f   %4.2f	%4.2f\n',2,year,month,day,timebin,time374(1),time374(end),time372(1),time372(end));
    	fprintf(fpout,'factorPMC   factorNoise	PMTsatCorrect ChopperCorrect\n');
	    fprintf(fpout,'%3.2f   %3.2f	%d	%d\n',factorPMC,factorNoise,PMTsaturation_correction_yn,Chopper_correction_yn);
	    fprintf(fpout,'PMCl(km)	PMCh(km)	PMCZcl(km) PMCZch(km)\n');
	    fprintf(fpout,'%7.3f	%7.3f	%7.3f	%7.3f\n',PMCl*Rawbinwid/1000,PMCh*Rawbinwid/1000,PMCZcl*Rawbinwid/1000,PMCZch*Rawbinwid/1000);
%	    fprintf(fpout,'Rmax	dRmax	ZRpk(km)	Bmax(10^-9 m-1sr-1)	dBmax	ZBpk(km)	TotalBeta(10^-6 sr-1)	dTotalBeta\n');
%	    fprintf(fpout,'%7.3f	%7.3f	%7.3f	%7.3f	%7.3f	%7.3f	%7.3f	%7.3f\n',Rmax,dRmax,ZRpk,Betamax,dBmax,ZBetapk,TotalBeta,dTotalBeta);	
%	    fprintf(fpout,'Zc(km)	dZc(km)	RMSwidth(km)	dRMSwidth(km)	FWHM_inner(km)	FWHM_outer(km)	FWHMl_outer	FWHMh_outer\n');
%	    fprintf(fpout,'%7.3f	%7.3f	%7.3f	%7.3f	%7.3f	%7.3f	%7.3f	%7.3f\n',Zc,dZc,SigmaRMS,dSigmaRMS,FWHM_inner,FWHM_outer,FWHMl_outer,FWHMh_outer);
	    fprintf(fpout,'Rmax	ZRpk(km)	Bmax(10^-9 m-1sr-1)	dBmax	ZBpk(km)	TotalBeta(10^-6 sr-1) deltaBtotal\n');
	    fprintf(fpout,'%7.3f	%7.3f	%7.3f	%7.3f	%7.3f	%7.3f	%7.3f\n',Rmax,ZRpk/1000,1e9*Beta_max,1e9*dBeta(Zpkbin),ZBpk/1000,1e6*Beta_Total,dBeta_Total);	
	    fprintf(fpout,'Zc(km)	RMSwidth(km)	FWHM(km)\n');
	    fprintf(fpout,'%7.3f	%7.3f	%7.3f\n',Zc/1000,SigmaRMS/1000,FWHM/1000);
	    fprintf(fpout,'Altitude(km)	BSRatio	BSRErr	Beta(10^-9 m-1sr-1)	BetaErr\n');
	    fprintf(fpout,'%7.3f	%7.3f	%7.3f	%7.3f	%7.3f\n',AltRdRBdB);
	    fclose(fpout);



end
guidata(hObject,handles);
uiresume


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in NoPMC.
function NoPMC_Callback(hObject, eventdata, handles)
% pressed it when the code find a PMC that in reality is non-existent. The
% function does not save the PMC data but goes straight to the next set
global factorPMC factorNoise
%reset factor PMC and Factor Noise
factorPMC = 1.5;
factorNoise = 1.5;
set(handles.FactorNoise,'string','Factor Noise');
set(handles.FactorPMC,'string','Factor PMC');
guidata(hObject,handles);
uiresume


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
%% Declare Global Variables
global PMCl PMCh Rawbinwid AltCorrection R dR Beta_PMC dBeta year smonth month day filenum timebin time374 time372 factorPMC factorNoise PMTsaturation_correction_yn
global Chopper_correction_yn Rmax ZRpk Beta_max ZBpk Beta_Total dBeta_Total SigmaRMS FWHM PMCZcl PMCZch Zpkbin Zc Altl Alth

%Get the boundary input from the user
upperbound = get(handles.upperbound, 'string');
lowerbound = get(handles.lowerbound, 'string');
%Get only numerical value within boundaries; Upper boundary has to be
%bigger than lower boundary
if isempty(str2num(upperbound)) | ~isnumeric(str2num(upperbound)) | str2num(upperbound) > 100
    set(handles.upperbound,'string','insert upper bound');
    warndlg('Upper bound must be numerical, less than 100 km, bigger than lower bound');
    return
end

if isempty(str2num(lowerbound)) | ~isnumeric(str2num(lowerbound)) | str2num(lowerbound) < 50 | str2num(lowerbound) > str2num(upperbound)
    set(handles.lowerbound,'string','insert lower bound');
    warndlg('Lower bound must be numerical, bigger than 50 km, lower than upper bound');
    return
end
%Transform upperbound and lowerbound to num
upperbound = str2num(upperbound);
lowerbound = str2num(lowerbound);

PMCZl = round(lowerbound*1e3/Rawbinwid);
PMCZh = round(upperbound*1e3/Rawbinwid);





% PMCZh = round(86e3/Rawbinwid);          % for 2 UT on 2011/Feb/12
%% Find the maximum peak in 75-90km range and its corresponding index
[Beta_max,ZpkIndex] = max(Beta_PMC(PMCZl:PMCZh));

Zpkbin = PMCZl-1+ZpkIndex;                  % bin number for the peak Beta
ZBpk = Rawbinwid*Zpkbin + AltCorrection;    % convert the index to altitude
%% !!!!!!!!!!!!!
dBeta(Zpkbin)                               % display for monitoring purpose
%% Judge whether Bmax is larger than 2*dBeta
if ((Beta_max - 2*abs(dBeta(Zpkbin))) > 0)
	BmaxLG2Err = 1;	% Bmax is larger than 2 times of error
else
	BmaxLG2Err = 0;
end

if (Beta_max > 0)
	BmaxLG0 = 1;
	FINDHfl = 0;
	m = 1;
	while (FINDHfl == 0 & m < (Zpkbin-PMCZl))		% to find the FWHM point at lower end
		if (Beta_PMC(Zpkbin-m) <= Beta_max/2 & Beta_PMC(Zpkbin-m+1) > Beta_max/2)
			ZHMlbin = Zpkbin-m;
			FINDHfl = 1;
		end
		m = m+1;
	end
	FINDHfh = 0;
	n = 1;	
	while (FINDHfh == 0 & n< (PMCZh-Zpkbin))		% to find the FWHM point at higher end
		if (Beta_PMC(Zpkbin+n-1) > Beta_max/2 & Beta_PMC(Zpkbin+n) <= Beta_max/2)
			ZHMhbin = Zpkbin+n;
			FINDHfh = 1;
		end
		n = n+1;
	end
	if (FINDHfl == 1 & FINDHfh == 1)
		FWHM = (ZHMhbin-ZHMlbin+1)*Rawbinwid;
		if (min(Beta_PMC(ZHMlbin:ZHMhbin)) > 0)
			BFWHMLG0 = 1;
			[Diff1,L1index] = min(Beta_PMC(ZHMlbin:ZHMhbin) - abs(dBeta(ZHMlbin:ZHMhbin)));
			if (Diff1>0)
				BLG1Err = 1;	% All Beta within FWHM range are larger than 1 time of error
			else
				BLG1Err = 0;
            end
           
            %factorPMC = 1.5;  % change to 1.5 times of error (for most days)
%             factorPMC = 1.2;  % change to 1.2 times of error (for a few hours in 12/27/03, 2 hours in 12/12/03, 1 hour on 2010/12/21)
%             factorPMC = 1.05;  % change to 1.1 times of error (for a few hours in 12/10/04, 1 hour on 2010/12/26)
			%[Diff2,L2index]=min(Beta_PMC(ZHMlbin:ZHMhbin)-2*dBeta(ZHMlbin:ZHMhbin));   
			%[Diff2,L2index]=min(Beta_PMC(ZHMlbin:ZHMhbin)-1.5*dBeta(ZHMlbin:ZHMhbin));  % change to 1.5 times of error (for most days)
			[Diff2,L2index] = min(Beta_PMC(ZHMlbin:ZHMhbin)-factorPMC*dBeta(ZHMlbin:ZHMhbin));  % change to 1.2 times of error (for a few hours in 12/27/03, 2 hours in 12/12/03)
			%[Diff2,L2index]=min(Beta_PMC(ZHMlbin:ZHMhbin)-1.1*dBeta(ZHMlbin:ZHMhbin)); 
			if (Diff2 > 0)	% This standard is too high!!!
				BLG2Err = 1;	% All Beta within FWHM range are larger than 2 times of error
            else
				BLG2Err = 0;
            end
			
		else
			BFWHMLG0 = 0;		% all bins within FWHM range should be larger than zero
			sprintf('there are negative bins within FWHM,although Bmax > 0')
			[Diff1,L1index] = min(Beta_PMC(ZHMlbin:ZHMhbin)-abs(dBeta(ZHMlbin:ZHMhbin)));
			if (Diff1 > 0)
				BLG1Err = 1;	% All Beta within FWHM range are larger than 1 time of error
			else
				BLG1Err = 0;
			end
            %factorPMC = 1.5;  % change to 1.5 times of error (for most
            %days)  ALSO UNCOMMENT IN CASE DO NOT NEED ADVANCED COMMANDS
            %ANYMORE
			[Diff2,L2index] = min(Beta_PMC(ZHMlbin:ZHMhbin)-factorPMC*dBeta(ZHMlbin:ZHMhbin));  % change to 1.2 times of error (for a few hours in 12/27/03, 2 hours in 12/12/03)
			if (Diff2 > 0)	% This standard is too high!!!
				BLG2Err = 1;	% All Beta within FWHM range are larger than 2 times of error
            else
				BLG2Err = 0;
            end
		end
	else
		sprintf('Cannot find FWHM points: FINDHfl, FINDHfh')
		BFWHMLG0 = NaN;
		BLG1Err = NaN;
        factorPMC = NaN;
	end
else
	BmaxLG0 = 0;		% Bmax should be larger than zero
	sprintf('All are negative photon counts (Bmax < 0), check background subtraction')
    PMCyn = NaN;
    FeInfo = NaN;
    return
end

%% Compute the standard deviation of noise (take range of 95-105km for now)
% change the range for background noise computation to 95-105 km
noiseZl = round(95e3/Rawbinwid);	% for most data
noiseZh = round(105e3/Rawbinwid);	% for most data
%noiseZl=round(90e3/Rawbinwid);	% for 29 Dec 2002 data
%noiseZh=round(100e3/Rawbinwid);	% for 29 Dec 2002 data
%noiseZl=round(90e3/Rawbinwid);	% for 19 Jan 2002 59UT data
%noiseZh=round(100e3/Rawbinwid);	% for 19 Jan 2002 59UT data
noisemean = mean(Beta_PMC(noiseZl:noiseZh));
stdnoise = std(Beta_PMC(noiseZl:noiseZh))		% divided by (N-1)
%factorNoise=2;
%factorNoise = 1.5;
% factorNoise=1.2;        % for 2UT on 2011/02/12
if (Beta_max > factorNoise*stdnoise)		% define as 2 times of std deviation
	BmaxLGnoise=1;
else
	BmaxLGnoise=0;
end

%% Judge whether the peak is qualified to be a PMC
%if (BmaxLG0==1 & BFWHMLG0==1 & BLG2Err==1 & BmaxLGnoise==1)
%if (BmaxLG0==1 & BFWHMLG0==1 & BLG1Err==1 & BmaxLG2Err==1 & BmaxLGnoise==1)
%if (BmaxLG0==1 & BFWHMLG0==1 & BLG1Err==1 & BLG2Err==1 & BmaxLG2Err==1)    % remove stdnoise condition for 17UT on 12/27/2003, 15UT on 12/12/03
if (BmaxLG0==1 & BFWHMLG0==1 & BLG1Err==1 & BLG2Err==1 & BmaxLG2Err==1 & BmaxLGnoise==1)    % BLG2Err is actually for 1.5 times of error
	PMCyn = 1;
    handles.PMCflag = 1;
else
	PMCyn = 0;
end






%% Compute the PMC parameters after qualifying a peak as PMC
if (PMCyn==1)
	FIND0l=0;
	m=1;
	while (FIND0l==0)		% to find the zero point at lower end
		if (Beta_PMC(Zpkbin-m)<=0 & Beta_PMC(Zpkbin-m+1)>0)
			PMCZcl=Zpkbin-m;		% PMC layer start bin
			FIND0l=1
		end
		m=m+1;
	end
	FIND0h=0;
	n=1;	
	while (FIND0h==0)		% to find the zero point at higher end
		if (Beta_PMC(Zpkbin+n-1)>0 & Beta_PMC(Zpkbin+n)<=0)
			PMCZch=Zpkbin+n;		% PMC layer stop bin
			FIND0h=1
		end
		n=n+1;
    end
    
	if (FIND0l==1 & FIND0h==1)
		if ((PMCZch-Zpkbin)*Rawbinwid>3.5e3)	% if half side of PMC layer is larger than 3.5 km
			rangeHfPMC = round(3e3/Rawbinwid);			% 3.0km corresponding bin number		
			[Bminh,minhIndex] = min(Beta_PMC(Zpkbin:Zpkbin+rangeHfPMC));
			if (abs(Bminh) < Beta_max/5)		% if the minimumm less than 1/10th of Bmax, take this minimum as zero point for PMC stop
				PMCZch = Zpkbin-1+minhIndex;
				refinedPMCZch = PMCZch*Rawbinwid
            else
                PMCZch=Zpkbin+rangeHfPMC;
			end
        end
        
		if ((Zpkbin-PMCZcl)*Rawbinwid>4e3)
			rangeHfPMC=round(3.5e3/Rawbinwid)			% 3.5km corresponding bin number		
			[Bminl,minlIndex]=min(Beta_PMC(Zpkbin-rangeHfPMC:Zpkbin))
			if (abs(Bminl)<Beta_max/5)		% if the minimumm less than 1/10th of Bmax, take this minimum as zero point for PMC start
				PMCZcl=Zpkbin-rangeHfPMC-1+minlIndex;
            else
                PMCZcl=Zpkbin-rangeHfPMC;
			end
		end
	else
		sprintf('Cannot find zero point for PMC layer start or stop')
	end
	% after find zero points for both start and stop, compute PMC parameters
	if (FIND0l==1 & FIND0h==1)
		% centroid altitude (km)
		Zc = sum(Beta_PMC(PMCZcl:PMCZch).*(PMCZcl:PMCZch)*Rawbinwid)/sum(Beta_PMC(PMCZcl:PMCZch));
		% rms width (km)
		SigmaRMS = sqrt(sum(((PMCZcl:PMCZch)*Rawbinwid-Zc).^2.*Beta_PMC(PMCZcl:PMCZch))/sum(Beta_PMC(PMCZcl:PMCZch)));
		% total backscatter coefficient
		Beta_Total = sum(Beta_PMC(PMCZcl:PMCZch).*Rawbinwid);
        dBeta_Total = sum(dBeta(PMCZcl:PMCZch))*1e6;
		% maximum backscatter ratio of PMC layer
		[Rmax,RmaxIndex] = max(R(PMCZcl:PMCZch));
		ZRpkbin = PMCZcl-1+RmaxIndex;
		ZRpk = ZRpkbin*Rawbinwid + AltCorrection;
        Zc = Zc + AltCorrection;
	end
end
%Plot on axis 1
    axes(handles.axes1)
    plot(([Altl:Alth-100]*Rawbinwid+AltCorrection)/1000,1e9*Beta_PMC([Altl:Alth-100]),'b',([Altl:Alth-100]*Rawbinwid+AltCorrection)/1000,1e9*(Beta_PMC([Altl:Alth-100])-dBeta([Altl:Alth-100])),'r')
     hold on
    xlabel('Altitude [km]')
    ylabel('Pure PMC signal Beta (x 10^-9 m^-1 sr^-1)')
    axis([75 110 -0.5 12 ])
    
    title(sprintf('374nm UT %4.2f-%4.2f [%02d/%02d/%04d]',time374(1),time374(end),month,day,year))
    if (FINDHfl==1 & FINDHfh==1)
	    title(sprintf('PMCyn= %d, BmaxLG0= %d, BFWHMLG0= %d, BLG1Err= %d, BLG2Err= %d, BmaxLG2Err= %d\n BmaxLGnoise=%d, ZBpk=%4.1fkm, FWHM=%3.1fkm, Bmax=%4.2f, ZHMl= %4.1fkm, ZHMh= %4.1fkm',PMCyn,BmaxLG0,BFWHMLG0,BLG1Err,BLG2Err,BmaxLG2Err,BmaxLGnoise,ZBpk/1000,FWHM/1000,1e9*Beta_max,ZHMlbin*Rawbinwid/1000,ZHMhbin*Rawbinwid/1000))
    else
	    title(sprintf('FINDHfl= %d, FINDHfh= %d, PMCyn= %d, BmaxLG0= %d, BmaxLG2Err= %d\n BmaxLGnoise=%d, ZBpk=%4.1fkm, Bmax=%4.2f',FINDHfl,FINDHfh,PMCyn,BmaxLG0,BmaxLG2Err,BmaxLGnoise,ZBpk/1000,1e9*Beta_max))
    end
    
    if (PMCyn==1)
	    text(97,1.5,sprintf('Zc=%4.2fkm\n Sigma=%3.2fkm\n Btotal=%5.2f\n Rmax=%5.2f\n PMCZcl=%4.2fkm\n PMCZch=%4.2fkm',Zc/1000,SigmaRMS/1000,1e6*Beta_Total,Rmax,PMCZcl*Rawbinwid/1000,PMCZch*Rawbinwid/1000))
    end
    pause(0.5)		% in order to show previous figures while processing the next hour
    hold off
    guidata(hObject, handles);



function FactorPMC_Callback(hObject, eventdata, handles)
global factorPMC
%% get factor PMC and use this value for it
 factorPMC = str2num(get(hObject, 'String'))
 if isempty(factorPMC) | ~isnumeric(factorPMC)
    factorPMC = 1.5; % reset PMC factor
    set(hObject,'string','Factor PMC');
    warndlg('Input must be numerical ');
    
    return
 end



% --- Executes during object creation, after setting all properties.
function FactorPMC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FactorPMC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FactorNoise_Callback(hObject, eventdata, handles)
    global factorNoise
    %% get factor noise and use this value for it
     factorNoise = str2double(get(hObject, 'String'));
     if isempty(factorNoise) | ~isnumeric(factorNoise)
        factorNoise = 1.5; % reset PMC factor
        set(hObject,'string','Factor Noise');
        warndlg('Input must be numerical ');
        return
     end



% --- Executes during object creation, after setting all properties.
function FactorNoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FactorNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global year smonth day month
if handles.contour == 0
    figure;
    I =imread(sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/OverviewStatistics/OverviewPlots/JPG/McMurdoFe372nm%04d%s%02dGoodNew.jpg',year,smonth(month,:),day(end)) );
    imshow(I);  
    handles.contour = 1;
end
guidata(hObject, handles);

      








