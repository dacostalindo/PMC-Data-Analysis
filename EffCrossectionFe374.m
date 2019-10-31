function [sigma_eff] = EffCrossectionFe374(freqlaser,temp,wind)

%% Claim global variables for constants and parameters
global pi c h Me Qe E0 kB AMU NA mFeMean mFeIsotope AbundanceFe DeltaE_diff 
global Aki_Fe372 gk_Fe372 gi_Fe372 RB_Fe372 Lambda_Center_Fe372 fosc_Fe372 IsotopeShift_Fe372
global Aki_Fe374 gk_Fe374 gi_Fe374 RB_Fe374 Lambda_Center_Fe374 fosc_Fe374 
global PAL374Eportion PAL372Eportion PAL374RMSwidth PAL372RMSwidth 
global PAL374Detune PAL372Detune PAL374WL PAL372WL dEoverKB

% actual laser frequency
freq = freqlaser;
% compute convoluted rms width
SigmaD = sqrt(kB*temp./mFeMean/(Lambda_Center_Fe374^2));
Sigmae1 = sqrt(SigmaD.^2 + PAL374RMSwidth(1)^2);
Sigmae2 = sqrt(SigmaD.^2 + PAL374RMSwidth(2)^2);
% compute effective cross-section
comon1 = fosc_Fe374*Qe^2/(sqrt(2*pi)*4*E0*Me*c)./Sigmae1;
comon2 = fosc_Fe374*Qe^2/(sqrt(2*pi)*4*E0*Me*c)./Sigmae2;
sigma_eff1 = comon1.*exp(-(freq - wind/Lambda_Center_Fe374).^2/2./(Sigmae1.^2));
                % portion determined by the narrow peak
sigma_eff2 = comon2.*exp(-(freq - wind/Lambda_Center_Fe374).^2/2./(Sigmae2.^2));
                % portion determined by the pedestal
sigma_eff = PAL374Eportion(1)*sigma_eff1 + PAL374Eportion(2)*sigma_eff2;
                % total effective cross-section contributed by the narrow
                % peak and the pedestal, 
                % PAL374Eportion(1) + PAL374Eportion(2) = 1 by definition
                
% Here, positive wind means the Fe atoms move along the same direction as
% the laser photon propagation, i.e., moving away from the laser source