clear all
close all
clc

%change 
%NOTE: this main is supposed to be run after main.m as it corrects the PMC
%findings by either adding or eliminating data sets in PMC Data and in plots.  

%PURPOSE: 
%1- screen by looking over the plots where

% %Season 2011-2012
% % %  year = [2012]; month = [12]; day = [08];
% year = [2011 2011 2012 2012];  month = [11 12 01 02]; day = [
%                                                              17 18 00 28 29 30 00 00 00 00 00 00 00 00 00 00; ...
%                                                              11 14 15 19 20 21 23 24 25 26 27 29 30 31 00 00; ... 
%                                                              05 00 10 11 12 15 16 00 00 29 00 31 00 00 00 00;  ...
%                                                              01 03 04 07 08 09 15 17 20 21 23 24 00 00 00 00; ... 
%                                                              ];
% INDEX = 1;

% % % % %Season 2012-2013
year = [2012 2012 2013 2013];  month = [11 12 01 02]; day = [
                                                             00 05 00 00 09 00 16 17 00 19 22 23 00 00 00 00; ...
                                                             08 09 17 27 28 29 31 00 00 00 00 00 00 00 00 00; ...
                                                             01 02 00 04 00 10 11 21 22 23 24 25 29 00 00 00; ...
                                                             09 20 24 26 00 00 00 00 00 00 00 00 00 00 00 00; 
                                                             ];    
INDEX = 2;
                                                       
% % %Season 2013-2014
% year = [2013 2013 2014 2014];  month = [11 12 01 02]; day = [
%                                                              12 14 00 18 19 00 26 29 00 00 00 00 00 00 00 00; ...
%                                                              01 09 00 11 12 00 23 24 00 27 28 29 00 31 00 00; ...
%                                                              10 14 19 20 27 28 00 00 00 00 00 00 00 00 00 00; ...
%                                                              05 06 00 11 12 21 22 26 27 00 00 00 00 00 00 00; ...
%                                                              ];                                                         
%                                                         
% INDEX = 3;                                                                                  
% % %                                                          
% % % %Season 2014-2015
% year = [2014 2014 2015 2015];  month = [11 12 01 02]; day = [
%                                                              01 02 03 10 11 00 16 21 00 00 00 00 00 00 00 00; ...
%                                                              00 02 04 05 08 09 10 11 00 13 17 18 19 20 28 31; ...
%                                                              01 11 21 22 23 26 27 28 29 30 00 00 00 00 00 00; ... 
%                                                              11 12 14 16 17 18 26 00 00 00 00 00 00 00 00 00
%                                                              ];                                                         
% INDEX = 4;                                                       




% % %Season 2015-2016
% year = [2015 2015 2016 2016];  month = [11 12 01 02]; day = [
%                                                                08 09 00 26 27 30 00 00 00 00 00 00 00 00 00 00;
%                                                                01 02 03 17 22 24 26 27 00 00 00 00 00 00 00 00; ...
%                                                                05 06 13 14 15 16 22 23 24 25 00 00 00 00 00 00; ...
%                                                                04 05 06 07 14 16 17 18 27 00 00 00 00 00 00 00; ...
%                                                                
%                                                              ];                                                         
% INDEX = 5;                                                       

% % % % %Season 2016-2017
% year = [2016 2016 2017 2017];  month = [11 12 01 02]; day = [
%                                                                00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00;
%                                                                03 04 12 13 14 00 00 00 00 00 00 00 00 00 00 00; ...
%                                                                02 03 07 10 11 16 18 19 20 24 25 26 28 29 30 31;
%                                                                01 02 03 05 06 11 12 14 15 16 18 20 00 00 00 00; ...
%                                                              
%                                                              ];                                                         
% INDEX = 6;                                                       

%                                                                   
 %TEST
year = 2013; month = 01; day = 04;
% %   
% year = [2014 2014 2015 2015];  month = [11 12 01 02]; day = [
%                                                              00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00; ...
%                                                              00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00; ...
%                                                              01 11 21 22 23 26 27 28 29 30 00 00 00 00 00 00; ... 
%                                                              00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
%                                                              ];                                                         
%                                                       
% 

%make a vector containing initials for each month                                                                
smonth = ['JA';'FB';'MR';'AR';'MY';'JN';'JL';'AG';'SP';'OT';'NV';'DC'];                                                                        
%take vector sizes for reading in data purposes
[yr,yc] = size(year);
[mr,mc] = size(month);
[dr,dc] = size(day);
%Open text file for 374nm statistics
% name_374 =  sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/OverviewStatistics/StatisticData/%s374.txt',smonth(month(1,:),:));
% name_372 =  sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/OverviewStatistics/StatisticData/%s372.txt',smonth(month(1,:),:));
% fpout_374 = fopen(name_374,'w');
% fpout_372 = fopen(name_372,'w');
% fprintf(fpout_374,'Year    Month    Day    Mean    Max    Min    STD \n');
% fprintf(fpout_372,'Year    Month    Day    Mean    Max    Min    STD \n');
%set count
count = 1;

if yc == mc && mc == dr 
    for i = 1:yc 
            for k = 1:dc
                if day(i,k)~= 0 && month(i)~= 0
                    
                    
                            
                              
                              
                               %Make filelist
%                                 path = ['/Users/asta10/Desktop/Research/SkyObservations/' num2str(year(i)) smonth(month(i),:) num2str(day(i,k),'%02i')];
%                                 makefilelist(path,1,1);
                    
                              %Eliminate Bad Data (Setting ScreenRawData = 1) 
                              % -eliminating data based on SNR and creating
                              % Goodfilelist
                          %    ScreenRawData = 1;
                           %   [SNR_std, SNR_mean,SNR_min,SNR_max] =  screenSNR(year(i),month(i),day(i,k),ScreenRawData);
                            %  close all
                   
                              %Overview
                             % overviewFunction(year(i),month(i),day(i,k))
                             % close all
                              
                              %Process PMC Data
                              mainfunctionprocess_retouch(year(i),month(i),day(i,k))  
                              close all
% % % % %                     
% % % % %                    
% % % % %                     %Statistics
% % % % %                     ScreenRawData = 0;
% % % % %                     [SNR_std, SNR_mean,SNR_min,SNR_max] =  statistic(year(i),month(i,j),day(i,k),ScreenRawData);
% % % % %                     close all
% % % % %                     DataInfo374(count,:) = [year(i) month(i,j) day(i,k) SNR_mean(1,1) ,SNR_max(1,1),SNR_min(1,1),SNR_std(1,1)];
% % % % %                     DataInfo372(count,:) = [year(i) month(i,j) day(i,k) SNR_mean(1,2), SNR_max(1,2),SNR_min(1,2),SNR_std(1,2)];
% % % % %                     count = count+ 1;

                             
                              
                end
             end
 
        end
     


else 
disp('enter same number of columns for year, month, day')
return
end

