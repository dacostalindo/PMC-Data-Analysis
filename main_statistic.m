clearvars -except Mean_R_max Mean_B_max Mean_B_tot Mean_Zc Mean_Sigma smdBTOT X smdZc smdRMS PMC_MATRIX MATRIX M D PMC_Hours Tot_hours
close all
clc
%PURPOSE: count the amount of hours of observation in a whole month and
%during the same months through several years
%NOTE: since matrices need to have same number of columns and rows, fill in
%the days and months and then fill with 0s the rest.
% %%%Example:
% % %NOVEMBER OF EACH YEAR
% year = [ 2011 2012 2013 2014 2015]; month = [ 11 11 11 11 11]; day =  [     17 18 00 28 29 30 00 00 00 00 00 00 00; ... 
%                                                                             00 05 00 00 09 00 16 17 00 19 22 23 00; ...
%                                                                             12 14 00 18 19 00 26 29 00 00 00 00 00 ; ...
%                                                                             01 02 03 10 11 00 16 21 00 00 00 00 00; ...
%                                                                             08 09 00 26 27 30 00 00 00 00 00 00 00 ];
% % %    
% % % % % % DECEMBER OF EACH YEAR
% year = [2011 2012 2013 2014 2015 2016]; month = [ 12 12 12 12 12 12]; day =  [      11 14 15 19 20 21 23 24 25 26 27 29 30 31 00 00; ... 
%                                                                                     08 09 17 27 28 29 31 00 00 00 00 00 00 00 00 00; ...
%                                                                                     01 09 00 11 12 22 23 24 25 27 28 29 00 31 00 00; ...
%                                                                                     02 04 05 08 09 10 11 12 13 17 18 19 20 28 00 31; ...
%                                                                                     01 02 03 17 22 24 26 27 28 00 00 00 00 00 00 00; ...
%                                                                                     03 04 12 13 14 00 00 00 00 00 00 00 00 00 00 00];
% 
% % %       
% % %JANUARY OF EACH YEAR
% year = [2012 2013 2014 2015 2016 2017]; month = [ 01 01 01 01 01 01]; day =  [05 00 10 11 12 15 16 00 00 29 00 31 00 00 00 00;  ... 
%                                                                                     01 02 00 04 00 10 11 21 22 23 24 25 29 00 00 00; ...
%                                                                                     10 14 19 20 27 28 00 00 00 00 00 00 00 00 00 00; ...
%                                                                                     01 11 21 22 23 26 27 28 29 30 00 00 00 00 00 00; ...
%                                                                                     05 06 13 14 15 16 22 23 24 25 00 00 00 00 00 00; ...
%                                                                                     02 03 07 10 11 16 18 19 20 24 25 26 28 29 30 31 ];    

% 
% % % % % %    %FEBRUARY OF EACH YEAR
%    year = [2012 2013 2014 2015]; month = [ 02 02 02 02]; day =  [                   01 03 04 07 08 09 15 17 20 21 23 24; ... 
%                                                                                     09 20 24 26 00 00 00 00 00 00 00 00; ...
%                                                                                     05 06 00 11 12 21 22 26 27 00 00 00; ...
%                                                                                     11 12 14 16 17 18 26 00 00 00 00 00];
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SEASON DATA
% % % % %Season 2011-2012
% year = [2011 2011 2012 2012];  month = [11 12 01 02]; day = [
%                                                              17 18 00 28 29 30 00 00 00 00 00 00 00 00 00 00; ...
%                                                              11 14 15 19 20 21 23 24 25 26 27 29 30 31 00 00; ... 
%                                                              05 00 10 11 12 15 16 00 00 29 00 31 00 00 00 00;  ...
%                                                              01 03 04 07 08 09 15 17 20 21 23 24 00 00 00 00; ... 
%                                                              ];
% INDEX = 1;

% %%Filtered
year = [2011 2011 2012 2012];  month = [11 12 01 02]; day = [
                                                             00 00 00 28 29 30 00 00 00 00 00 00 00 00 00 00; ...
                                                             11 14 15 19 20 21 23 24 25 26 27 29 30 31 00 00; ... 
                                                             05 00 10 11 12 15 16 00 00 29 00 31 00 00 00 00;  ...
                                                             01 03 04 07 08 09 00 00 00 00 00 00 00 00 00 00; ... 
                                                             ];                                          
INDEX = 1;

% % % % % % % % % % % %Season 2012-2013
% year = [2012 2012 2013 2013];  month = [11 12 01 02]; day = [
%                                                              00 05 00 00 09 00 16 17 00 19 22 23 00 00 00 00; ...
%                                                              08 09 17 27 28 29 31 00 00 00 00 00 00 00 00 00; ...
%                                                              01 02 00 04 00 10 11 21 22 23 24 25 29 00 00 00; ...
%                                                              09 20 24 26 00 00 00 00 00 00 00 00 00 00 00 00; 
%                                                              ];    
% INDEX = 2;

% %%Filtered
% 
% year = [2012 2012 2013 2013];  month = [11 12 01 02]; day = [
%                                                              00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00; ...
%                                                              08 09 17 27 28 29 31 00 00 00 00 00 00 00 00 00; ...
%                                                              01 02 00 04 00 10 11 21 22 23 24 25 29 00 00 00; ...
%                                                              09 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00; 
%                                                              ];    
% INDEX = 2;

% % % % % % % %                                                        
% % % % % % %Season 2013-2014
% year = [2013 2013 2014 2014];  month = [11 12 01 02]; day = [
%                                                              12 14 00 18 19 00 26 29 00 00 00 00 00 00 00 00; ...
%                                                              01 09 00 11 12 00 23 24 00 27 28 29 00 31 00 00; ...
%                                                              10 14 19 20 27 28 00 00 00 00 00 00 00 00 00 00; ...
%                                                              05 06 00 11 12 21 22 26 27 00 00 00 00 00 00 00; ...
%                                                              ];                                                         
%                                                         
% INDEX = 3; 
% 
% 
% %% Filtered
% year = [2013 2013 2014 2014];  month = [11 12 01 02]; day = [
%                                                              00 00 00 00 00 00 26 29 00 00 00 00 00 00 00 00; ...
%                                                              01 09 00 11 12 00 23 24 00 27 28 29 00 31 00 00; ...
%                                                              10 14 19 20 27 28 00 00 00 00 00 00 00 00 00 00; ...
%                                                              05 06 00 11 12 00 00 00 00 00 00 00 00 00 00 00; ...
%                                                              ];                                                         
%                                                         
% INDEX = 3;  
% %                                                          
% % % % % % %Season 2014-2015
% year = [2014 2014 2015 2015];  month = [11 12 01 02]; day = [
%                                                              01 02 03 10 11 00 16 21 00 00 00 00 00 00 00 00; ...
%                                                              00 02 04 05 08 09 10 11 00 13 17 18 19 20 28 31; ...
%                                                              01 11 21 22 23 26 27 28 29 30 00 00 00 00 00 00; ... 
%                                                              11 12 14 16 17 18 26 00 00 00 00 00 00 00 00 00
%                                                              ];                                                         
% INDEX = 4;                                                       

% %%FILTERED
% year = [2014 2014 2015 2015];  month = [11 12 01 02]; day = [
%                                                              00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00; ...
%                                                              00 02 04 05 08 09 10 11 00 13 17 18 19 20 28 31; ...
%                                                              01 11 21 22 23 26 27 28 29 30 00 00 00 00 00 00; ... 
%                                                              11 12 14 16 00 00 00 00 00 00 00 00 00 00 00 00
%                                                              ];                                                         
% INDEX = 4;                                                       

% 
% % 
% % % % % % %Season 2015-2016
% year = [2015 2015 2016 2016];  month = [11 12 01 02]; day = [
%                                                                08 09 00 26 27 30 00 00 00 00 00 00 00 00 00 00;
%                                                                01 02 03 17 22 24 26 27 00 00 00 00 00 00 00 00; ...
%                                                                05 06 13 14 15 16 22 23 24 25 00 00 00 00 00 00; ...
%                                                                04 05 06 07 14 16 17 18 27 00 00 00 00 00 00 00; ...
%                                                                
%                                                              ];                                                         
% INDEX = 5;                                                       
% %%filtered
% year = [2015 2015 2016 2016];  month = [11 12 01 02]; day = [
%                                                                00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00;
%                                                                00 00 00 00 00 24 26 27 00 00 00 00 00 00 00 00; ...
%                                                                05 06 13 14 15 16 22 23 24 25 00 00 00 00 00 00; ...
%                                                                04 05 06 07 14 16 17 18 00 00 00 00 00 00 00 00; ...
%                                                                
%                                                              ];                                                         
%  INDEX = 5;  
% % % % % %Season 2016-2017
% year = [2016 2016 2017 2017];  month = [11 12 01 02]; day = [
%                                                                00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00;
%                                                                03 04 12 13 14 00 00 00 00 00 00 00 00 00 00 00; ...
%                                                                02 03 07 10 11 16 18 19 20 24 25 26 28 29 30 31;
%                                                                01 02 03 05 06 11 12 14 15 16 18 20 00 00 00 00; ...
%                                                              
%                                                              ];                                                         
% INDEX = 6;                                                       

% %Filtered
% year = [2016 2016 2017 2017];  month = [11 12 01 02]; day = [
%                                                                00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00;
%                                                                03 04 12 13 14 00 00 00 00 00 00 00 00 00 00 00; ...
%                                                                02 03 07 10 11 16 18 19 20 24 25 26 28 29 30 31;
%                                                                01 02 03 05 06 11 12 14 00 00 00 00 00 00 00 00; ...
%                                                              
%                                                              ];                                                         
% INDEX = 6;   
% % % % % % % 
% 
%  



%PARAMETERS TO CHANGE
Interannual_Analysis = false;
plot_interannual_variation = false; % Set to true only after having run all the seasons with Interannual_Analysis == true

Seasonal_Analysis = false;
plot_seasonal_variation = false;  % Set to true only after having run all the seasons with Seasonal_Analysis == true

Diurnal_Analysis = true;
plot_diurnal_analysis = true; %Set to true only after having run all the seasons with Diurnal_Analysis == true    

plot_overall_Histogram = 1;

%make a vector containing initials for each month                                                                
smonth = ['JA';'FB';'MR';'AR';'MY';'JN';'JL';'AG';'SP';'OT';'NV';'DC'];                                                                        
%take vector sizes for reading in data purposes
[yr,yc] = size(year);
[mr,mc] = size(month);
[dr,dc] = size(day);
%%%%%%%%%%%%%%%%%%%%%%%%DIURNAL ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Diurnal_Analysis == true

%Create MATRIX 
name =  sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/SeasonalAnalysis/Season_%i_%i_best_daysdiurnal.txt',year(1),year(4));
name_xls =  sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/SeasonalAnalysis/Season_%i_%i_best_daysdiurnal.csv',year(1),year(4));
fpout = fopen(name,'w');
%set count
count1 = zeros(length(year),1);
count2 = zeros(length(year),1);

if yc == mc && mc == dr

%Give path to folder
path = '/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/pmcData/';
path2 = '/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/FeDensityData/';

%initialize count
Cnt = 1;
Cnt2 = 1;
    for i = 1:yc 
            
       
            for k = 1:dc
                if day(i,k)~= 0 && month(i)~= 0
                    %Get number of files in both folders
                    folders = dir([path '*'  num2str(year(i)) smonth(month(i),:) num2str(day(i,k),'%02i') '*']);
                    folders2 = dir([path2 '*'  num2str(year(i)) smonth(month(i),:) num2str(day(i,k),'%02i') '*']);
                    if length(folders)== 0 
                      fprintf('%d/%d/%d does not have any data \n',year(i),month(i), day(i,k))
                    end
                    if length(folders2)== 0 
                      fprintf('%d/%d/%d does not have any data \n',year(i),month(i), day(i,k))
                    end
                    [rows,cols] = size(folders);
                    [rows2,cols2] = size(folders2);
                    
                    %Calculate hours of detection of PMC
                    for ii = 1:rows
                        
                    filename = folders(ii).name;
                    [year(i) month(i) day(i,k) UT_Hour zuluday(Cnt) Rmax(Cnt) Bmax(Cnt) dBmax(Cnt) Btotal(Cnt) dBtotal(Cnt) ZcPMC(Cnt) RMSwidth(Cnt) FWHM(Cnt) ZBpk(Cnt) ZRpk(Cnt) factorPMC(Cnt) factorNoise(Cnt) duration(Cnt)]= retrievePMCData([path filename]);
                    %create Matrix to print to txt file
                     Matrix_PMC(Cnt,:) = [year(i) month(i) day(i,k) UT_Hour zuluday(Cnt) Rmax(Cnt) Bmax(Cnt) dBmax(Cnt) Btotal(Cnt) dBtotal(Cnt) ZcPMC(Cnt) RMSwidth(Cnt) FWHM(Cnt) ZBpk(Cnt) ZRpk(Cnt) factorPMC(Cnt) factorNoise(Cnt)];
                   
                     Cnt = Cnt+1;
                     
                    end
                    
                    
                    %Calculate hours of observation
                    for ii = 1:rows2
                        
                    filename2 = folders2(ii).name;
                    [Duration(Cnt2),ZuluD(Cnt2),UT_H(Cnt2)]= retrieveFeDensity([path2 filename2]);
                    %create Matrix to print to txt file
                     Matrix(Cnt2,:) = [ZuluD(Cnt2) Duration(Cnt2) UT_H(Cnt2)];
                   
                     Cnt2 = Cnt2+1;
                     
                    end
                    
                    
                    
                    
                                      
                end
            end
 
       
     

    end
    


else 
disp('enter same number of columns for year, month, day')
return


end

siz = size(Matrix_PMC)

for c = 1:siz(1)
    
    if Matrix_PMC(i,4) > 24
       Matrix_PMC(i,4) = mod(Matrix_PMC(i,4), 24)    
    end
    
%     if UT_Hour(i) > 24
%        UT_Hour(i) = mod(UT_Hour(i), 24)    
%     end
%     
       
end

% 
% %Write data to txt file
fprintf(fpout,'Year Month Day UT_Hour zuluday Rmax Bmax dBmax Btotal dBtotal ZcPMC RMSwidth FWHM ZBpk ZRpk factorPMC factorNoise \n');
title_string = {'Year','Month' 'Day', 'UT_Hour', 'zuluday', 'Rmax', 'Bmax', 'dBmax', 'Btotal','dBtotal' ,  'ZcPMC', 'RMSwidth', 'FWHM', 'ZBpk', 'ZRpk', 'factorPMC', 'factorNoise'};
dlmwrite(name,Matrix_PMC,'-append','delimiter',' ');
%write to excel file
 fid = fopen(name_xls, 'w') ;
 fprintf(fid, '%s,', title_string{1,1:end-1}) ;
 fprintf(fid, '%s\n', title_string{1,end}) ;
 fclose(fid) ;
 dlmwrite(name_xls, Matrix_PMC(1:end,:), '-append') ;


%Create seasonal variations versus day number relative to summer solstice
%create summer solsitce table
if year(1) ==2010
sumSol = 21+9.38/24; %[in days (and hours)]
elseif year(1) == 2011
sumSol = 22+5.30/24; %[in days (and hours)]
elseif year(1) == 2012
sumSol = 21+11.12/24; %[in days (and hours)]
elseif year(1) == 2013
sumSol = 21+17.11/24; %[in days (and hours)]   
elseif year(1) == 2014
sumSol = 21+23.03/24; %[in days (and hours)]
elseif year(1) == 2015
sumSol = 22+4.48/24; %[in days (and hours)]
elseif year(1) == 2016
sumSol = 21+10.44/24; %[in days (and hours)]
end
    
zulumonth = [31;28;31;30;31;30;31;31;30;31;30;31];		% day numbers in each month through whole year
zulu_sumSol_day = -(sum(zulumonth(1:12-1))+sumSol-365);	% convert to zulu day 

% plot(Matrix_PMC(:,5)+zulu_sumSol_day,Matrix_PMC(:,6),'x')


% %     %Occurrence probability 
    Matrix_PMC(:,5) = Matrix_PMC(:,5)+zulu_sumSol_day;
    Matrix(:,1) = Matrix(:,1)+zulu_sumSol_day;
    
    
    PMC_MATRIX{INDEX} = Matrix_PMC;
    MATRIX{INDEX} = Matrix; 
    
    %% Plot Observation hours and occurrence probability
    for i = 0:23
        ii = find(Matrix_PMC(:,4) == i);
        pmc_found(i+1) = sum(length(ii));
        j = find(Matrix(:,3)== i);
        observationHours(i+1) = sum(length(j));
        occurrence(i+1) = pmc_found(i+1)/observationHours(i+1);
    end
    
    %get number of hours
    PMC_Hours(INDEX) = Cnt-1;
    Tot_hours(INDEX) = Cnt2-1;
    
    
    scrsz = get(0,'ScreenSize');
            figure('Position',[200 70 scrsz(3)*0.3 scrsz(4)*0.9])
            orient portrait; set(gcf,'PaperPositionMode','auto')
            
            %Subplot 1 - Peak Backscatter Ratio
            subplot(2,1,1);
    histogram(Matrix_PMC(:,4),24);
    ax = gca;
    ysize = ax.YLim 
    ylabel('PMC Observed per Hour')
     xlabel('UT\_Hour')
     t = title(['Season ' num2str(year(2)) '-' num2str(year(3))]);
     set(t, 'FontSize', 18);
       text(18,ysize(2)-3,sprintf('PMC Hours= %i\n Tot Hours= %i\n',PMC_Hours(INDEX),Tot_hours(INDEX)))
     
          subplot(2,1,2);
    plot(0:23,occurrence,'o');
    ylim([0 1])
    ylabel('Occurrence Probability')
     xlabel('UT\_Hour') 
     saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/DiurnalAnalysis/PercentageOccurrenceseason%04d_%04d.jpg',year(2),year(3)));          
   
     
    %% Plot mean diurnal parameters
    if plot_diurnal_analysis == true
       
        
        
        I = 6; %This corresponds to the number of seasons
        
% % %         %%Real one
%         if add_season_first_season == true
%             season1_PMC_Matrix = xlsread('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/SeasonalAnalysis/McMurdoPMCstat1stSeasontxt.xlsx');
%             I = 1;
%             PMC_MATRIX{I} = season1_PMC_Matrix;
%         end   
%       
% %         %%experimenting stupid
%         if add_season_first_season == true
%             clear PMC_MATRIX
%             I = 1;  
%             PMC_MATRIX{I} = Matrix_PMC;
%         end   
    
%          if plot_only_season1 == true
%             I = 1;
%             clear PMC_MATRIX
%             PMC_MATRIX{I} = season1_PMC_Matrix;
%          end   
%          
%         if Test_Best_days == true
%             clear PMC_MATRIX
%             I = 1;
%             PMC_MATRIX{I} = xlsread('/Users/asta10/Desktop/Research/TESTDAY/BESTDAYS3'); 
%         end
        
        B_tot_sum = zeros(24,1);
        B_max_sum = zeros(24,1);
        Zc_sum = zeros(24,1);
        R_max_sum = zeros(24,1);
        n_terms = zeros(24,1);
         h = [0:23];
         %Do averages over the 6 seasons of data
       for k = 1:length(h)  
           for i = 1:I
              
              loc = find(PMC_MATRIX{i}(:,4)== h(k)); %%  & PMC_MATRIX{i}(:,7)> 2
              n_terms(k) = length(loc) +n_terms(k);
              B_max_sum(k) = sum(PMC_MATRIX{i}(loc,7))+B_max_sum(k);
              B_tot_sum(k) = sum(PMC_MATRIX{i}(loc,9))+B_tot_sum(k);
              Zc_sum(k) =   sum(PMC_MATRIX{i}(loc,11))+Zc_sum(k);
              R_max_sum(k) = sum(PMC_MATRIX{i}(loc,6))+ R_max_sum(k);
           end 
           
       end 
% %        %% avoid NaN loop
%        for var = 1:24
%            if n_terms(var) == 0
%               n_terms(var)=1; 
%            end
%                
%            
%        end


      
            
            %Subplot 1 - Peak Backscatter Ratio
    figure;
    bar(1:24,n_terms);
    ax =gca;
    ysize = ax.YLim
    text(18,ysize(2)-3,sprintf('PMC Hours= %i\n Tot Hours= %i\n',sum(PMC_Hours),sum(Tot_hours)))
    ylabel('Total Observation Hour')
     xlabel('UT\_Hour')
     t = title('Overall Season Stats');
     set(t, 'FontSize', 18);
       
       
       
       
       %calculate the mean values
       B_tot_mean = B_tot_sum./n_terms;
       Zc_mean = Zc_sum./n_terms;
       B_max_mean = B_max_sum./n_terms;
       R_max_mean = R_max_sum./n_terms;
       
       
        %% Do a fitting to data
        timeday = 0:23;
        timeday = timeday';
        ftype=fittype({'1','sin(2*pi*x*1/12)','cos(2*pi*x*1/12)','sin(2*pi*x*1/24)','cos(2*pi*x*1/24)'},'coeff',{'a0','a1','a2','b1','b2'});

        
    
    %%%plot mean values
     %Plot histograms for season variations
            scrsz = get(0,'ScreenSize');
            figure('Position',[200 70 scrsz(3)*0.3 scrsz(4)*0.9])
            orient portrait; set(gcf,'PaperPositionMode','auto')
          
           
             % Subplot 1- in peak backscatter ratio vs time 
             [f1,goodness,output]= fit(timeday,R_max_mean ,ftype);
             yfit =f1.a0+f1.a1*sin(2*pi*timeday*1/12)+f1.a2*cos(2*pi*timeday*1/12)+f1.b1*sin(2*pi*timeday*1/24)+f1.b2*cos(2*pi*timeday*1/24);
             A0_R = f1.a0
             A12_R = sqrt((f1.a1)^2+(f1.a2)^2)
             A24_R = sqrt((f1.b1)^2+(f1.b2)^2)
             
             UT12_RTemp = acos(f1.a2/A12_R)/pi*180;
             UT12_R_checkaTemp = asin(f1.a1/A12_R)/pi*180;
             
             UT24_RTemp = acos(f1.b2/A24_R)/pi*180;
             UT24_R_checkTemp = asin(f1.b1/A24_R)/pi*180;
             
             UT12_R = acos(f1.a2/A12_R)*12/(2*pi)
             UT12_R_checka = asin(f1.a1/A12_R)*12/(2*pi)
             
             UT24_R = asin(f1.b1/A24_R)*24/(2*pi)
             UT24_R_check = asin(f1.b1/A24_R)*24/(2*pi)
             
             
             corrcoef(yfit,R_max_mean')
             yfit1 =A0_R+A12_R*cos(2*pi/12*(timeday-UT12_R))+A24_R*cos(2*pi/24*(timeday-UT24_R));
            
             subplot(4,1,1); plot(0:1:23,R_max_mean,'s') 
             hold on
            % plot(timeday,yfit,'-k', 'LineWidth',2);
            % plot(timeday,yfit1,'--r','LineWidth',2)
             
             hold off
             hold off
             ylabel('R_{max}')
             xlabel('UT\_Hour')
             t = title('Diurnal Variations of Mean Parameters');
             set(t, 'FontSize', 17);
              x0=[f1.a0,sqrt(f1.a1^2+f1.a2^2),UT12_R,sqrt(f1.b1^2+f1.b2^2),UT24_R];
[beta1,resnorm,resid,exitflag,output,lambda,J]=lsqcurvefit(@HarmonicFitting,x0,timeday,R_max_mean);
disp(beta1);
ci1=nlparci(beta1,resid,'jacobian',J)
error1=(ci1(:,2)-ci1(:,1))/2;
disp(error1)

             
            
            
             % Subplot 2 - Volume Backscatter  Coefficient vs time 
             [f1,goodness,output]= fit(timeday,B_max_mean ,ftype);
             yfit =f1.a0+f1.a1*sin(2*pi*timeday*1/12)+f1.a2*cos(2*pi*timeday*1/12)+f1.b1*sin(2*pi*timeday*1/24)+f1.b2*cos(2*pi*timeday*1/24);
             A0_R = f1.a0
             A12_R = sqrt((f1.a1)^2+(f1.a2)^2)
             A24_R = sqrt((f1.b1)^2+(f1.b2)^2)
             
             UT12_RTemp = acos(f1.a2/A12_R)/pi*180;
             UT12_R_checkaTemp = asin(f1.a1/A12_R)/pi*180;
             
             UT24_RTemp = acos(f1.b2/A24_R)/pi*180;
             UT24_R_checkTemp = asin(f1.b1/A24_R)/pi*180;
             
             UT12_R = acos(f1.a2/A12_R)*12/(2*pi)
             UT12_R_checka = asin(f1.a1/A12_R)*12/(2*pi)
             
             UT24_R = asin(f1.b1/A24_R)*24/(2*pi)
             UT24_R_check = asin(f1.b1/A24_R)*24/(2*pi)
             
             
             corrcoef(yfit,B_max_mean')
             yfit1 =A0_R+A12_R*cos(2*pi/12*(timeday-UT12_R))+A24_R*cos(2*pi/24*(timeday-UT24_R));
            
             subplot(4,1,2); plot(0:1:23,B_max_mean,'s') 
             hold on
             %plot(timeday,yfit,'-k', 'LineWidth',2);
            % plot(timeday,yfit1,'--r','LineWidth',2)
             
             hold off
             ylabel('\beta_{max} (10^{-9} m^{-1}sr^{-1})')
             xlabel('UT\_Hour')
              x0=[f1.a0,sqrt(f1.a1^2+f1.a2^2),UT12_R,sqrt(f1.b1^2+f1.b2^2),UT24_R];
[beta1,resnorm,resid,exitflag,output,lambda,J]=lsqcurvefit(@HarmonicFitting,x0,timeday,B_max_mean);
disp(beta1);
ci1=nlparci(beta1,resid,'jacobian',J)
error1=(ci1(:,2)-ci1(:,1))/2;
disp(error1)

     
             % Subplot 3 - Total Volume Backscatter  Coefficient vs time  
             [f1,goodness,output]= fit(timeday,B_tot_mean ,ftype);
             yfit =f1.a0+f1.a1*sin(2*pi*timeday*1/12)+f1.a2*cos(2*pi*timeday*1/12)+f1.b1*sin(2*pi*timeday*1/24)+f1.b2*cos(2*pi*timeday*1/24);
             A0_R = f1.a0
             A12_R = sqrt((f1.a1)^2+(f1.a2)^2)
             A24_R = sqrt((f1.b1)^2+(f1.b2)^2)
             
             UT12_RTemp = acos(f1.a2/A12_R)/pi*180;
             UT12_R_checkaTemp = asin(f1.a1/A12_R)/pi*180;
             
             UT24_RTemp = acos(f1.b2/A24_R)/pi*180;
             UT24_R_checkTemp = asin(f1.b1/A24_R)/pi*180;
             
             UT12_R = acos(f1.a2/A12_R)*12/(2*pi)
             UT12_R_checka = asin(f1.a1/A12_R)*12/(2*pi)
             
             UT24_R = asin(f1.b1/A24_R)*24/(2*pi)
             UT24_R_check = asin(f1.b1/A24_R)*24/(2*pi)
             
             
             corrcoef(yfit,B_tot_mean')
             yfit1 =A0_R+A12_R*cos(2*pi/12*(timeday-UT12_R))+A24_R*cos(2*pi/24*(timeday-UT24_R));
            
             subplot(4,1,3); plot(0:1:23,B_tot_mean,'s') 
             hold on
            % plot(timeday,yfit,'-k', 'LineWidth',2);
            % plot(timeday,yfit1,'--r','LineWidth',2)
             
             hold off
             ylabel('\beta_{Total} (10^{-6} sr^{-1})')
             xlabel('UT\_Hour')
           
             x0=[f1.a0,sqrt(f1.a1^2+f1.a2^2),UT12_R,sqrt(f1.b1^2+f1.b2^2),UT24_R];
[beta1,resnorm,resid,exitflag,output,lambda,J]=lsqcurvefit(@HarmonicFitting,x0,timeday,B_tot_mean);
disp(beta1);
ci1=nlparci(beta1,resid,'jacobian',J)
error1=(ci1(:,2)-ci1(:,1))/2;
disp(error1)

             
             
             % Subplot 4 - Centroid Altitude vs time
             [f1,goodness,output]= fit(timeday,Zc_mean ,ftype);
             yfit =f1.a0+f1.a1*sin(2*pi*timeday*1/12)+f1.a2*cos(2*pi*timeday*1/12)+f1.b1*sin(2*pi*timeday*1/24)+f1.b2*cos(2*pi*timeday*1/24);
             A0_R = f1.a0
             A12_R = sqrt((f1.a1)^2+(f1.a2)^2)
             A24_R = sqrt((f1.b1)^2+(f1.b2)^2)
             
             UT12_RTemp = acos(f1.a2/A12_R)/pi*180;
             UT12_R_checkaTemp = asin(f1.a1/A12_R)/pi*180;
             
             UT24_RTemp = acos(f1.b2/A24_R)/pi*180;
             UT24_R_checkTemp = asin(f1.b1/A24_R)/pi*180;
             
             UT12_R = acos(f1.a2/A12_R)*12/(2*pi)
             UT12_R_checka = asin(f1.a1/A12_R)*12/(2*pi)
             
             UT24_R = acos(f1.b2/A24_R)*24/(2*pi)
             UT24_R_check = asin(f1.b1/A24_R)*24/(2*pi)
             
             
             corrcoef(yfit,Zc_mean')
             yfit1 =A0_R+A12_R*cos(2*pi/12*(timeday-UT12_R_checka))+A24_R*cos(2*pi/24*(timeday-UT24_R));
            
             subplot(4,1,4); plot(0:1:23,Zc_mean,'s') 
             hold on
            % plot(timeday,yfit,'-k', 'LineWidth',2);
            % plot(timeday,yfit1,'--r','LineWidth',2)
             
             hold off
             ylabel('Z_c (km)')
             xlabel('UT\_Hour')
             
             
x0=[f1.a0,sqrt(f1.a1^2+f1.a2^2),UT12_R_checka,sqrt(f1.b1^2+f1.b2^2),UT24_R];
[beta1,resnorm,resid,exitflag,output,lambda,J]=lsqcurvefit(@HarmonicFitting,x0,timeday,Zc_mean);
disp(beta1);
ci1=nlparci(beta1,resid,'jacobian',J)
error1=(ci1(:,2)-ci1(:,1))/2;

             
 saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/DiurnalAnalysis/mean_diurnal_parameters.jpg'));          
   
    
    
   
end 



end
% %%%%%%%%%%%%%%%%%%%%%%SEASONAL ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Seasonal_Analysis == true

    if add_season_first_season == false
%Create MATRIX 
name =  sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/SeasonalAnalysis/Season_%i_%i.txt',year(1),year(4));
name_xls =  sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/SeasonalAnalysis/Season_%i_%i.csv',year(1),year(4));
fpout = fopen(name,'w');
%set count
count1 = zeros(length(year),1);
count2 = zeros(length(year),1);

if yc == mc && mc == dr

%Give path to folder
path = '/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/pmcData/';
path2 = '/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/FeDensityData/';

%initialize count
Cnt = 1;
Cnt2 = 1;
    for i = 1:yc 
            
       
            for k = 1:dc
                if day(i,k)~= 0 && month(i)~= 0
                    %Get number of files in both folders
                    folders = dir([path '*'  num2str(year(i)) smonth(month(i),:) num2str(day(i,k),'%02i') '*']);
                    folders2 = dir([path2 '*'  num2str(year(i)) smonth(month(i),:) num2str(day(i,k),'%02i') '*']);
                    if length(folders)== 0 
                      fprintf('%d/%d/%d does not have any data \n',year(i),month(i), day(i,k))
                    end
                    if length(folders2)== 0 
                      fprintf('%d/%d/%d does not have any data \n',year(i),month(i), day(i,k))
                    end
                    [rows,cols] = size(folders);
                    [rows2,cols2] = size(folders2);
                    
                    %Calculate hours of detection of PMC
                    for ii = 1:rows
                        
                    filename = folders(ii).name;
                    [year(i) month(i) day(i,k) UT_Hour zuluday(Cnt) Rmax(Cnt) Bmax(Cnt) dBmax(Cnt) Btotal(Cnt) dBtotal(Cnt) ZcPMC(Cnt) RMSwidth(Cnt) FWHM(Cnt) ZBpk(Cnt) ZRpk(Cnt) factorPMC(Cnt) factorNoise(Cnt) duration(Cnt)]= retrievePMCData([path filename]);
                    %create Matrix to print to txt file
                     Matrix_PMC(Cnt,:) = [year(i) month(i) day(i,k) UT_Hour zuluday(Cnt) Rmax(Cnt) Bmax(Cnt) dBmax(Cnt) Btotal(Cnt) dBtotal(Cnt) ZcPMC(Cnt) RMSwidth(Cnt) FWHM(Cnt) ZBpk(Cnt) ZRpk(Cnt) factorPMC(Cnt) factorNoise(Cnt)];
                   
                     Cnt = Cnt+1;
                     
                    end
                    
                    
                    %Calculate hours of observation
                    for ii = 1:rows2
                        
                    filename2 = folders2(ii).name;
                    [Duration(Cnt2),ZuluD(Cnt2)]= retrieveFeDensity([path2 filename2]);
                    %create Matrix to print to txt file
                     Matrix(Cnt2,:) = [ZuluD(Cnt2) Duration(Cnt2)];
                   
                     Cnt2 = Cnt2+1;
                     
                    end
                    
                    
                    
                    
                                      
                end
            end
 
       
     

    end
    


else 
disp('enter same number of columns for year, month, day')
return


end


%Write data to txt file
fprintf(fpout,'Year Month Day UT_Hour zuluday Rmax Bmax dBmax Btotal dBtotal ZcPMC RMSwidth FWHM ZBpk ZRpk factorPMC factorNoise \n');
title_string = {'Year','Month' 'Day', 'UT_Hour', 'zuluday', 'Rmax', 'Bmax', 'dBmax', 'Btotal','dBtotal' ,  'ZcPMC', 'RMSwidth', 'FWHM', 'ZBpk', 'ZRpk', 'factorPMC', 'factorNoise'};
dlmwrite(name,Matrix_PMC,'-append','delimiter',' ');
%write to excel file
 fid = fopen(name_xls, 'w') ;
 fprintf(fid, '%s,', title_string{1,1:end-1}) ;
 fprintf(fid, '%s\n', title_string{1,end}) ;
 fclose(fid) ;
 dlmwrite(name_xls, Matrix_PMC(1:end,:), '-append') ;

    end 
    
   %add the data from season 1
    if add_season_first_season == true
    INDEX = 1;
   Matrix_PMC = xlsread('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/SeasonalAnalysis/McMurdoPMCstat1stSeasontxt.xlsx');
   UT_Hour = Matrix_PMC(:,4);
   year(1) = 2010;
   year(2) = 2010;
   year(3) = 2011;
   year(4) = 2011;
    end    
    
    
%Create seasonal variations versus day number relative to summer solstice
%create summer solsitce table
if year(1) == 2010
sumSol = 21+23.38/24; %[in days (and hours)]    
elseif year(1) == 2011
sumSol = 22+5.30/24; %[in days (and hours)]
elseif year(1) == 2012
sumSol = 21+11.12/24; %[in days (and hours)]
elseif year(1) == 2013
sumSol = 21+17.11/24; %[in days (and hours)]   
elseif year(1) == 2014
sumSol = 21+23.03/24; %[in days (and hours)]
elseif year(1) == 2015
sumSol = 22+4.48/24; %[in days (and hours)]
elseif year(1) == 2016
sumSol = 21+10.44/24; %[in days (and hours)]
end
    
zulumonth = [31;28;31;30;31;30;31;31;30;31;30;31];		% day numbers in each month through whole year
zulu_sumSol_day = -(sum(zulumonth(1:12-1))+sumSol+UT_Hour/24-365);	% convert to zulu day 



% plot(Matrix_PMC(:,5)+zulu_sumSol_day,Matrix_PMC(:,6),'x')

% Plot graphs for season variations
    scrsz = get(0,'ScreenSize');
            figure('Position',[200 70 scrsz(3)*0.3 scrsz(4)*0.9])
            orient portrait; set(gcf,'PaperPositionMode','auto')
            
            %Subplot 1 - Peak Backscatter Ratio
            subplot(5,1,1); plot(Matrix_PMC(:,5)+zulu_sumSol_day,Matrix_PMC(:,6),'x')
            ylabel('R_{max}')
            t = title([num2str(year(2)) '-' num2str(year(3)) ' Season '])
            set(t, 'FontSize', 18);
%             Mean_R_max(INDEX) = mean(DataInfo374_PMC(:,4));
%             STD_R_max = std(DataInfo374_PMC(:,4));
%             Standard_error_R_max = STD_R_max/sqrt(length(DataInfo374_PMC(:,4)));
           % ylim([0 200])
        
            
            
             %Subplot 2 - Peak Volume Backscatter  Coefficient \Beta_{max}
            subplot(5,1,2); plot(Matrix_PMC(:,5)+zulu_sumSol_day,Matrix_PMC(:,7),'x')
            ylabel('\beta_{max} (10^{-6} sr^{-1})')
            
             %Subplot 3 - Peak Volume Backscatter  Coefficient \Beta_{tot}
            
            subplot(5,1,3); plot(Matrix_PMC(:,5)+zulu_sumSol_day,Matrix_PMC(:,9),'x')
            ylabel(' \beta_{Total} (10^{-6} sr^{-1})')

             %Subplot 4 - Centroid Altitude Z_c (km)
            
            subplot(5,1,4); plot(Matrix_PMC(:,5)+zulu_sumSol_day,Matrix_PMC(:,11),'x')
            ylabel('Z_c (km)')
        
             %Subplot 5 - RMS Width \sigma_{rms} (km)
            
            subplot(5,1,5); plot(Matrix_PMC(:,5)+zulu_sumSol_day,Matrix_PMC(:,12),'x')
            ylabel('\sigma_{rms} (km)')
           
 scrsz = get(0,'ScreenSize');
            figure('Position',[200 70 scrsz(3)*0.3 scrsz(4)*0.9])
            orient portrait; set(gcf,'PaperPositionMode','auto')
            
            %Subplot 1 - Peak Backscatter Ratio
           
            ;
%             Mean_R_max(INDEX) = mean(DataInfo374_PMC(:,4));
%             STD_R_max = std(DataInfo374_PMC(:,4));
%             Standard_error_R_max = STD_R_max/sqrt(length(DataInfo374_PMC(:,4)));
           % ylim([0 200])
        
            
         
             %Subplot 3 - Peak Volume Backscatter  Coefficient \Beta_{tot}
            
            subplot(3,1,1); plot(Matrix_PMC(:,5)+zulu_sumSol_day,Matrix_PMC(:,9),'x')
            ylabel(' \beta_{Total} (10^{-6} sr^{-1})')
            t = title([num2str(year(2)) '-' num2str(year(3)) ' Season '])
            set(t, 'FontSize', 18)

             %Subplot 4 - Centroid Altitude Z_c (km)
            
            subplot(3,1,2); plot(Matrix_PMC(:,5)+zulu_sumSol_day,Matrix_PMC(:,11),'x')
            ylabel('Z_c (km)')
        
             %Subplot 5 - RMS Width \sigma_{rms} (km)
            
            subplot(3,1,3); plot(Matrix_PMC(:,5)+zulu_sumSol_day,Matrix_PMC(:,12),'x')
            ylabel('\sigma_{rms} (km)')
    saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/SeasonalAnalysis/FIG/season%04d_%04d.fig',year(2),year(3)));          
    print('-djpeg',sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/SeasonalAnalysis/season%04d_%04d',year(2),year(3))); 
    
  if add_season_first_season == false  
%     %Occurrence probability 
    Matrix_PMC(:,5) = Matrix_PMC(:,5)+zulu_sumSol_day;
    Matrix(:,1) = Matrix(:,1)+zulu_sumSol_day;
    [startSmooth_pos] = Matrix(find(Matrix(:,1)>=0),1);
    [startSmooth_neg] = Matrix(find(Matrix(:,1)< 0),1);
    
    range_pos = startSmooth_pos(1):11:startSmooth_pos(end);
    range_neg = startSmooth_neg(end):-11:startSmooth_neg(1);
    
    if range_pos(end) ~= startSmooth_pos(end)
        range_pos(end+1) = startSmooth_pos(end);
    end    
    
    if range_neg(end) ~= startSmooth_neg(1)
        range_neg(end+1) = startSmooth_neg(1);
    end 
    
    for i = 1:(length(range_pos)-1)
        indexes_pos{i} = find(Matrix(:,1) >= range_pos(i) & Matrix(:,1)<= range_pos(i+1));
        indexes_pos_PMC{i} = find(Matrix_PMC(:,5) >= range_pos(i) & Matrix_PMC(:,5)<= range_pos(i+1));
        occurence_pos(i) = sum(duration(indexes_pos_PMC{i}))/sum(Matrix(indexes_pos{i},2));
        medium_pos(i) = (range_pos(i)+ range_pos(i+1))/2;
        error_pos(i) =  range_pos(i+1)-range_pos(i);
    end
    
     for i = 1:(length(range_neg)-1)
        indexes_neg{i} = find(Matrix(:,1) <= range_neg(i) & Matrix(:,1)>= range_neg(i+1));
        indexes_neg_PMC{i} = find(Matrix_PMC(:,5) <= range_neg(i) & Matrix_PMC(:,5)>= range_neg(i+1));
        occurence_neg(i) = sum(duration(indexes_neg_PMC{i}))/sum(Matrix(indexes_neg{i},2));
        medium_neg(i) = (range_neg(i)+ range_neg(i+1))/2;
        error_neg(i) =  abs(range_neg(i+1)-range_neg(i));
     end
     
     occurrence = [fliplr(occurence_neg) occurence_pos];
     for i = 1:length(occurrence)
        if  occurrence(i) < 0.05 |  occurrence(i) > 0.96
            occurrence(i) = NaN;
        end
     end
     medium = [fliplr(medium_neg) medium_pos];
     error = [fliplr(error_neg)/2 error_pos/2];
     time0 = [fliplr(range_neg) range_pos];
     figure;
     herrorbar(medium,occurrence,error,'o')
     ylabel('Occurrence Probability')
     xlabel('Day Number Relative to Summer Solstice')
     ylim([0 1])
     t = title(['Season ' num2str(year(2)) '-' num2str(year(3))]);
     set(t, 'FontSize', 18);
      saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/SeasonalAnalysis/PercentageOccurrenceseason%04d_%04d.jpg',year(2),year(3)));          
    %print('-djpeg',sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/barCharts/season%04d_%04d',year(2),year(3))); 

    % plot(medium,occurrence,'o')

%      for i = 1:(length(time0)-1)
%      indices{i} =  find(Matrix(:,1) >= time0(i) & Matrix(:,1)<= time0(i+1));
%      indices_PMC{i} = find(Matrix_PMC(:,5) >= time0(i) & Matrix_PMC(:,5)<= time0(i+1));
%      occurrence(i) = sum(duration(indices_PMC{i}))/sum(Matrix(indices{i},2));
%      medium_point(i) = (time0(i)+time0(i+1))/2;
%      end
%    
  end
     
     
     %% smooth data using 11 days window and Hamming Smooth
FW = 11; %[days]
[smdBTOT{INDEX}] = HammingSmth(FW,Matrix_PMC(:,5),Matrix_PMC(:,5),Matrix_PMC(:,9));
[smdZc{INDEX}] = HammingSmth(FW,Matrix_PMC(:,5),Matrix_PMC(:,5),Matrix_PMC(:,11));
[smdRMS{INDEX}] = HammingSmth(FW,Matrix_PMC(:,5),Matrix_PMC(:,5),Matrix_PMC(:,12));
X{INDEX} = Matrix_PMC(:,5);

%% Plot smoothed PMC parameters vs day number
if plot_seasonal_variation == true
%scrsz = get(0,'ScreenSize');
%            figure('Position',[200 70 scrsz(3)*0.3 scrsz(4)*0.9])
 %           orient portrait; set(gcf,'PaperPositionMode','auto')
            
            %Plot 1 - Smoothed Total Backscatter Coefficient
            figure;
            hold on
            plot(X{1},smdBTOT{1},'-o')
            plot(X{2},smdBTOT{2},'-o')
            plot(X{3},smdBTOT{3},'-o')
            plot(X{4},smdBTOT{4},'-o')
            plot(X{5},smdBTOT{5},'-o')
            plot(X{6},smdBTOT{6},'-o')
            plot(X{7},smdBTOT{7},'-o')
            ylabel('Smoothed Total Backscatter Coefficient \beta_{Total} (10^{-6} sr^{-1})')
            t = title('Smoothed Total Backscatter Coefficient  vs day number')
            set(t, 'FontSize', 18);
           legend('2010-2011 Season','2011-2012 Season', '2012-2013 Season','2013-2014 Season', '2014-2015 Season','2015-2016 Season','2016-2017 Season')
            hold off
           saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/MeanStatisticPlots/SmoothedBackscatterParameters.jpg'));          
    
            %Plot 2 - Smoothed Centroid Altitude
            figure;
            hold on
            plot(X{1},smdZc{1},'-o')
            plot(X{2},smdZc{2},'-o')
            plot(X{3},smdZc{3},'-o')
            plot(X{4},smdZc{4},'-o')
            plot(X{5},smdZc{5},'-o')
            plot(X{6},smdZc{6},'-o')
            plot(X{7},smdZc{7},'-o')
            ylabel('Smoothed Centroid Altitude Z_c (km)')
            t = title('Smoothed Centroid Altitude vs day number')
            set(t, 'FontSize', 18);
           legend('2010-2011 Season','2011-2012 Season', '2012-2013 Season','2013-2014 Season', '2014-2015 Season','2015-2016 Season','2016-2017 Season')
            hold off
            saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/MeanStatisticPlots/SmoothedCentroidParameters.jpg'));          
    
            %Plot 3 - Smoothed RMS Width (km)
            figure;
            hold on
            plot(X{1},smdRMS{1},'-o')
            plot(X{2},smdRMS{2},'-o')
            plot(X{3},smdRMS{3},'-o')
            plot(X{4},smdRMS{4},'-o')
            plot(X{5},smdRMS{5},'-o')
            plot(X{6},smdRMS{6},'-o')
            plot(X{7},smdRMS{7},'-o')
            ylabel('Smoothed RMS Width \sigma_{rms}(km)')
            t = title('Smoothed RMS Width vs day number')
            set(t, 'FontSize', 18);
           legend('2010-2011 Season','2011-2012 Season', '2012-2013 Season','2013-2014 Season', '2014-2015 Season','2015-2016 Season','2016-2017 Season')
            hold off
            saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/MeanStatisticPlots/SmoothedRMSParameters.jpg'));          
    
             for i = 1:length(smdZc)
                 
             corrcoef(smdZc{i},smdBTOT{i},'rows','complete')
             end



end







end

%%%%%%%%%%%%%%%%%%%%%%%%%%INTERANNUAL ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Interannual_Analysis == true
%set count
count1 = zeros(length(year),1);
count2 = zeros(length(year),1);

if yc == mc && mc == dr
    
if add_season_first_season == false
%Give path to folder
path = '/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/FeDensityData/';
path_PMC = '/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/pmcData/';

%initialize count
Cnt = 1;
Cnt_PMC =1;
    for i = 1:yc 
            monthlycnt1(i) = 0;
            monthlycnt2(i) = 0;
            monthlycnt1_PMC(i) = 0;
            monthlycnt2_PMC(i) = 0;
            
       
            for k = 1:dc
                if day(i,k)~= 0 && month(i)~= 0
                    %Get number of files in both folders
                    folders = dir([path '*'  num2str(year(i)) smonth(month(i),:) num2str(day(i,k),'%02i') '*']);
                    folders_PMC = dir([path_PMC '*'  num2str(year(i)) smonth(month(i),:) num2str(day(i,k),'%02i') '*']);
                    if length(folders)== 0 
                      fprintf('%d/%d/%d does not have any data \n',year(i),month(i), day(i,k))
                    end
                    [rows,cols] = size(folders);
                    [rows_PMC,cols_PMC] = size(folders_PMC);
                    %Calculate hours of observation
                    for ii = 1:rows
                        
                    filename = folders(ii).name;
                    [count1(Cnt),count2(Cnt),startHour1(Cnt),endHour1(Cnt),startHour2(Cnt),endHour2(Cnt)]= countHours([path filename]);
                    %create Matrix to print to txt file
                     DataInfo374(Cnt,:) = [year(i) month(i) day(i,k) count1(Cnt) startHour1(Cnt) endHour1(Cnt)];
                     DataInfo372(Cnt,:) = [year(i) month(i) day(i,k) count2(Cnt) startHour2(Cnt) endHour2(Cnt)];
                     
                     monthlycnt1(i) = monthlycnt1(i)+count1(Cnt);
                     monthlycnt2(i) = monthlycnt2(i)+count2(Cnt);
                     Cnt = Cnt+1;
                     
                    end
                    %Calculate hours of PMC observation
                   for ii = 1:rows_PMC
                        
                    filename = folders_PMC(ii).name;
                    %retrieve start and end of PMC time
                    [count1_PMC(Cnt_PMC),count2_PMC(Cnt_PMC),startHour1_PMC(Cnt_PMC),endHour1_PMC(Cnt_PMC),startHour2_PMC(Cnt_PMC),endHour2_PMC(Cnt_PMC)]= countHours([path_PMC filename]);
                    %retrieve Zc, backscatter ratio.....
                    [R_MAX(Cnt_PMC),Beta_MAX(Cnt_PMC),Beta_TOT(Cnt_PMC),Z_c(Cnt_PMC),Sigma_RMS(Cnt_PMC)] = retrievePMCData_interannual([path_PMC filename]);
                    
                    %create Matrix to print to txt file
                     DataInfo374_PMC(Cnt_PMC,:) = [year(i) month(i) day(i,k) R_MAX(Cnt_PMC),Beta_MAX(Cnt_PMC),Beta_TOT(Cnt_PMC),Z_c(Cnt_PMC),Sigma_RMS(Cnt_PMC)];
%                      DataInfo372_PMC(Cnt_PMC,:) = [year(i) month(i) day(i,k) count2(Cnt_PMC) startHour2(Cnt_PMC) endHour2(Cnt_PMC)];
                     
                     monthlycnt1_PMC(i) = monthlycnt1_PMC(i)+count1_PMC(Cnt_PMC);
                     monthlycnt2_PMC(i) = monthlycnt2_PMC(i)+count2_PMC(Cnt_PMC);
                     Cnt_PMC = Cnt_PMC+1;
                     
                    end
       
                    
                end
            end
 
       
     

    end
    
count1_total = sum(count1);
count2_total = sum(count2);

end


else 
disp('enter same number of columns for year, month, day')
return


end

% %%Calculate statistic on monthly occurence
if add_season_first_season == false
Percentage_occur_monthly = (monthlycnt1_PMC./monthlycnt1)*100;
tot = sum(monthlycnt1)
pmc= sum(monthlycnt1_PMC)
Percentage_occur_yearly = (sum(monthlycnt1_PMC)./sum(monthlycnt1))*100
%%%%%%%%%%%%%%Plot pie chart
% figure;
% h = pie(Percentage_occur_monthly);
% 
% hText = findobj(h,'Type','text'); % text object handles
% percentValues = get(hText,'String'); % percent values
% txt = {'November: ';'December: ';'January: ';'February: '}; % strings
% combinedtxt = strcat(txt,percentValues); % strings and percent values
% oldExtents_cell = get(hText,'Extent'); % cell array
% oldExtents = cell2mat(oldExtents_cell); % numeric array
% 
% hText(1).String = combinedtxt(1);
% hText(2).String = combinedtxt(2);
% hText(3).String = combinedtxt(3);
% hText(4).String = combinedtxt(4);
% 
% 
% newExtents_cell = get(hText,'Extent'); % cell array
% newExtents = cell2mat(newExtents_cell); % numeric array 
% width_change = newExtents(:,4)-oldExtents(:,4);
% signValues = sign(oldExtents(:,1));
% offset = signValues.*(width_change/2);
% textPositions_cell = get(hText,{'Position'}); % cell array
% textPositions = cell2mat(textPositions_cell); % numeric array
% textPositions(:,1) = textPositions(:,1) + offset; % add offset 
% 
% hText(1).Position = textPositions(1,:);
% hText(2).Position = textPositions(2,:);
% hText(3).Position = textPositions(3,:);
% hText(4).Position = textPositions(4,:);
% 
% t = title(['Season ' num2str(year(2)) '-' num2str(year(3)) ' PMC occurrence']);
% set(t, 'FontSize', 18);
%% Plot bar chart
figure
h= bar(Percentage_occur_monthly,0.5)
t = title(['Season ' num2str(year(2)) '-' num2str(year(3)) ' PMC occurrence']);
set(t, 'FontSize', 18);
txt = {'November ';'December ';'January ';'February '}; % strings
 set(gca,'XTick',1:4,'XTickLabel',txt)
 ylabel('Percentage of PMC occurrence')

 saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/barCharts/fig/season%04d_%04d.fig',year(2),year(3)));          
    print('-djpeg',sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/barCharts/season%04d_%04d',year(2),year(3))); 

else
   INDEX = 1;
   DataInfo374 = xlsread('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/SeasonalAnalysis/McMurdoPMCstat1stSeasontxt.xlsx');
   DataInfo374_PMC = [DataInfo374(:,1) DataInfo374(:,2) DataInfo374(:,3) DataInfo374(:,6) DataInfo374(:,7) DataInfo374(:,9) DataInfo374(:,11) DataInfo374(:,12)];
   year(2) = 2010;
   year(3) = 2011;
end    
    
    %% Plot histograms for season variations
    scrsz = get(0,'ScreenSize');
            figure('Position',[200 70 scrsz(3)*0.3 scrsz(4)*0.9])
            orient portrait; set(gcf,'PaperPositionMode','auto')
            
            %Subplot 1 - Peak Backscatter Ratio
            subplot(5,1,1); histogram( DataInfo374_PMC(:,4),[0:10:200])
            ylabel('Number')
            xlabel('Peak Backscatter Ratio R_{max}')
            t = title([num2str(year(2)) '-' num2str(year(3)) ' Season '])
            set(t, 'FontSize', 18);
            Mean_R_max(INDEX) = mean(DataInfo374_PMC(:,4));
            STD_R_max = std(DataInfo374_PMC(:,4));
            Standard_error_R_max = STD_R_max/sqrt(length(DataInfo374_PMC(:,4)));
            ylim([0 50])
            yl = get(gca,'ylim');
            text(140,yl(2)-10,sprintf(['Mean = %.02f ' char(177) ' %.01f \n Std = %.02f '],Mean_R_max(INDEX),Standard_error_R_max,STD_R_max))
            
             %Subplot 2 - Peak Volume Backscatter  Coefficient \Beta_{max}
            subplot(5,1,2); histogram( DataInfo374_PMC(:,5),[0:0.5:14])
            ylabel('Number')
            xlabel('Peak Volume Backscatter  Coefficient \beta_{max} (10^{-6} sr^{-1})')
            Mean_B_max(INDEX) = mean(DataInfo374_PMC(:,5));
            STD_B_max = std(DataInfo374_PMC(:,5));
            Standard_error_B_max = STD_B_max/sqrt(length(DataInfo374_PMC(:,5)));
            ylim([0 40])
            yl = get(gca,'ylim');
            text(10,yl(2)-10,sprintf(['Mean = %.02f ' char(177) ' %.01f \n Std = %.02f '],Mean_B_max(INDEX),Standard_error_B_max,STD_B_max))
            
             %Subplot 3 - Peak Volume Backscatter  Coefficient \Beta_{tot}
            
            subplot(5,1,3); histogram( DataInfo374_PMC(:,6),[0:20])
            ylabel('Number')
            xlabel('Total Backscatter Coefficient \beta_{Total} (10^{-6} sr^{-1})')
            Mean_B_tot(INDEX) = mean(DataInfo374_PMC(:,6));
            STD_B_tot = std(DataInfo374_PMC(:,6));
            Standard_error_B_tot = STD_B_tot/sqrt(length(DataInfo374_PMC(:,6)));
            ylim([0 50])
            yl = get(gca,'ylim');
            text(14.5,yl(2)-10,sprintf(['Mean = %.02f ' char(177) ' %.01f \n Std = %.02f '],Mean_B_tot(INDEX),Standard_error_B_tot,STD_B_tot))
            
             %Subplot 4 - Centroid Altitude Z_c (km)
            
            subplot(5,1,4); histogram( DataInfo374_PMC(:,7),[81:0.5:89])
            ylabel('Number')
            xlabel('Centroid Altitude Z_c (km)')
            Mean_Zc(INDEX) = mean( DataInfo374_PMC(:,7));
            STD_Zc = std(DataInfo374_PMC(:,7));
            Standard_error_Zc = STD_Zc/sqrt(length(DataInfo374_PMC(:,7)));
            ylim([0 70])
            yl = get(gca,'ylim');
            text(86.8,yl(2)-15,sprintf(['Mean = %.02f ' char(177) ' %.01f \n Std = %.02f '],Mean_Zc(INDEX),Standard_error_Zc,STD_Zc))
         
    
             %Subplot 5 - RMS Width \sigma_{rms} (km)
            
            subplot(5,1,5); histogram( DataInfo374_PMC(:,8),[0:0.1:3])
            ylabel('Number')
            xlabel('RMS Width \sigma_{rms} (km)')
            Mean_Sigma(INDEX) = mean( DataInfo374_PMC(:,8));
            STD_Sigma = std(DataInfo374_PMC(:,8));
            Standard_error_Sigma = STD_Sigma/sqrt(length(DataInfo374_PMC(:,8)));
            ylim([0 40])
            yl = get(gca,'ylim');
            text(2.2,yl(2)-10,sprintf(['Mean = %.02f ' char(177) ' %.02f \n Std = %.02f '],Mean_Sigma(INDEX),Standard_error_Sigma,STD_Sigma))
    
    saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/InterannualAnalysis/FIG/season%04d_%04d.png',year(2),year(3)));          
    print('-djpeg',sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/InterannualAnalysis/season%04d_%04d',year(2),year(3)));
M{INDEX} = DataInfo374_PMC;
D{INDEX} = DataInfo374;
if plot_interannual_variation == true
%       Plot histograms for season variations
    %% Plot histograms for season variations
if plot_overall_Histogram == 1
     clear DataInfo374_PMC 
     DataInfo374_PMC = [M{1};M{2};M{3};M{4};M{5};M{6};M{7};]
     DataInfo374_overall = [D{1}(:,1:6);D{2};D{3};D{4};D{5};D{6};D{7};]
    scrsz = get(0,'ScreenSize');
            figure('Position',[200 70 scrsz(3)*0.3 scrsz(4)*0.9])
            orient portrait; set(gcf,'PaperPositionMode','auto')
            
            %Subplot 1 - Peak Backscatter Ratio
            subplot(5,1,1); histogram( DataInfo374_PMC(:,4),[0:10:200])
            ylabel('Number')
            xlabel('Peak Backscatter Ratio R_{max}')
            t = title('Medium-Strong PMC')
            set(t, 'FontSize', 18);
            Mean_R_max(INDEX) = mean(DataInfo374_PMC(:,4));
            STD_R_max = std(DataInfo374_PMC(:,4));
            Standard_error_R_max = STD_R_max/sqrt(length(DataInfo374_PMC(:,4)));
            yl = get(gca,'ylim');
            text(140,yl(2)-10,sprintf(['Mean = %.02f ' char(177) ' %.01f \n Std = %.02f '],Mean_R_max(INDEX),Standard_error_R_max,STD_R_max))
            
             %Subplot 2 - Peak Volume Backscatter  Coefficient \Beta_{max}
            subplot(5,1,2); histogram( DataInfo374_PMC(:,5),[0:0.5:14])
            ylabel('Number')
            xlabel('Peak Volume Backscatter  Coefficient \beta_{max} (10^{-6} sr^{-1})')
            Mean_B_max(INDEX) = mean(DataInfo374_PMC(:,5));
            STD_B_max = std(DataInfo374_PMC(:,5));
            Standard_error_B_max = STD_B_max/sqrt(length(DataInfo374_PMC(:,5)));
           
            yl = get(gca,'ylim');
            text(10,yl(2)-10,sprintf(['Mean = %.02f ' char(177) ' %.01f \n Std = %.02f '],Mean_B_max(INDEX),Standard_error_B_max,STD_B_max))
            
             %Subplot 3 - Peak Volume Backscatter  Coefficient \Beta_{tot}
            
            subplot(5,1,3); histogram( DataInfo374_PMC(:,6),[0:20])
            ylabel('Number')
            xlabel('Total Backscatter Coefficient \beta_{Total} (10^{-6} sr^{-1})')
            Mean_B_tot(INDEX) = mean(DataInfo374_PMC(:,6));
            STD_B_tot = std(DataInfo374_PMC(:,6));
            Standard_error_B_tot = STD_B_tot/sqrt(length(DataInfo374_PMC(:,6)));
            yl = get(gca,'ylim');
            text(14.5,yl(2)-10,sprintf(['Mean = %.02f ' char(177) ' %.01f \n Std = %.02f '],Mean_B_tot(INDEX),Standard_error_B_tot,STD_B_tot))
            
             %Subplot 4 - Centroid Altitude Z_c (km)
            
            subplot(5,1,4); histogram( DataInfo374_PMC(:,7),[81:0.5:89])
            ylabel('Number')
            xlabel('Centroid Altitude Z_c (km)')
            Mean_Zc(INDEX) = mean( DataInfo374_PMC(:,7));
            STD_Zc = std(DataInfo374_PMC(:,7));
            Standard_error_Zc = STD_Zc/sqrt(length(DataInfo374_PMC(:,7)));
            yl = get(gca,'ylim');
            text(86.8,yl(2)-15,sprintf(['Mean = %.02f ' char(177) ' %.01f \n Std = %.02f '],Mean_Zc(INDEX),Standard_error_Zc,STD_Zc))
         
    
             %Subplot 5 - RMS Width \sigma_{rms} (km)
            
            subplot(5,1,5); histogram( DataInfo374_PMC(:,8),[0:0.1:3])
            ylabel('Number')
            xlabel('RMS Width \sigma_{rms} (km)')
            Mean_Sigma(INDEX) = mean( DataInfo374_PMC(:,8));
            STD_Sigma = std(DataInfo374_PMC(:,8));
            Standard_error_Sigma = STD_Sigma/sqrt(length(DataInfo374_PMC(:,8)));
            yl = get(gca,'ylim');
            text(2.2,yl(2)-10,sprintf(['Mean = %.02f ' char(177) ' %.02f \n Std = %.02f '],Mean_Sigma(INDEX),Standard_error_Sigma,STD_Sigma))
    
    saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/InterannualAnalysis/OverallHistMediumStrong.png',year(2),year(3)));          
   % print('-djpeg',sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/InterannualAnalysis/season%04d_%04d',year(2),year(3))); 

%%Compute overall Mean
R_max_overall = mean(DataInfo374_PMC(:,4))
std_R_max_overall = std(DataInfo374_PMC(:,4))
Standard_error_Sigma_veral = std_R_max_overall /sqrt(length(DataInfo374_PMC(:,4)))

B_max_overall = mean(DataInfo374_PMC(:,5))
std_B_max_overall = std(DataInfo374_PMC(:,5))
Standard_error_Bmax_veral = std_B_max_overall /sqrt(length(DataInfo374_PMC(:,5)))

B_tot_overall = mean(DataInfo374_PMC(:,6))
std_B_tot_overall = std(DataInfo374_PMC(:,6))
Standard_error_Btot_veral = std_B_tot_overall /sqrt(length(DataInfo374_PMC(:,6)))

Zc_overall = mean(DataInfo374_PMC(:,7))
std_Zc_overall = std(DataInfo374_PMC(:,7))
Standard_error_Zc_veral = std_Zc_overall /sqrt(length(DataInfo374_PMC(:,7)))

sigma_overall = mean(DataInfo374_PMC(:,8))
std_sigma_overall = std(DataInfo374_PMC(:,8))
Standard_error_sigma_veral = std_sigma_overall /sqrt(length(DataInfo374_PMC(:,8)))

Percentage_occurrence_overall = length(DataInfo374_PMC)/length(DataInfo374_overall)

end
   
            scrsz = get(0,'ScreenSize');
            figure('Position',[200 70 scrsz(3)*0.3 scrsz(4)*0.9])
            orient portrait; set(gcf,'PaperPositionMode','auto')
            %%Change based on how many seasons of data
            X_labels = ['2010-2011';'2011-2012'; '2012-2013';'2013-2014';'2014-2015';'2015-2016';'2016-2017'] ;
            
             % Subplot 1- in peak backscatter ratio vs time 
             subplot(5,1,1); plot(1:size(X_labels,1),Mean_R_max) 
             hold on
             plot(1:size(X_labels,1),Mean_R_max,'o') 
             set(gca,'XTick',1:7,'XTickLabel',X_labels,'FontSize',8)
             title('R_{max} vs Season Mean') 
             ylabel('R_{max}')
             hold off
             
             % Subplot 2 - Peak Volume Backscatter  Coefficient vs time 
             subplot(5,1,2); plot(1:size(X_labels,1),Mean_B_max) 
             hold on
             plot(1:size(X_labels,1),Mean_B_max,'o') 
             set(gca,'XTick',1:7,'XTickLabel',X_labels,'FontSize',8)
             title('B_{max} vs Season Mean') 
             ylabel('\beta_{max} (10^{-6} sr^{-1})')
             hold off
     
             % Subplot 3 - Total Volume Backscatter  Coefficient vs time 
             subplot(5,1,3); plot(1:size(X_labels,1),Mean_B_tot) 
             hold on
             plot(1:size(X_labels,1),Mean_B_tot,'o') 
            set(gca,'XTick',1:7,'XTickLabel',X_labels,'FontSize',8)
             title('B_{total} vs Season Mean')
             ylabel('\beta_{Total} (10^{-6} sr^{-1})')
             hold off
             
             % Subplot 4 - Centroid Altitude vs time 
             subplot(5,1,4); plot(1:size(X_labels,1),Mean_Zc) 
             hold on
             plot(1:size(X_labels,1),Mean_Zc,'o') 
             set(gca,'XTick',1:7,'XTickLabel',X_labels,'FontSize',8)
             title('Z_c vs Season Mean') 
             ylabel('Z_c (km)')
             hold off
             
             % Subplot 5 - RMS Width vs time 
             subplot(5,1,5); plot(1:size(X_labels,1),Mean_Sigma) 
             hold on
             plot(1:size(X_labels,1),Mean_Sigma,'o') 
             set(gca,'XTick',1:7,'XTickLabel',X_labels,'FontSize',8)
             title('\sigma_{RMS} vs Season Mean') 
             ylabel('\sigma_{RMS}')
             hold off
     Mean_Zc = Mean_Zc';
     Mean_B_tot = Mean_B_tot';
     corrcoef(Mean_Zc,Mean_B_tot')
     
             
    saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/MeanStatisticPlots/mean_season_parameters.jpg'));          
    %print('-djpeg',sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/MeanStatisticPlots/mean_season%04d_%04d',year(2),year(3))); 
       
   end
   
end

% 
% if Seasonal_Analysis == true
% 
% 
% %Open text file for 374nm and 372nm statistics
% name_374 =  sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/ObservationHours/%s374_hours.txt',smonth(month(1,:),:));
% name_372 =  sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/ObservationHours/%s372_hours.txt',smonth(month(1,:),:));
% fpout_374 = fopen(name_374,'w');
% fpout_372 = fopen(name_372,'w');
% name_374_PMC =  sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/ObservationHours/%s374_PMC_hours.txt',smonth(month(1,:),:));
% name_372_PMC =  sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/ObservationHours/%s372_PMC_hours.txt',smonth(month(1,:),:));
% fpout_374_PMC = fopen(name_374_PMC,'w');
% fpout_372_PMC = fopen(name_372_PMC,'w');
% %set count
% count1 = zeros(length(year),1);
% count2 = zeros(length(year),1);
% if yc == mc && mc == dr
% 
% %Give path to folder
% path = '/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/FeDensityData/';
% path_PMC = '/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/pmcData/';
% 
% %initialize count
% Cnt = 1;
% Cnt_PMC =1;
%     for i = 1:yc 
%             monthlycnt1(i) = 0;
%             monthlycnt2(i) = 0;
%             monthlycnt1_PMC(i) = 0;
%             monthlycnt2_PMC(i) = 0;
%             
%        
%             for k = 1:dc
%                 if day(i,k)~= 0 && month(i)~= 0
%                     %Get number of files in both folders
%                     folders = dir([path '*'  num2str(year(i)) smonth(month(i),:) num2str(day(i,k),'%02i') '*']);
%                     folders_PMC = dir([path_PMC '*'  num2str(year(i)) smonth(month(i),:) num2str(day(i,k),'%02i') '*']);
%                     if length(folders)== 0 
%                       fprintf('%d/%d/%d does not have any data \n',year(i),month(i), day(i,k))
%                     end
%                     [rows,cols] = size(folders);
%                     [rows_PMC,cols_PMC] = size(folders_PMC);
%                     %Calculate hours of observation
%                     for ii = 1:rows
%                         
%                     filename = folders(ii).name;
%                     [count1(Cnt),count2(Cnt),startHour1(Cnt),endHour1(Cnt),startHour2(Cnt),endHour2(Cnt)]= countHours([path filename]);
%                     %create Matrix to print to txt file
%                      DataInfo374(Cnt,:) = [year(i) month(i) day(i,k) count1(Cnt) startHour1(Cnt) endHour1(Cnt)];
%                      DataInfo372(Cnt,:) = [year(i) month(i) day(i,k) count2(Cnt) startHour2(Cnt) endHour2(Cnt)];
%                      
%                      monthlycnt1(i) = monthlycnt1(i)+count1(Cnt);
%                      monthlycnt2(i) = monthlycnt2(i)+count2(Cnt);
%                      Cnt = Cnt+1;
%                      
%                     end
%                     %Calculate hours of PMC observation
%                    for ii = 1:rows_PMC
%                         
%                     filename = folders_PMC(ii).name;
%                     [count1_PMC(Cnt_PMC),count2_PMC(Cnt_PMC),startHour1_PMC(Cnt_PMC),endHour1_PMC(Cnt_PMC),startHour2_PMC(Cnt_PMC),endHour2_PMC(Cnt_PMC)]= countHours([path_PMC filename]);
%                     %create Matrix to print to txt file
%                      DataInfo374_PMC(Cnt_PMC,:) = [year(i) month(i) day(i,k) count1_PMC(Cnt_PMC) startHour1_PMC(Cnt_PMC) endHour1(Cnt_PMC)];
%                      DataInfo372_PMC(Cnt_PMC,:) = [year(i) month(i) day(i,k) count2(Cnt_PMC) startHour2(Cnt_PMC) endHour2(Cnt_PMC)];
%                      
%                      monthlycnt1_PMC(i) = monthlycnt1_PMC(i)+count1_PMC(Cnt_PMC);
%                      monthlycnt2_PMC(i) = monthlycnt2_PMC(i)+count2_PMC(Cnt_PMC);
%                      Cnt_PMC = Cnt_PMC+1;
%                      
%                     end
%        
%                     
%                 end
%             end
%  
%        
%      
% 
%     end
%     
% count1_total = sum(count1);
% count2_total = sum(count2);
% 
% 
% else 
% disp('enter same number of columns for year, month, day')
% return
% 
% 
% end
% 
% %Calculate statistic on monthly occurence
% Percentage_occur_monthly = (monthlycnt1_PMC./monthlycnt1)*100;
% 
% 
% 
% %Write data to txt file
% fprintf(fpout_374,'Year    Month    Day    Tot_hour    startHour    endHour \n');
% fprintf(fpout_372,'Year    Month    Day    Tot_hour    startHour    endHour \n');
% fprintf(fpout_374_PMC,'Year    Month    Day    Tot_hour    startHour    endHour \n');
% fprintf(fpout_372_PMC,'Year    Month    Day    Tot_hour    startHour    endHour \n');
% 
% dlmwrite(name_374,DataInfo374,'-append','delimiter',' ');
% dlmwrite(name_372,DataInfo372,'-append','delimiter',' ');
% dlmwrite(name_374_PMC,DataInfo374_PMC,'-append','delimiter',' ');
% dlmwrite(name_372_PMC,DataInfo372_PMC,'-append','delimiter',' ');
% 
% fclose(fpout_374);
% fclose(fpout_372);
% fclose(fpout_374_PMC);
% fclose(fpout_372_PMC);
% 
% end







