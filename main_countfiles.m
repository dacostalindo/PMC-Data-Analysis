clear all
close all
clc
%PURPOSE: count the amount of datafiles in a whole month
%change 
%NOTE: since matrices need to have same number of columns and rows, fill in
%the days and months and then fill with 0s the rest.
%%%Example:
% % %NOVEMBER OF EACH YEAR
year = [ 2011 2012 2013 2014 2015]; month = [ 11 11 11 11 11]; day =  [ 17 18 0 28 29 30 0 0 0 0 0 0 0; ... 
                                                                            0 5 7 8 9 10 16 17 18 19 22 23 24; ...
                                                                            12 14 17 18 19 20 26 29 0 0 0 0 0 ; ...
                                                                            1 2 3 10 11 15 16 21 0 0 0 0 0; ...
                                                                            8 9 10 26 27 30 0 0 0 0 0 0 0 ];
   
% % % % % DECEMBER OF EACH YEAR
% year = [2011 2012 2013 2014 2015 2016]; month = [ 12; 12; 12; 12; 12; 12]; day =  [ 11 14 15 19 20 21 23 24 25 26 27 29 30 31 00 00; ... 
%                                                                                     08 09 17 27 28 29 31 00 00 00 00 00 00 00 00 00; ...
%                                                                                     01 09 00 11 12 22 23 24 25 27 28 29 00 31 00 00; ...
%                                                                                     02 04 05 08 09 10 11 12 13 17 18 19 20 28 00 31; ...
%                                                                                     01 02 03 17 22 24 26 27 28 00 00 00 00 00 00 00; ...
%                                                                                     03 04 12 13 14 00 00 00 00 00 00 00 00 00 00 00];
% % 
% % %    %FEBRUARY OF EACH YEAR
%    year = [2012 2013 2014 2015]; month = [ 02; 02; 02; 02;]; day =  [               01 03 04 07 08 09 15 17 20 21 23 24; ... 
%                                                                                     09 20 24 26 00 00 00 00 00 00 00 00; ...
%                                                                                     05 06 10 11 12 21 22 26 27 00 00 00; ...
%                                                                                     11 12 14 16 17 18 26 00 00 00 00 00];
% % % %       
% %JANUARY OF EACH YEAR
% year = [2012 2013 2014 2015 2016 2017]; month = [ 01; 01; 01; 01; 01; 01;]; day =  [05 06 10 11 12 15 16 17 26 29 30 31 00 00 00 00;  ... 
%                                                                                     01 02 03 04 09 10 11 21 22 23 24 25 29 00 00 00; ...
%                                                                                     10 14 19 20 27 28 00 00 00 00 00 00 00 00 00 00; ...
%                                                                                     01 11 21 22 23 26 27 28 29 30 00 00 00 00 00 00; ...
%                                                                                     05 06 13 14 15 16 22 23 24 25 00 00 00 00 00 00; ...
%                                                                                     02 03 07 10 11 16 18 19 20 24 25 26 28 29 30 31 ];    
%    
                                                                               
                                                                                
                                                                        

%make a vector containing initials for each month                                                                
smonth = ['JA';'FB';'MR';'AR';'MY';'JN';'JL';'AG';'SP';'OT';'NV';'DC'];                                                                        
%take vector sizes for reading in data purposes
[yr,yc] = size(year);
[mr,mc] = size(month);
[dr,dc] = size(day);

%set count
count1 = zeros(length(year),1);
count2 = zeros(length(year),1);
if yc == mc && mc == dr 
    for i = 1:yc 
        
            for k = 1:dc
                if day(i,k)~= 0 && month(i)~= 0
                    %Make filelist
                    path = ['/Users/asta10/Desktop/Research/SkyObservations/' num2str(year(i)) smonth(month(i),:) num2str(day(i,k),'%02i')];
                    [count11,count22]= countfiles(path,1,1);
                    %Count total
                    count1(i) = count1(i)+ count11;
                    count2(i) = count2(i) +count22;
                end
            end
    end
    
count1_total = sum(count1);
count2_total = sum(count2);


else 
disp('enter same number of columns for year, month, day')
return


end






