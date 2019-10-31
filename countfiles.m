function [count1,count2] = countfiles(pathToDay,rx1exists,rx2exists)
%function creates a text file called "filelist.txt" containing all list of
%name of all the files in the folder. Pass 1 for rx1exists and rx2exists if
%both folders exist.
%EXAMPLE: makefilelist('/Users/asta10/Desktop/Research/SkyObservations/2012JA31',1,1)
% 

if rx1exists ==1
rx1 = [pathToDay '/RX1/GoodfilelistUpdated.txt'];
C = textread(rx1, '%s','delimiter', '\n');
count1 = str2num(C{1});
%List files in RX2


if rx2exists ==1
rx2 = [pathToDay '/RX2/GoodfilelistUpdated.txt'];
C = textread(rx2, '%s','delimiter', '\n');
count2 = str2num(C{1});
%List files in RX2
  
 

end


end

% %List files in RX1
% 
% if rx1exists ==1
% rx1 = [pathToDay '/RX1/'];
%   count1 = 0;
%   year = pathToDay(end-7:end-4);
%   folders = dir([rx1 '*' year '*']);
%   
%   
%   for i = 1:length(folders)
%  if folders(i).bytes > 100
%      count1 = count1+1;
%  end
%   end
%   
%   
%  
% %List files in RX2
% 
% 
% if rx2exists ==1
% rx2 = [pathToDay '/RX2/'];
%   count2 = 0;
%   year = pathToDay(end-7:end-4);
%   folders = dir([rx2 '*' year '*']);
%   
% 
%   for i = 1:length(folders)
%  if folders(i).bytes > 100
%      count2 = count2+1;
%  end
%   end
%   
%   
%  
% 
% end

