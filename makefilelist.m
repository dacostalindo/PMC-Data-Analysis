function [ ] = makefilelist(pathToDay,rx1exists,rx2exists)
%function creates a text file called "filelist.txt" containing all list of
%name of all the files in the folder. Pass 1 for rx1exists and rx2exists if
%both folders exist.
%EXAMPLE: makefilelist('/Users/asta10/Desktop/Research/SkyObservations/2012JA31',1,1)
% 

%List files in RX1

if rx1exists ==1
rx1 = [pathToDay '/RX1/'];
  count = 0;
  year = pathToDay(end-7:end-4);
  folders = dir([rx1 '*' year '*']);
  
  
  for i = 1:length(folders)
 if folders(i).bytes > 100
     count = count+1;
 end
  end
  
  
  fpout = fopen([rx1 'filelist.txt'],'w');
  fprintf(fpout,'%d \n',count);
  
 
 for i = 1: length(folders)
     if folders(i).bytes >= 100
 fprintf(fpout,[folders(i).name '\n']);
     end
 end

fclose(fpout);  
end

%List files in RX2


if rx2exists ==1
rx2 = [pathToDay '/RX2/'];
  count = 0;
  year = pathToDay(end-7:end-4);
  folders = dir([rx2 '*' year '*']);
  

  for i = 1:length(folders)
 if folders(i).bytes > 100
     count = count+1;
 end
  end
  
  
 
 fpout = fopen([rx2 'filelist.txt'],'w');
  fprintf(fpout,'%d \n',count);
 for i = 1: length(folders)
     if folders(i).bytes > 100
 fprintf(fpout,[folders(i).name '\n']);
     end
end
fclose(fpout);  

end
end



