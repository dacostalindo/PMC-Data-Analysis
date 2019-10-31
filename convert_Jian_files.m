load('/Users/asta10/Desktop/profile_2011JA28_RX2.mat');
sz = size(FileList);

%mkdir '/Users/asta10/Desktop/Research/SkyObservations/2011JA28'
%mkdir '/Users/asta10/Desktop/Research/SkyObservations/2011JA28/RX2'

 

for i = 1:sz(2)
    
    name =  sprintf('/Users/asta10/Desktop/temp/%s', FileList{i});
    fpout = fopen(name,'w');
    fprintf(fpout,'Date: \n');
    
    if FileList{i}(6:7) == 'FB'
        k = 2;
    elseif FileList{i}(6:7) == 'JA'
        k =1;
    end
    
      fprintf(fpout,'%02d %02d %04d \n',str2num(FileList{i}(8:9)),k,str2num(FileList{i}(2:5)));
   
      fprintf(fpout,'Day of Year:\n');
      fprintf(fpout,'%d \n',DOY);
      
      fprintf(fpout,'Time:\n');
      fprintf(fpout,'%f \n',UTHour(i)); 
      
      fprintf(fpout,'Base Altitude: \n0.20 \n');
      
      fprintf(fpout,'Bin Width:\n');
      fprintf(fpout,'%d\n',BinWidth(i));
      
      fprintf(fpout,'Shots:\n');
      fprintf(fpout,'%d\n',ShotNum(i));

      fprintf(fpout,'Bin Numbers:\n');
      fprintf(fpout,'4096\n');
      
      fprintf(fpout,'Wavelength:\n');
      fprintf(fpout,'374.00000\n');

      fprintf(fpout,'Power of Laser:\n');
      fprintf(fpout,'0.00\n');
      
      fprintf(fpout,'Latitude:\n');
      fprintf(fpout,'0.000\n');

      fprintf(fpout,'Longitude:\n');
      fprintf(fpout,'0.000\n');

      fprintf(fpout,'Azimuth:\n');
      fprintf(fpout,'0.0\n');

      fprintf(fpout,'Elevation:\n');
      fprintf(fpout,'90.0\n');
      
       fprintf(fpout,'\n \n');
      
      for j= 1:4096
      fprintf(fpout,[num2str(PhotonProfile(i,j)) '\n']);
      end


    
    
    fclose(fpout);
end

folders = dir('/Users/asta10/Desktop/temp/')
