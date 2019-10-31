%plot first season only
clearvars -except Mean_R_max Mean_B_max Mean_B_tot Mean_Zc Mean_Sigma smdBTOT X smdZc smdRMS PMC_MATRIX MATRIX M D
close all
clc


I = 1;
season1_PMC_Matrix = xlsread('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/SeasonalAnalysis/McMurdoPMCstat1stSeasontxt.xlsx');
PMC_MATRIX{I} = season1_PMC_Matrix;

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
            plot(timeday,yfit,'-k', 'LineWidth',2);
            plot(timeday,yfit1,'--r','LineWidth',2)
             
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
             plot(timeday,yfit,'-k', 'LineWidth',2);
             plot(timeday,yfit1,'--r','LineWidth',2)
             
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
             plot(timeday,yfit,'-k', 'LineWidth',2);
             plot(timeday,yfit1,'--r','LineWidth',2)
             
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
             plot(timeday,yfit,'-k', 'LineWidth',2);
             plot(timeday,yfit1,'--r','LineWidth',2)
             
             hold off
             ylabel('Z_c (km)')
             xlabel('UT\_Hour')
             
             
x0=[f1.a0,sqrt(f1.a1^2+f1.a2^2),UT12_R_checka,sqrt(f1.b1^2+f1.b2^2),UT24_R];
[beta1,resnorm,resid,exitflag,output,lambda,J]=lsqcurvefit(@HarmonicFitting,x0,timeday,Zc_mean);
disp(beta1);
ci1=nlparci(beta1,resid,'jacobian',J)
error1=(ci1(:,2)-ci1(:,1))/2;

             
 saveas(gcf,sprintf('/Users/asta10/Desktop/Research/Results/McMurdo/ProcessedData/DiurnalAnalysis/mean_diurnal_parameters.jpg'));          
   
    
    