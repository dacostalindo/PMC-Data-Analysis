% launch GUI
 figureH = interactive(year, month, day,BinNum,Rawbinwid,Normdpuresig374,time374,AltCorrection,Beta_PMC,Altl,Alth,dBeta,FINDHfl,FINDHfh,PMCyn,BmaxLG0,BFWHMLG0,BLG1Err,BLG2Err,BmaxLG2Err,BmaxLGnoise,ZBpk,FWHM,Beta_max,ZHMlbin,ZHMhbin);
 % get the handle to the GUI and wait for it to be closed
 hGui = findobj('Tag','tformChoiceGui');
 waitfor(hGui);