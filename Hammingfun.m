function[w]=Hammingfun(FW,x0,x)

w=0.54-0.46*cos(2*pi*(x-x0+FW/2)/FW);

% Here FW is the full-wdith for the hamming smooth;
% x0 is the center of the data point where the smoothed result is assigned to
% x is the altitude or time as the x-axis

% the original Hamming function is defined as
% w[k+1]=0.54-0.46cos(2pi*k/(n-1)), k=0,1,...,n-1