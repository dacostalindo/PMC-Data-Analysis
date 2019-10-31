function F=HarmonicFitting(x,xdata)
F=x(1)+x(2).*cos(2*pi/12.*(xdata-x(3)))+x(4).*cos(2*pi/24.*(xdata-x(5)));