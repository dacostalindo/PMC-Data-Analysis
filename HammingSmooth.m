function[y]=HammingSmooth(x,m)

% Hamming window smooth, where m is bin # for FWHM smoothing window

hamming= .54 - .46*cos(2*pi*(0:m*2)'/(m*2));
w=hamming/sum(hamming);	% Normalized Hamming weights (it is a column)
%w=userhamming(m*2+1)/sum(userhamming(m*2+1));	% Normalized Hamming weights (it is a column)
y=x;
for i=m+1:length(x)-m
	y(i)=x(i-m:i+m)*w;
end
