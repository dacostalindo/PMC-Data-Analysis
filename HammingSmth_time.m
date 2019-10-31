function[ysmd]=HammingSmth(FW,x0,x,y)
%
for nn=1:length(x0)
   ysmd(nn)=0;
   wsum(nn)=0;
   count(nn)=0;
   for jj=1:length(x)
      if (abs(x(jj)-x0(nn))<=FW/2 & ~isnan(y(jj))==1)   
            % for temperature, using the data points that are not NaN
            weight(jj)=Hammingfun(FW,x0(nn),x(jj));
			ysmd(nn)=ysmd(nn)+y(jj)*weight(jj);
            wsum(nn)=wsum(nn)+weight(jj);
            count(nn) = count(nn)+1;
      end
   end
   if (wsum(nn)~=0) && count(nn) >= (FW/2)*60
		ysmd(nn)=ysmd(nn)/wsum(nn);
   else
       ysmd(nn)=NaN;
   end
end