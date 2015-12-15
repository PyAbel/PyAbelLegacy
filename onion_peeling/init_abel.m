function [val1,val2]=init_abel(xcg,ycg)
xc = double(xcg);
yc = double(ycg);
val1=zeros(xc+1,xc+1);
val2=zeros(xc+1,yc+1);

for ii = 1: xc+1
    for  jj = ii: xc+1
        val1(ii,jj)=asin(ii/jj)-asin((ii-1)/jj);
    end
end
for idist = 1:xc+1
    for  jdist = 1:yc+1
        val2(idist,jdist)=sqrt(idist^2+(jdist-1)^2)/idist;
    end
end
end