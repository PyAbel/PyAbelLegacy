close all;
clear all;

img=zeros(451,451);
img(126,126)=2*sqrt(100.5^2-10000);
img(126,326)=2*sqrt(100.5^2-10000);
img(326,126)=2*sqrt(100.5^2-10000);
img(326,326)=2*sqrt(100.5^2-10000);
for i=127:325
   img(126,i)=2*(sqrt(100.5^2-(i-226)^2)-sqrt(99.5^2-(i-226)^2));
   img(326,i)=2*(sqrt(100.5^2-(i-226)^2)-sqrt(99.5^2-(i-226)^2));
end;

f=fopen('f:\TestImages\pr4pnts.asc','w');
for i=1:451
    for j=1:451
        fprintf(f,'%d,  %d  %2.5f\n',j,i,img(i,j));
    end;
end;
fclose(f);

return