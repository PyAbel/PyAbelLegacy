close all;
clear all;
%raw = load('e:\dimer_scan\.asc');
%H=max(raw(:,1))-min(raw(:,1))+1;
%W=max(raw(:,2))-min(raw(:,2))+1;
%ims = reshape(raw(:,3),H,W);
ims=load('e:\dimer_scan\ion_images\processing\303eci1.txt');
ims=ims';
[W H]=size(ims);
mindel=10^100;
i=1;
w=300;
for c=fix(W/2)-w:fix(W/2)+w
    Im=ims((c-200):(c+200),(fix(H/2)-100):(fix(H/2)+100));
    Iminv=ims((c+200):-1:(c-200),(fix(H/2)-100):(fix(H/2)+100));
    delta=abs(Iminv-Im);
    del=median(sum(delta));
    %%%%%%%%%%%%%%%%%%%%
    m1(i)=del;
    c1(i)=c;
    i=i+1;
    %%%%%%%%%%%%%%%%%%%%
    if del<mindel
       mindel=del;
       xc=c;
    end;
end;
xc
mindel=10^100;
i=1;
for c=fix(H/2)-w:fix(H/2)+w
    Im1=ims((xc-100):(xc+100),(c-200):(c+200));
    Iminv1=ims((xc-100):(xc+100),(c+200):-1:(c-200));
    delta=abs(Iminv1-Im1);
    del=median(sum(delta'));
    %%%%%%%%%%%%%%%%%%%%
    m2(i)=del;
    c2(i)=c;
    i=i+1;
    %%%%%%%%%%%%%%%%%%%%
    if del<mindel
       mindel=del;
       zc=c;
    end;
end;
zc
figure
plot(c1,m1/max(m1),'o');
hold on;
plot(c2,m2/max(m2),'ro');
title('Xc(blue) & Zc(red)');

xc=input('Input the center point in horizontal direction xc=');
zc=input('Input the center point in vertical direction zc=');
outsizeX=input('Input the size (in pixels) of the output file in horizontal direction: ');
outsizeZ=input('Input the size (in pixels) of the output file in vertical direction: ');

IMobrez=zeros(outsizeZ,outsizeX);
IMobrez=(ims(xc-fix(outsizeX/2):xc+fix(outsizeX/2),zc-fix(outsizeZ/2):zc+fix(outsizeZ/2)))';
f=fopen('e:\dimer_scan\ion_images\processing\303eci1.asc','w');
for i=1:outsizeZ
    for j=1:outsizeX
        fprintf(f,'%d,  %d  %2.5f\n',j,i,IMobrez(i,j));
    end;
end;
fclose(f); 

return
