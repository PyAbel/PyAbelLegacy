close all;
clear all;
raw = load('e:\basex2\data.asc');
H=max(raw(:,1))-min(raw(:,1))+1;
W=max(raw(:,2))-min(raw(:,2))+1;
ims = reshape(raw(:,3),H,W);
[W H]=size(ims);
mindel=10^100;
i=1;
w=20;
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
    Im1=ims((fix(W/2)-100):(fix(W/2)+100),(c-200):(c+200));
    Iminv1=ims((fix(W/2)-100):(fix(W/2)+100),(c+200):-1:(c-200));
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

outsize=421;
IMobrez=zeros(outsize,outsize);
IMobrez=(ims(xc-fix(outsize/2):xc+fix(outsize/2),zc-fix(outsize/2):zc+fix(outsize/2)))';
f=fopen('i:\NO_dimer\220nm_reprocessed\ctemp.asc','w');
for i=1:outsize
    for j=1:outsize
        fprintf(f,'%d,  %d  %2.5f\n',j,i,IMobrez(i,j));
    end;
end;
fclose(f); 

return
