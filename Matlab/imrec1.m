%function imrec(fname)

close all;
clear all;

f=fopen('path.txt','r');
fname=fscanf(f,'%s');
Vsyms=input('Symmetrize the image with respect to the vertical axis? Y/N :  ','s');
Vsymv=0;Vsymv=(Vsyms=='Y')|(Vsyms=='y');

Hsyms=input('Symmetrize the image with respect to the horizontal axis? Y/N :  ','s');
Hsymv=0;Hsymv=(Hsyms=='Y')|(Hsyms=='y');

f=fopen(fname,'rb');
fseek(f,-8,'eof');
SZ=fread(f,2,'float');
fseek(f,0,'bof');
ims = reshape(fread(f,prod(SZ),'float'),SZ(1),SZ(2))';
fclose(f);

[H W]=size(ims);

if(Vsymv==1)
    W_2=fix(W/2)+1;
    ims(:,1:W_2)=(ims(:,1:W_2)+ims(:,W:-1:W_2))/2;
    ims(:,W_2:W)=ims(:,W_2:-1:1);
end;

if(Hsymv==1)
    H_2=fix(H/2)+1;
    ims(1:H_2,:)=(ims(1:H_2,:)+ims(H:-1:H_2,:))/2;
    ims(H_2:H,:)=ims(H_2:-1:1,:);
end;

f = fopen('basissize.txt');
aux=fscanf(f,'%d %d');
N=aux(1);NBF=aux(2);
fclose(f);

f=fopen('basisset.bin','rb');
Mc=fread(f,[1001 NBF],'float');
fclose(f);
f=fopen('basisset_pr.bin','rb');
M=fread(f,[1001 NBF],'float');
fclose(f);

N_2=fix(N/2);
W_2=fix(W/2);
H_2=fix(H/2);
im = zeros(N,N);
im(N_2+1-H_2:N_2+H_2+1,N_2+1-W_2:N_2+W_2+1)=ims;

NBFc=NBF;
NBF=NBF-1;
M1=M;
M(:,1:end-1)=M(:,2:end);
M(:,end)=[];

RANK = NBF-0;
[u s v]=svd(M);
ds = diag(s);
dsi = [1./ds(1:RANK); zeros(NBF-RANK,1)];
si = [diag(dsi);zeros(N-NBF,NBF)];
invM=(u*si*v')';


RANK = NBFc-0;
[u s v]=svd(Mc);
ds = diag(s);
dsi = [1./ds(1:RANK); zeros(NBFc-RANK,1)];
si = [diag(dsi);zeros(N-NBFc,NBFc)];
invMc=(u*si*v')';


[u s v]=svd(M1);
ds = diag(s);
dsi = [1./ds(1:RANK); zeros(NBFc-RANK,1)];
si = [diag(dsi);zeros(N-NBFc,NBFc)];
invM1=(u*si*v')';

return
Ci=invMc*(im)*invM';
Ci1=invMc*(im)*invM1';

P = Mc*Ci1*M1';
Mcc=Mc(:,2:end);
IM=Mc*Ci*Mcc';
IM(:,N_2+1)=IM(:,N_2);

IMobrez=IM(N_2+1-H_2:N_2+1+H_2,N_2+1-W_2:N_2+1+W_2);
IMobrez = IMobrez.*(IMobrez>=0);
RES=ims-P(N_2+1-H_2:N_2+1+H_2,N_2+1-W_2:N_2+1+W_2);

if W_2<=H_2
   vel=zeros(W_2,1);
   for m=1:W_2
      for k=1:m
         z=sqrt(m^2-k^2);
         vel(m)=vel(m)+(z-fix(z))*IMobrez(H_2+fix(z)+2,W_2+1+k)*k+...
         (1-z+fix(z))*IMobrez(H_2+fix(z)+1,W_2+1+k)*k;
      end;
   end;
else
   vel=zeros(H_2,1);
   for m=1:H_2
      for k=1:m
         z=sqrt(m^2-k^2);
         vel(m)=vel(m)+(z-fix(z))*IMobrez(H_2+fix(z)+2,W_2+1+k)*k+...
         (1-z+fix(z))*IMobrez(H_2+fix(z)+1,W_2+1+k)*k;
      end;
   end;
end;

vel=vel/max(vel);
figure
subplot(2,2,1);
imagesc(IMobrez);
title('Reconstructed image');
axis square
colorbar
subplot(2,2,2);
plot(vel);
title('Speed distribution');

axis square
subplot(2,2,3);
imagesc(ims);
title('Original image (projection)');
axis square
colorbar
subplot(2,2,4);
imagesc(RES);
title('Residuals of projection expansion');

axis square
colorbar

[pathstr,name,ext,versn] = fileparts(fname);
fname=fullfile(pathstr,[name '_3D' ext],versn);
fvel=fullfile(pathstr,[name '_vel.dat'],versn);


disp('Saving results.');

f=fopen(fname,'w');
for i=1:H
    for j=1:W
        fprintf(f,'%d,  %d  %2.5f\n',j,i,IMobrez(i,j));
    end;
end;
fclose(f); 
 
f = fopen(fvel,'w');
for i=1:length(vel)
    fprintf(f,'%d %d\n',i,vel(i));    
end;
fclose(f);
disp('Done');
disp('Close the figure to exit the program.');


return

