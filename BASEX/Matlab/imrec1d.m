close all;
clear all;

fname=input('Input the name of the data file: ','s');

Vsyms=input('Symmetrize the image with respect to vertical axis? Y/N :  ','s');
Vsymv=0;Vsymv=(Vsyms=='Y')|(Vsyms=='y');

Hsyms=input('Symmetrize the image with respect to horizontal axis? Y/N :  ','s');
Hsymv=0;Hsymv=(Hsyms=='Y')|(Hsyms=='y');

raw = load(fname);
ww = length(raw(1,:));
if  ww == 3 
    H=max(raw(:,1))-min(raw(:,1))+1;
    W=max(raw(:,2))-min(raw(:,2))+1;
    ims = reshape(raw(:,3),H,W);
    [W H]=size(ims);
    ims = ims';
else 
    ims=raw;
    [H W]=size(raw);
end;

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


N = 1001;
N_2=fix(N/2);
W_2=fix(W/2);
H_2=fix(H/2);
im = zeros(N,N);
im(N_2+1-H_2:N_2+H_2+1,N_2+1-W_2:N_2+W_2+1)=ims;

M = load('basis1000pr_1.bst');
Mc = load('basis1000_1.bst');
NBF=size(M,2);

RANK = NBF-0;
[u s v]=svd(M);
ds = diag(s);
dsi = [1./ds(1:RANK); zeros(NBF-RANK,1)];
si = [diag(dsi);zeros(N-NBF,NBF)];
invM=(u*si*v')';

E = zeros(NBF,NBF);
for n=1:NBF
   E(n,n)=5;
end;

Ci=(im)*M*inv(M'*M+E);
P = Ci*M';
IM=Ci*Mc';

IMobrez=IM(N_2+1-H_2:N_2+1+H_2,N_2+1-W_2:N_2+1+W_2);
IMobrez = IMobrez.*(IMobrez>=0);
RES=ims-P(N_2+1-H_2:N_2+1+H_2,N_2+1-W_2:N_2+1+W_2);
if W_2<=H_2
   vel=zeros(W_2,1);
   for m=1:W_2-1
      for k=1:m
         z=m*cos(pi*k/(2*m));
         x=m*sin(pi*k/(2*m));
         vel(m)=vel(m)+(z-fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+2,W_2+fix(x)+2)*x+...
         (1-z+fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+2)*x+...
         (z-fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+2,W_2+fix(x)+1)*x+...
         (1-z+fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+1)*x;
      end;
      X=[zeros(90,1) ones(90,1)];
      Y=zeros(90,1);
      for k=1:90
          X(k,1)=1.5*(cos(pi*k/180))^2-0.5;
          z=m*cos(pi*k/180);
          x=m*sin(pi*k/180);
          Y(k)=(z-fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+2,W_2+fix(x)+2)+...
         (1-z+fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+2)+...
         (z-fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+2,W_2+fix(x)+1)+...
         (1-z+fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+1);
      end;
      b=regress(Y,X);
      beta(m)=b(1)/b(2);
      if beta(m)>2
          beta(m)=2;
      end;
      if beta(m)<-1
          beta(m)=-1;
      end;
   end;
else
   vel=zeros(H_2,1);
   for m=1:H_2-1
      for k=1:m
         z=m*cos(pi*k/(2*m));
         x=m*sin(pi*k/(2*m));
         vel(m)=vel(m)+(z-fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+2,W_2+fix(x)+2)*x+...
         (1-z+fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+2)*x+...
         (z-fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+2,W_2+fix(x)+1)*x+...
         (1-z+fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+1)*x;
      end;
         X=[zeros(90,1) ones(90,1)];
      Y=zeros(90,1);
      for k=1:90
          X(k,1)=1.5*(cos(pi*k/180))^2-0.5;
          z=m*cos(pi*k/180);
          x=m*sin(pi*k/180);
          Y(k)=(z-fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+2,W_2+fix(x)+2)+...
         (1-z+fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+2)+...
         (z-fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+2,W_2+fix(x)+1)+...
         (1-z+fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+1);
      end;
      b=regress(Y,X);
      beta(m)=b(1)/b(2);
      if beta(m)>2
          beta(m)=2;
      end;
      if beta(m)<-1
          beta(m)=-1;
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
hold on;
plot(beta,'r');
title('Speed distribution(blue) & anisotropy(red)');
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
fvel=fullfile(pathstr,[name '_spd.dat'],versn);
fanis=fullfile(pathstr,[name '_anis.dat'],versn);
if ww == 3
f=fopen(fname,'w');
for i=1:H
    for j=1:W
        fprintf(f,'%d,  %d  %2.5f\n',j,i,IMobrez(i,j));
    end;
end;
fclose(f); 
else 
    save(fname,'IMobrez', '-ASCII');
 end;
 
 m=(1:length(vel))';
 aux=[m vel(:)];
 save(fvel,'aux', '-ASCII');

 m=(1:length(beta))';
 aux=[m beta(:)];
 save(fanis,'aux', '-ASCII');
return

