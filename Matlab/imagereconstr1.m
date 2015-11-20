close all;
clear all;
%%%%%%%%% Reading data file (ASCII-2 format) and other input parameters %%%%%%%%%%
fname=input('Input the name of the data file: ','s');
q=input('Choose the regularization parameter:');
step=input('Choose step size for speed distribution calculations:');
step1=input('Choose step size for angular distribution calculations:');
%%%%%%%%%%%%%%% Done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Converting image into a matrix %%%%%%%%
raw = load(fname);
H = 451;
W = 451;
ims = raw;
%ww = length(raw(1,:));
%H=max(raw(:,1))-min(raw(:,1))+1;
%W=max(raw(:,2))-min(raw(:,2))+1;
%ims = reshape(raw(:,3),H,W);
%[W H]=size(ims);
%ims = ims';
%%%%%%%%%%%%%% Done %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Appending the image with zeroes up to the size of 1001*1001 pixels %%%%%%%%% 
N = 1001;
N_2=fix(N/2);
W_2=fix(W/2);
H_2=fix(H/2);
im = zeros(N,N);
im(N_2+1-H_2:N_2+H_2+1,N_2+1-W_2:N_2+W_2+1)=ims;
%%%%%%%%%%%%%%%% Done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Loading basis sets %%%%%%%%%%%
M = load('basis1000pr_1.txt');
Mc = load('basis1000_1.txt');
sigma=1;
NBF=size(M,2);
%%%%%%%%%%%%%%%% Done %%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% Creating diagonal matrix for regularization %%%%%%%%%%%%%%
E = zeros(NBF,NBF);
for n=1:NBF
   E(n,n)=q;
end;
%%%%%%%%%%%%%%%%% Done %%%%%%%%%%%

%%%%%%%%%%%%%%% Reconstructing image %%%%%%%%%%%%%%%%%%%%%%%
Ci=inv(Mc'*Mc)*Mc'*(im)*M*inv(M'*M+E);
P = Mc*Ci*M';
IM=Mc*Ci*Mc';
%%%%%%%%%%%%%%%%%%% Done %%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Cutting the image to the size of the original data %%%%%%%%%%% 
IMcut=IM(N_2+1-H_2:N_2+1+H_2,N_2+1-W_2:N_2+1+W_2);
IMcut = IMcut.*(IMcut>=0);
%%%%%%%%%%%%%%%%%%%%%%% Done %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Calculating residuals of the fit %%%%%%%%%%%%%%%%%%%%
RES=ims-P(N_2+1-H_2:N_2+1+H_2,N_2+1-W_2:N_2+1+W_2);
%%%%%%%%%% Done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% Calculating the speed distribution %%%%%%%%%%%%%%%%%%%%%%%
if W_2<=H_2
   Rmax=W_2;
else
   Rmax=H_2;
end;
spd=zeros(fix(Rmax/step));
B = load('matrT.txt', '-ASCII');
for R=step:step:(20-step)
   aux=Ci(1,1)*B(1,1)*exp(2*log(R)-R^2/sigma^2);
    for m=2:fix(R/sigma+11)
       aux=aux+Ci(m,1)*B(1,m)*R^2*exp((m-1)^2-R^2/sigma^2+((m-1)^2)*log(R^2/(sigma^2*((m-1)^2))));
    end;
   for k=2:fix(R/sigma+11)
      for m=1:(fix(sqrt((R+20)^2/sigma^2-(k-1)^2))+1)
            aux=aux+Ci(m,k)*B(k,m)*R^2*exp((k-1)^2+(m-1)^2-R^2/sigma^2+((k-1)^2+(m-1)^2)*log(R^2/(sigma^2*((k-1)^2+(m-1)^2))));
     end;
   end;
   if aux>0
      spd(R/step)=aux;
   end;
end;
for R=20:step:(Rmax-20)
   aux=0;
   for k=1:fix(R/sigma+11)
      for m=(fix(sqrt((R-20)^2/sigma^2-(k-1)^2))+2):(fix(sqrt((R+20)^2/sigma^2-(k-1)^2))+1)
         aux=aux+Ci(m,k)*B(k,m)*R^2*exp((k-1)^2+(m-1)^2-R^2/sigma^2+((k-1)^2+(m-1)^2)*log(R^2/(sigma^2*((k-1)^2+(m-1)^2))));
     end;
   end;
   if aux>0
      spd(R/step)=aux;
   end;
end;
for R=(Rmax-20+step):step:Rmax
   aux=0;
   for k=1:Rmax/sigma+1
      for m=(fix(sqrt((R-20)^2/sigma^2-(k-1)^2))+2):(fix(sqrt(Rmax^2/sigma^2-(k-1)^2))+1)
         aux=aux+Ci(m,k)*B(k,m)*R^2*exp((k-1)^2+(m-1)^2-R^2/sigma^2+((k-1)^2+(m-1)^2)*log(R^2/(sigma^2*((k-1)^2+(m-1)^2))));
     end;
   end;
   if aux>0
      spd(R/step)=aux;
   end;
end;
%spd=spd/max(spd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Calculating anisotropy parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
astep=0.5;
Np=fix(80/astep);
X=[ones(Np,1) zeros(Np,3)];
beta=zeros(fix(Rmax/step1),4);
for n=1:Np   
   X(n,2)=1.5*cos(pi*(n*astep+9)/180)^2-0.5;
   X(n,3)=(35*cos(pi*(n*astep+9)/180)^4-30*cos(pi*(n*astep+9)/180)^2+3)/8;
   X(n,4)=(231*cos(pi*(n*astep+9)/180)^6-315*cos(pi*(n*astep+9)/180)^4+105*cos(pi*(n*astep+9)/180)^2-5)/16;
end;
baux=inv(X'*X)*X';
for R=step1:step1:Rmax
   Y=zeros(Np,1);
   for n=1:Np
      aux=0;
      x=R*sin(pi*(n*astep+9)/180);
      z=R*cos(pi*(n*astep+9)/180);
      maxM=min(NBF,fix(z/sigma)+20);
      minM=max(1,fix(z/sigma)-20);
      maxK=min(NBF,fix(x/sigma)+20);
      minK=max(1,fix(x/sigma)-20);
      for m=minM:maxM
          if m==1
              mlnm=0;
          else
              mlnm=2*(m-1)^2*log(m-1);
          end;

         for k=minK:maxK
             if k==1
                 klnk=0;
             else 
                 klnk=2*(k-1)^2*log(k-1);
             end;
            aux=aux+Ci(m,k)*exp((k-1)^2+(m-1)^2-R^2/sigma^2+2*(k-1)^2*log(x/sigma)+2*(m-1)^2*log(z/sigma)-klnk-mlnm);
         end;
      end;
      Y(n)=aux;
   end;
   baux1=baux*Y;
%   fnm = ['Y', int2str(fix(R*2)) ,'.dat']
%   save(fnm,'Y', '-ASCII');
   beta(fix(R/step1),1)=R;
   beta(fix(R/step1),2)=baux1(2)/baux1(1);
   beta(fix(R/step1),3)=baux1(3)/baux1(1);
   beta(fix(R/step1),4)=baux1(4)/baux1(1);
   pause(0.1)
%   R
end;
%%%%%%%%%%%%%%%%%%%%%%%%%% Done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Plotting results %%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,2,1);
imagesc(IMcut);
title('Reconstructed image');
axis square
colorbar
subplot(2,2,2);
plot(spd);
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
%%%%%%%%%%%%%%%%%%%% Done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Saving output files %%%%%%%%%%%%%%%%%%%%%%%%
[pathstr,name,ext,versn] = fileparts(fname);
fname=fullfile(pathstr,[name '_3D' ext],versn);
fanis=fullfile(pathstr,[name '_anis.dat'],versn);
fspd=fullfile(pathstr,[name '_spd.dat'],versn);
f=fopen(fname,'w');
%%%%%%%% Image in ASCII-2 format %%%%%%%%%%%%%
for i=1:H
    for j=1:W
        fprintf(f,'%d,  %d  %2.5f\n',j,i,IMcut(i,j));
    end;
end;
fclose(f); 

%%%%%%%%%%% Anisotropy parameters %%%%%%%%%%%
save(fanis,'beta', '-ASCII');

%%%%%%%%%%%%%% Speed distribution %%%%%%%%%%%%%%%%
m=(1:length(spd(:,1)))';
aux=[m*step spd(:,1)];
save(fspd,'aux', '-ASCII');
%%%%%%%%%%%%%%%%% Done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
return


