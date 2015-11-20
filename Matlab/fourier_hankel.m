close all;
clear all;

fname=input('Input the name of the data file: ','s');

Vsyms=input('Symmetrize the image with respect to vertical axis? Y/N :  ','s');
Vsymv=0;Vsymv=(Vsyms=='Y')|(Vsyms=='y');

Hsyms=input('Symmetrize the image with respect to horizontal axis? Y/N :  ','s');
Hsymv=0;Hsymv=(Hsyms=='Y')|(Hsyms=='y');

raw = load(fname);
H=max(raw(:,1))-min(raw(:,1))+1;
W=max(raw(:,2))-min(raw(:,2))+1;
ims = reshape(raw(:,3),H,W);
[W H]=size(ims);
ims = ims';

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

W_2=fix(W/2)+1;
H_2=fix(H/2)+1;

temp=fft(ims');
fc=exp(-2*sqrt(-1)*pi*(W_2-1)/W*(0:W-1))';
for n=1:H
   temp(:,n)=temp(:,n).*fc;
end;
F=real(temp');
FF=F(:,1:W_2);
    for i=1:W_2
        for j=1:W_2
            %q=0.5*(j-1)/W_2;
            q=(j-1)/W;
            B(j,i)=q*besselj(0,2*pi*q*(i-1));
        end;
    end
Hkl=FF*B;
IMobrez=[Hkl(:,W_2-1:-1:1) Hkl];
IMobrez = IMobrez.*(IMobrez>=0);
if W_2<=H_2
   vel=zeros(W_2,1);
   for m=1:W_2-5
      for k=1:m
         z=m*cos(pi*k/(2*m));
         x=m*sin(pi*k/(2*m));
         vel(m)=vel(m)+(z-fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+1)*x+...
         (1-z+fix(z))*(x-fix(x))*IMobrez(H_2+fix(z),W_2+fix(x)+1)*x+...
         (z-fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x))*x+...
         (1-z+fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z),W_2+fix(x))*x;
      end;
      X=[zeros(90,1) ones(90,1)];
      Y=zeros(90,1);
      for k=1:90
          X(k,1)=1.5*(cos(pi*k/180))^2-0.5;
          z=m*cos(pi*k/180);
          x=m*sin(pi*k/180);
          Y(k)=(z-fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+1)+...
         (1-z+fix(z))*(x-fix(x))*IMobrez(H_2+fix(z),W_2+fix(x)+1)+...
         (z-fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x))+...
         (1-z+fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z),W_2+fix(x));
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
         vel(m)=vel(m)+(z-fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+1)*x+...
         (1-z+fix(z))*(x-fix(x))*IMobrez(H_2+fix(z),W_2+fix(x)+1)*x+...
         (z-fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x))*x+...
         (1-z+fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z),W_2+fix(x))*x;
      end;
         X=[zeros(90,1) ones(90,1)];
      Y=zeros(90,1);
      for k=1:90
          X(k,1)=1.5*(cos(pi*k/180))^2-0.5;
          z=m*cos(pi*k/180);
          x=m*sin(pi*k/180);
          Y(k)=(z-fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+1)+...
         (1-z+fix(z))*(x-fix(x))*IMobrez(H_2+fix(z),W_2+fix(x)+1)+...
         (z-fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x))+...
         (1-z+fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z),W_2+fix(x));
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
subplot(2,1,1);
imagesc(IMobrez);
title('Reconstructed image');
axis square
colorbar
subplot(2,1,2);
plot(vel);
hold on;
plot(beta,'r');
title('Speed distribution(blue) & anisotropy(red)');
axis square

[pathstr,name,ext,versn] = fileparts(fname);
fname=fullfile(pathstr,[name '_3D' ext],versn);
fvel=fullfile(pathstr,[name '_spd.dat'],versn);
fanis=fullfile(pathstr,[name '_anis.dat'],versn);
f=fopen(fname,'w');
for i=1:H
    for j=1:W
        fprintf(f,'%d,  %d  %2.5f\n',j,i,IMobrez(i,j));
    end;
end;
fclose(f); 
 
 m=(1:length(vel))';
 aux=[m vel(:)];
 save(fvel,'aux', '-ASCII');

 m=(1:length(beta))';
 aux=[m beta(:)];
 save(fanis,'aux', '-ASCII');
return

