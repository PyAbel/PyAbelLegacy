close all;
clear all;
raw=load('u:\users\scott\project\tests\blanket noise\fourier-abel\inverse.asc');
H=max(raw(:,1))-min(raw(:,1))+1;
W=max(raw(:,2))-min(raw(:,2))+1;
ims = reshape(raw(:,3),H,W);
[W H]=size(ims);
IM = ims';
W_2=fix(W/2);
H_2=fix(H/2);
if (W/2-W_2)==0
   IMobrez=[IM(:,1:W_2) IM(:,W_2) IM(:,W_2+1:W)];
else
   IMobrez=IM;
end;

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
plot(vel);
hold on;
plot(beta,'r');
title('Speed distribution(blue) & anisotropy(red)');

m=(1:length(vel))';
aux=[m vel(:)];
save('u:\users\scott\project\tests\blanket noise\fourier-abel\spdinverse.dat','aux', '-ASCII');

m=(1:length(beta))';
aux=[m beta(:)];
save('u:\users\scott\project\tests\center shift-2left\fourier-abel\anis2left.dat','aux', '-ASCII');

return
