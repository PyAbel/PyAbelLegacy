close all;
clear all;
theta1=30;
theta=theta1*pi/180;
raw=load('u:\users\scott\sum02c.asc');
H=max(raw(:,1))-min(raw(:,1))+1;
W=max(raw(:,2))-min(raw(:,2))+1;
h=ceil(H/2);
w=ceil(W/2);
l=length(raw);
[M N]=size(raw);
sect=zeros(M,N);
sect(:,1)=raw(:,1)-min(raw(:,1))+1;
sect(:,2)=raw(:,2)-min(raw(:,2))+1;
oth=sect;
image=sect;

for i=1:l
   if sect(i,2)<=-H/(W*tan(theta))*sect(i,1)+h*(1+cot(theta));
      sect(i,3)=raw(i,3);
   end;
   if sect(i,1)<h;
      sect(i,3)=0;
   end;
   if sect(i,2)>w;
     sect(i,3)=0;
   end;
   if oth(i,2)>=-H/(W*tan(theta))*oth(i,1)+h*(1+cot(theta));
      oth(i,3)=raw(i,3);
   end;
   if oth(i,1)>h;
      oth(i,3)=0;
   end;
   if oth(i,2)<w;
     oth(i,3)=0;
   end;
   if image(i,2)>=(W*tan(-theta)/H)*image(i,1)+W/2-(H/2)*tan(-theta);
      image(i,3)=raw(i,3);
   end;
   if image(i,1)<h;
      image(i,3)=0;
   end;
   if image(i,2)>w;
     image(i,3)=0;
   end;
end;

save 'u:\users\scott\project\fourier-hankel\upper.dat' sect -ASCII;
save 'u:\users\scott\project\fourier-hankel\lower.dat' oth -ASCII;
save 'u:\users\scott\project\fourier-hankel\test.dat' image -ASCII;

H=max(sect(:,1))-min(sect(:,1))+1;
W=max(sect(:,2))-min(sect(:,2))+1;
ims = reshape(sect(:,3),H,W);
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
   end;
end;
vel=vel/max(vel);

H=max(oth(:,1))-min(oth(:,1))+1;
W=max(oth(:,2))-min(oth(:,2))+1;
ims = reshape(oth(:,3),H,W);
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
   velo=zeros(W_2,1);
   for m=1:W_2-1
      for k=1:m
         z=m*cos(pi*k/(2*m));
         x=m*sin(pi*k/(2*m));
         velo(m)=vel(m)+(z-fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+2,W_2+fix(x)+2)*x+...
         (1-z+fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+2)*x+...
         (z-fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+2,W_2+fix(x)+1)*x+...
         (1-z+fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+1)*x;
      end;
   end;
else
   velo=zeros(H_2,1);
   for m=1:H_2-1
      for k=1:m
         z=m*cos(pi*k/(2*m));
         x=m*sin(pi*k/(2*m));
         velo(m)=vel(m)+(z-fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+2,W_2+fix(x)+2)*x+...
         (1-z+fix(z))*(x-fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+2)*x+...
         (z-fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+2,W_2+fix(x)+1)*x+...
         (1-z+fix(z))*(1-x+fix(x))*IMobrez(H_2+fix(z)+1,W_2+fix(x)+1)*x;
      end;
   end;
end;
vel=vel/max(vel)
velo=velo/max(velo);
figure
plot(vel);
hold on;
plot(velo,'r');
title('Speed distribution in Upper Quadrant(blue) & Lower Quadrant(red)');

m=(1:length(vel))';
aux=[m vel(:)];
save('u:\users\scott\upper.dat','aux', '-ASCII');

m=(1:length(velo))';
aux=[m velo(:)];
save('u:\users\scott\lower.dat','aux', '-ASCII');

return
