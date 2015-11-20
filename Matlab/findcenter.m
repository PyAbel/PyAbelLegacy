close all;
clear all;
raw = load('d:\temp\n11f213s.asc');
H=max(raw(:,1))-min(raw(:,1))+1;
W=max(raw(:,2))-min(raw(:,2))+1;
ims = reshape(raw(:,3),H,W);
[W H]=size(ims);
ims = ims';
Im=zeros(3*H-2,3*W-2);
Iminv=zeros(3*H-2,3*W-2);
Im(H:2*H-1,W:2*W-1)=ims;
dw=10;
Xi=fix((3*W-2)/2)-dw;
Zi=fix((3*H-2)/2)-dw;
xmin=0;
zmin=0;
mindel=10^100;
for dx=1:2*dw
	for dz=1:2*dw
%dx=50;
%dz=50;
      for z=1:H
         for x=1:W
            Iminv(2*(Zi+dz)-H-z+1,2*(Xi+dx)-W-x+1)=ims(z,x);
         end;
      end;
      delta=(Iminv-Im).^2;
      del=sum((sum(delta))');
      if del<mindel
         mindel=del;
         xmin=Xi+dx-W+1;
         zmin=Zi+dz-H+1;
      end;
    end;
 end;
 xmin
 zmin
return;
      
