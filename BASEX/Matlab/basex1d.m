close all;
clear all;

fname=input('Input the name of the data file: ','s');
q=input('Choose the regularization parameter:');

Vsyms=input('Symmetrize the image with respect to vertical axis? Y/N :  ','s');
Vsymv=0;Vsymv=(Vsyms=='Y')|(Vsyms=='y');

Hsyms=input('Symmetrize the image with respect to horizontal axis? Y/N :  ','s');
Hsymv=0;Hsymv=(Hsyms=='Y')|(Hsyms=='y');

ims = load(fname);
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


N = 1001;
N_2=fix(N/2);
W_2=fix(W/2);
H_2=fix(H/2);
im = zeros(N,N);
im(N_2+1-H_2:N_2+H_2+1,N_2+1-W_2:N_2+W_2+1)=ims;

M = load('basis1000pr_1.txt');
Mc = load('basis1000_1.txt');
NBF=size(M,2);

E = zeros(NBF,NBF);
for n=1:NBF
   E(n,n)=q;
end;

Ci=(im)*M*inv(M'*M+E);
P = Ci*M';
IM=Ci*Mc';

IMobrez=IM(N_2+1-H_2:N_2+1+H_2,N_2+1-W_2:N_2+1+W_2);
IMobrez = IMobrez.*(IMobrez>=0);
Pobrez = P(N_2+1-H_2:N_2+1+H_2,N_2+1-W_2:N_2+1+W_2);
Pobrez = Pobrez.*(Pobrez>=0);
RES=ims-Pobrez;

figure
subplot(2,2,1);
imagesc(IMobrez);
title('Reconstructed image');
axis square
colorbar
subplot(2,2,2);
imagesc(ims);
title('Original image (projection)');
axis square
colorbar
subplot(2,2,3);
imagesc(Pobrez);
title('Projection of the reconstructed image');
axis square
colorbar
subplot(2,2,4);
imagesc(RES);
title('Residue');
axis square
colorbar

[pathstr,name,ext,versn] = fileparts(fname);
fname=fullfile(pathstr,[name '_3D' ext],versn);
fname2=fullfile(pathstr,[name '_proj' ext],versn);

save(fname,'IMobrez', '-ASCII');
save(fname2,'Pobrez', '-ASCII');
return

