close all;
clear all;
raw = load('u:\users\scott\sections\circle.asc');
H=max(raw(:,1))-min(raw(:,1))+1;
W=max(raw(:,2))-min(raw(:,2))+1;
ims = reshape(raw(:,3),H,W);
[W H]=size(ims);
ims = ims';
figure
imagesc(ims)
save 'u:\users\scott\project\iterative\dataset.dat' ims -ASCII
return