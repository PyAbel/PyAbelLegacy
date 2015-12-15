close all;
clear all;

raw=load('u:\users\scott\project\iterative\pofe.dat');
h=length(raw(:,1));
P=zeros(h,2);

P(:,1)=sqrt(raw(:,1))-0.5;
P(:,2)=raw(:,2).*P(:,1);

P(:,2)=P(:,2)/max(P(:,2));
plot(P(:,1),P(:,2));

save('u:\users\scott\project\tests\spd.dat','P', '-ASCII');

