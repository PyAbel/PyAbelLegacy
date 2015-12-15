close all;
clear all;

R1=1;
step=1;
Psize=fix(256/step);
P=zeros(Psize,2);

for n=1:Psize
   P(n,1)=R1+(n-1)*step;
   P(n,2)=(P(n,1)^2)*(2000*(28/3*exp(-((P(n,1)-10)^2)/4)+6*exp(-((P(n,1)-15)^2)/4)+10/3*exp(-((P(n,1)-20)^2)/4))+...
      200*(2*exp(-((P(n,1)-70)^2)/4)+4/3*exp(-((P(n,1)-85)^2)/4)+4/3*exp(-((P(n,1)-100)^2)/4))+...
      50*(8/3*exp(-((P(n,1)-145)^2)/4)+2*exp(-((P(n,1)-150)^2)/4)+2*exp(-((P(n,1)-155)^2)/4))+...
      40*exp(-((P(n,1)-45)^2)/3600));
end;
P(:,2)=P(:,2)/max(P(:,2));
save('u:\users\scott\project\tests\exactspd_fa.dat','P', '-ASCII');
return

