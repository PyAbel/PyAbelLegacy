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
      maxM=min(NBF,z/sigma+20);
      minM=max(1,z/sigma-20);
      maxK=min(NBF,x/sigma+20);
      minK=max(1,x/sigma-20);
      for m=minM:maxM
         for k=minK:maxK
            aux=aux+Ci(m,k)*exp((k-1)^2+(m-1)^2-R^2/sigma^2+2*(k-1)^2*log(x/((k-1)*sigma))+2*(m-1)^2*log(z/((m-1)*sigma)));
         end;
      end;
      Y(n)=aux;
   end;
   baux1=baux*Y;
   beta(fix(R/step1),1)=R;
   beta(fix(R/step1),2)=baux(2)/baux(1);
   beta(fix(R/step1),3)=baux(3)/baux(1);
   beta(fix(R/step1),4)=baux(4)/baux(1);
end;

