clear all;
astep=0.5;
Np=fix(80/astep);
X=[ones(Np,1) zeros(Np,3)];
for n=1:Np   
   X(n,2)=1.5*cos(pi*(n*astep+9)/180)^2-0.5;
   X(n,3)=(35*cos(pi*(n*astep+9)/180)^4-30*cos(pi*(n*astep+9)/180)^2+3)/8;
   X(n,4)=(231*cos(pi*(n*astep+9)/180)^6-315*cos(pi*(n*astep+9)/180)^4+105*cos(pi*(n*astep+9)/180)^2-5)/16;
end;
matrX=inv(X'*X)*X';

f = fopen('matrX.bin','wb');
fwrite(f,matrX','float');
fclose(f);
return

