close all;
clear all;
B2=zeros(501,501);
B2(1,1)=2;
kk(1)=0;
for i=2:250000
    kk(i)=kk(i-1)+i*log(i)-(i-1)*log(i-1)-log(i);
end;
i
mm(1)=log(2);
for j=2:500000
    mm(j)=mm(j-1)+j*log(j)-(j-1)*log(j-1)-log(j-0.5);
end;
j
for m=1:500
    B2(1,m+1)=1/(m^2+0.5);
end;
for k=1:500
    i=1:k^2;
    B2(k+1,1)=exp(sum(log(i))-sum(log(i-0.5)))/(k^2+0.5);
end;

for m=1:500
    for k=1:500
        B2(k+1,m+1)=exp(mm(k^2+m^2)-kk(k^2)-mm(m^2))/(m^2+k^2+0.5);
    end;
    m
end;
save('e:\matrT2.asc','B2', '-ASCII');
return