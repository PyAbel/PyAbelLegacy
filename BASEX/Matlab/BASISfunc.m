close all;
clear all;
nn=10;
N=501;
Rm=fix(N/2)+1;
I=(1:N)';
R2  = 0.01*((I-Rm).^2);
R  = 0.1*(I-Rm);
M = zeros(N,3);
M(:,1)=R;
logR2=log(R2);
    
    n=nn*nn;
    pot=n*log(R2)-R2+n-n*log(n);
    pot(Rm)=double((nn==0));
    M(:,2)=exp(pot);
    M(Rm,2)=double((nn==0));
    M(:,2)=fix(10^100*M(:,2))/10^100;
    aux=zeros(N,1);
    fn=n-n*log(n);
    for k=1:n-1
        lnCnk=sum(log((k+1):n))-sum(log(1:(n-k)));
        m=1:k;
        lnP=sum(log((2*m-1)/2));
        u=-R2+fn +(n-k)*logR2+lnCnk+lnP;
        aux=aux + 2*exp(u);
    end;
    aux(Rm)=0;

    m=1:n;
    lnP=sum(log((2*m-1)/2));
    u=-R2+fn+lnP;
    aux=aux + 2*exp(u);  
    M(:,3)=2*M(:,2)+aux;
    M(:,3)=fix(10^100*M(:,3))/10^100;
figure
plot(M(:,2));
figure
plot(M(:,3));
save 'e:\bf10.dat' M -ASCII;
return
