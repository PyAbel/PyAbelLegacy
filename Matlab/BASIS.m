Rm = 501;
N=1001;
NBF=500;
I=(1:N)';
R2  = ((I-Rm).^2);
R  = I-Rm;
M = zeros(N,NBF);
Mc =zeros(N,NBF);
Mc(:,1)=exp(-R2);
M(:,1)=2*exp(-R2);
logR2=log(R2);
for nn=1:(NBF-1)
    
    n=nn*nn;
    pot=n*log(R2)-R2+n-n*log(n);
    pot(Rm)=double((nn==0));
    Mc(:,nn+1)=exp(pot);
    Mc(Rm,nn+1)=double((nn==0));
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
        
    M(:,nn+1)=2*Mc(:,nn+1)+aux;
    nn
end;

save 'd:\temp\basis1000_1.bst' Mc -ASCII;
save 'd:\temp\basis1000pr_1.bst' M -ASCII;
return
