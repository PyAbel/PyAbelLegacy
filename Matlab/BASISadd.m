Rm = 501;
N=1001;
NBF=250;
I=(1:N)';
R2  = ((I-Rm).^2);
R  = I-Rm;
M = zeros(N,NBF);
Mc =zeros(N,NBF);
sigma=4;
R2_sigma=R2./sigma;
logR2=log(R2);
logsigma=log(sigma);
Mc(:,1)=exp(-R2_sigma);
M(:,1)=2*sqrt(sigma)*exp(-R2_sigma);
for nn=240:(NBF-1)
    
    n=nn*nn;
    pot=n*log(R2)-(sigma^(-1))*R2+n-n*log(n*sigma);
    pot(Rm)=double((nn==0));
    Mc(:,nn+1)=exp(pot);
    Mc(Rm,nn+1)=double((nn==0));
    aux=zeros(N,1);
    fn=n-n*log(sigma*n);
    for k=1:n-1
        lnCnk=sum(log((k+1):n))-sum(log(1:(n-k)));
        m=1:k;
        lnP=sum(log((2*m-1)/2));
        u=-R2_sigma+fn +(n-k)*logR2+lnCnk+(k+0.5)*logsigma+lnP;
        aux=aux + 2*exp(u);
    end;
    aux(Rm)=0;

    m=1:n;
    lnP=sum(log((2*m-1)/2));
    u=-R2_sigma+fn+(n+0.5)*logsigma+lnP;
    aux=aux + 2*exp(u);
        
    M(:,nn+1)=2*sqrt(sigma)*Mc(:,nn+1)+aux;
    nn
end;

save 'e:\add_2.bst' Mc -ASCII;
save 'e:\addpr_2.bst' M -ASCII;
return
