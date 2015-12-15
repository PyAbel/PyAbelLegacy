close all;
clear all;
Rm = 501;
N=1001;
NBF=500;
M = zeros(N,NBF);
Mc =zeros(N,NBF);
Mc(Rm,1)=1;
M(Rm,:)=2;
for n=2:NBF
    Mc(Rm-n+1,n)=1;
    Mc(Rm+n-1,n)=1;
    for i=1:(n-1)
        M(Rm-i,n)=2*(sqrt(n^2-i^2)-sqrt((n-1)^2-i^2));
        M(Rm+i,n)=M(Rm-i,n);
    end;
    
end;

save 'd:\step1000_1.bst' Mc -ASCII;
save 'd:\step1000pr_1.bst' M -ASCII;
return
