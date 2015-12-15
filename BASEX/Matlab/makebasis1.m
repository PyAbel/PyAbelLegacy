clear all;
M = load('basis1000pr_1.bst');
Mc = load('basis1000_1.bst');
sigma=1;
Rm = 501;
N=1001;
I=(1:N)';
R2  = ((I-Rm).^2);
Mc(:,1)=exp(-R2/sigma^2);
M(:,1)=2*sigma*exp(-R2/sigma^2);

f = fopen('e:\basisset1.bin','wb');
fwrite(f,Mc','float');
fclose(f);

NBF=size(Mc,2);
NBFc=NBF;

RANK = NBF-0;
[u s v]=svd(M);
ds = diag(s);
dsi = [1./ds(1:RANK); zeros(NBF-RANK,1)];
si = [diag(dsi);zeros(N-NBF,NBF)];
invM=(u*si*v')';

RANK = NBF-0;
[u s v]=svd(Mc);
ds = diag(s);
dsi = [1./ds(1:RANK); zeros(NBF-RANK,1)];
si = [diag(dsi);zeros(N-NBF,NBF)];
invMc=(u*si*v')';

f = fopen('e:\invbasis_pr1.bin','wb');
fwrite(f,invM','float');
fclose(f);

E = zeros(NBF,NBF);
for n=1:NBF
   E(n,n)=5;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
invM=M*inv(M'*M+E);

f = fopen('e:\Rinvbasis_pr1.bin','wb');
fwrite(f,invM,'float');
fclose(f);

f = fopen('e:\invbasis1.bin','wb');
fwrite(f,invMc','float');
fclose(f);

f = fopen('e:\basissize1.txt','w');
fprintf(f,'%d %d',size(Mc));
fclose(f);

return
