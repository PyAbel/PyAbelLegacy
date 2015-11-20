NBF = 300;
T = zeros(NBF,NBF);
for m=0:NBF-1
   T(1,m+1)=1/(m^2+0.5);
end;
for k=0:NBF-1
   T(k+1,1)=10^(-10)*QUADL(@tau0,-1,1,10^(-10),[],k);
end;

for m = 1:NBF-1
    for k = 1:NBF-1
         T(k+1,m+1)=10^(-10)*QUADL(@tau,-1,1,10^(-10),[],k,m);    
     end;
     m
 end;
 save 'matrT.asc' T -ASCII;
 return
    
    