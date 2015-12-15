close all;
clear all;

raw = load('e:\dimer_scan\ion_images\processing\235sum12.txt');
W = length(raw(1,:));
H = length(raw(:,1));
f=fopen('e:\dimer_scan\ion_images\processing\235sum12.asc','w');
for i=1:H
    for j=1:W
        fprintf(f,'%d,  %d  %2.5f\n',j,i,raw(i,j));
    end;
end;
fclose(f);
return