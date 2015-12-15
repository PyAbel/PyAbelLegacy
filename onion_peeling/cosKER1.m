%function cosKER1(Im)

% Clear variables
clear all;
% Close all open windows
close all;

% Load the datafile into variable Im
%filename=['C:\Users\Eric Wells\Documents\Data\Research\Projects\Pulse Shaping\GA VMI SLM 2015\August 2015 D3p control\Best Images\run10_D3p_best_500V_10s_10shots.mod'];
filename=['run10_D3p_best_500V_10s_10shots.dat'];
Im=load(filename);
tic;
Im=augie_image_processing_function(Im);
dim=size(Im);
% halfheight=dim(1)/2; % number of rows (y-coordinate)
% halfwidth=dim(2)/2; % number of columns (x-coordinate)
radius = min(dim)/2; % Radius, based on smaller of height/width

angle=[0:360];
angle=angle/180*pi;
volt=0.5;
charge=1;
thres = 0.01;  % resolution in theta
rres = 0.05; % radial resolution (in eV)
k=88.0;  %depends on spectrometer

% Define center 
centerx=radius;
centery=radius;

logIm=Im+1;
logIm=log(logIm);


% % colorbar;
center = [centerx,centery];
maxradius = radius;
minradius=0;
param = 2*((volt*charge/k)^2);

sx=size(Im,2);  %get xsize
sy=size(Im,1);  % get ysize
[ydim,xdim] = size(Im);
rdim = floor(param * (sx^2+sy^2)/4 / rres);  %radial size
thdim = floor(2 / thres);
result = zeros(thdim,rdim);
num = zeros(thdim,rdim);

X=(1:sx)-center(1);
Y=(1:sy)-center(2);

theta=[0:359];
theta=theta/180*pi; %change to radians

R=minradius+1:maxradius;
%KER = 2*((volt*charge/k)^2).*(R.*R); % make axis for KER vs. cos(theta)

for i=1:length(theta)
   [xprime(:,i),yprime(:,i)]=pol2cart(theta(i)+pi/2,R);
end

ang_distribution = interp2(X,Y,Im,xprime,yprime,'*linear');
ang_distribution=[ang_distribution(:,1:180),ang_distribution(:,181:end)];

% % Plot the figure
%logang_distribution=ang_distribution+1;
%logang_distribution=log(logang_distribution);
Tang_distribution = ang_distribution';
%Tlogang_distribution = logang_distribution';
%zmax=max(max(Tang_distribution));
%zlmax = max(max(Tlogang_distribution));

figure(2);
imagesc(Tang_distribution);
shading interp;
grid off;
axis tight;
axis square;
colormap(jet(128));
grid off;
xlabel('radius');
ylabel('theta');
hold on;
% % colorbar;

%save .ver file
newfilename=[filename(1:end-3),'ver'];
save(newfilename,'ang_distribution','-ASCII');

%make cos vs. KER

zzx = radius; %center x
zzy = radius; % center y
[xx,yy] = meshgrid((-size(Im,2)/2):(size(Im,2)/2 - 1),(-size(Im,1)/2 + 1):(size(Im,1)/2));
[Th,R] = cart2pol(xx,-yy);
for y=1:ydim
    for x=1:xdim
        rcoord =((x-zzx)^2+(y-zzy)^2);
        thcoord = (y-zzy)/sqrt(rcoord);
        rcoord = param*rcoord;
        if thcoord < 0
            thcoord = thcoord + 2;
        end
        rcoord = fix(rcoord / rres);
        thcoord = fix(thcoord / thres);
        rcoord = max(rcoord,1);
        thcoord = max(thcoord,1);
        %result(thcoord,rcoord) = Im(y,x);
        result(thcoord,rcoord) = result(thcoord,rcoord) + Im(y,x);
        num(thcoord,rcoord) = num(thcoord,rcoord)+1;
    end
end
num = max(num,1);
result = result ./ num;

figure(3);
imagesc(result);
% shading interp;
grid off;
axis tight;
axis square;
colormap(jet(128));
grid off;
xlabel('KER');
ylabel('cos(theta)');
hold on;
toc;