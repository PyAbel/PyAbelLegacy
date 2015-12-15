function abel_arr=augie_image_processing_function(work_im)
%%% Image processing with cartesian Abel-inversion
%%% Abel-inversion algorithm from: Rev. Sci. Instrum. 67, 2257 (1996).
%%% Original Author: Boris Bergues
%%% Current Revision Author: Chris Rallis
%%% last changes: 06.18.10
%%% info about this implementation:  Rallis et al., Rev. Sci. Instrum. 85, 113105 (2014) 

% filename=['D:\Physics Research\GA-VMI\added mod files\','addedC+_46.mod'];
% work_im=load(filename);

% sym_lr=0;
% sym_ud=0;

% make dimensions even
dim=size(work_im);
if mod(dim(1),2)~=0,work_im=[work_im;(1:dim(2))*0];end
dim=size(work_im);
if mod(dim(2),2)~=0,work_im=[work_im,(1:dim(1))'*0];end

% raw_im=work_im;
% disp(['filename: ',filename])

% Zero negative values prior
work_im=max(work_im,0);

% Symmetrize Left and Right
work_im=0.5*(work_im+flipdim(work_im,2));
sym_lr=1;

% Symmetrize Up and Down
work_im=0.5*(work_im+flipdim(work_im,1));
sym_ud=1;

% find dimensions:
dim=size(work_im);
ymax=dim(1); % number of rows (y-coordinate)

if sym_ud == 1
    work_im = work_im(1:fix(ymax/2),:);
end

% update dimensions:
dim=size(work_im);
ymax=dim(1); % number of rows (y-coordinate)
xmax=dim(2); % number of columns (x-coordinate)

% update center:

xcl=fix(xmax/2);
xcr=fix(xmax/2)+1;
if sym_ud~=1,yc=fix(ymax/2);end
if sym_ud==1,yc=fix(ymax);end

% Initialize abel routine
[val1,val2]=init_abel(xcl,yc);
abel_arr=work_im*0;
rest_arr=work_im;
vect=(1:ymax)*0;

% process right half of the image
for col_index = xmax:-1:xcr
    idist = col_index-xcr;
    rest_row = rest_arr(:,col_index);
    normfac = 1/val1(idist+1,idist+1);
    for i = idist:-1:0
        ic = xcr + i;
        shadow=rest_row*normfac*val1(i+1,idist+1);
        rest_arr(:,ic)=rest_arr(:,ic)-shadow;
        rest_arr(:,ic)=max(rest_arr(:,ic),0);
    end
    for row_index =1:ymax
        jdist=abs(row_index-yc);
        vect(row_index)=val2(idist+1,jdist+1);
    end
    abel_arr(:,col_index)=normfac*rest_row.*vect';
end

% process left half of the image
if sym_lr ~=1
for col_index = 1:xcl
    idist = xcl-col_index;
    rest_row = rest_arr(:,col_index);
    normfac = 1/val1(idist+1,idist+1);
    for i = idist:-1:0
        ic = xcl-i;
        shadow=rest_row*normfac*val1(i+1,idist+1);
        rest_arr(:,ic)=rest_arr(:,ic)-shadow;
        rest_arr(:,ic)=max(rest_arr(:,ic),0);
    end
    for row_index = 1:ymax
        jdist=abs(row_index - yc);
        vect(row_index)=val2(idist+1,jdist+1);
    end
    abel_arr(:,col_index+1)=normfac*rest_row.*vect';
end
end

if sym_ud ==1,abel_arr=[abel_arr;flipdim(abel_arr,1)];end
if sym_lr ==1,abel_arr=abel_arr+flipdim(abel_arr,2);end
imagesc(abel_arr);