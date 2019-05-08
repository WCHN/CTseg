function y = resize_def(y,dm,subsmp,bs)
% FORMAT y = resize_def(y,dm_s,subsmp,bs)
% y      - Input deformation
% dm     - Size of output deformation
% subsmp - Struct with sub-sampling info
% bs     - B-spline interpolation parameters
%
% Resize a deformation using b-spline interpolation.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 4, bs = [4 4 4 1 1 1]; end

% Mapping from voxels in input grid to voxels in output grid
iMT = inv(subsmp.MT);   

% Create interpolation grid
[x1,y1,z1] = ndgrid(1:dm(1),1:dm(2),1:dm(3));
x1         = x1*iMT(1,1) + iMT(1,4);
y1         = y1*iMT(2,2) + iMT(2,4);
z1         = z1*iMT(3,3) + iMT(3,4);
if dm(3) == 1
    % Data is 2D
    z1     = ones(size(z1));
end

% Subtract identity from deformation
% Because the deformation is y = y + id, and we want to interpolate y
% before id is added
dm0       = size(y);
id        = cell(3,1);
[id{1:3}] = ndgrid(single(1:dm0(1)),single(1:dm0(2)),single(1:dm0(3)));
id        = cat(4,id{:});
y         = y - id;
clear id

% Do interpolation
ny = zeros([dm 3],'single');
for i=1:3
    coeff       = spm_bsplinc(y(:,:,:,i),bs);
    ny(:,:,:,i) = spm_bsplins(coeff,x1,y1,z1,bs);            
end
clear x1 y1 z1   
ny(~isfinite(ny)) = 0;

% Rescale interpolated deformations with sub-sampling factor, then add back
% identity
id        = cell(3,1);
[id{1:3}] = ndgrid(single(1:dm(1)),single(1:dm(2)),single(1:dm(3)));
id        = cat(4,id{:});
ny        = bsxfun(@times,ny,subsmp.sk4);
y         = ny + id; 
%==========================================================================