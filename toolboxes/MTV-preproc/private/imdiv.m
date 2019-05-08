%==========================================================================
function D = imdiv(X,vs,type)  
% Computes the divergence of an image (with voxel size)
% FORMAT div = spm_imbasics('div',img,(vs),(type))
% img   - Image in "gradient space" (Nx * Nz * Nz * ... * Ndim * Ntype)
% vs    - Voxel size [1]
% type  - Finite difference type '+'/'-'/['+-']/'-+'
% div   - Divergence (Nx * Nz * Nz * ...)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% -------------------------------------------------------------------------
% Default parameters
if nargin < 2, vs = ones(1,'like',X); end
if nargin < 3, type = '+-'; end

dim = size(X);
if numel(type) > 1
    dimt = numel(dim);
    dim = dim(1:end-1);
else
    dimt = numel(dim) + 1;
end
dimd = numel(dim);
dim  = dim(1:end-1);
vs   = padarray(vs(:)', [0 max(0,3-numel(dim))], 'replicate', 'post');

% -------------------------------------------------------------------------
% Sum finite differences
D = zeros(dim,'like', X);
for i=1:numel(size(D))
    for t=1:numel(type)
        Xit = select_slices(X, [dimd dimt], {i t});
        switch type(t)
            case '+'
                D = D + cat(i,      -select_slices(Xit, i, 1), ...
                               -diff(select_slices(Xit, i, 1:(dim(i)-1)), 1, i), ...
                                     select_slices(Xit, i, dim(i)-1)) ./ vs(i);
            case '-'
                D = D + cat(i,      -select_slices(Xit, i, 2), ...
                               -diff(select_slices(Xit, i, 2:dim(i)), 1, i), ...
                                     select_slices(Xit, i, dim(i))) ./ vs(i);
        end
    end
end
D = D/sqrt(numel(type));
%==========================================================================

% function D = imdiv(imx,imy,imz,vx)  
% % Compute image divergence
% % _______________________________________________________________________
% %  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
% 
% if nargin < 4, vx = ones([1 3],'like',imx); end
% 
% Du = cat(2, -imx(:,1,:), -diff(imx(:,1:end-1,:),1,2), imx(:,end-1,:)); 
% Dv = cat(1, -imy(1,:,:), -diff(imy(1:end-1,:,:),1,1), imy(end-1,:,:));
% Dw = cat(3, -imz(:,:,1), -diff(imz(:,:,1:end-1),1,3), imz(:,:,end-1));
% 
% D = Du./vx(1) + Dv./vx(2) + Dw./vx(3);
% %==========================================================================