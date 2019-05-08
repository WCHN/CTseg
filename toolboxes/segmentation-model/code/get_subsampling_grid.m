function [subsmp,grd] = get_subsampling_grid(dm,vs,samp)
% FORMAT [subsmp,grd] = get_subsampling_grid(dm,vs,samp)
% dm     - Full lattice dimensions
% vs     - Full lattice voxel size
% samp   - Sampling distance in mm
% subsmp - Subsampling info structure, with fields:
%          * sk  - Rounded sampling distance
%          * sk4 - Reshaped version: first 3 dimensions are singletons
%          * MT  - Inverse mapping (sub -> full) as an affine matrix
%          * dm  - Dimensions of subsampled grid
% grd    - Subsampling grid structure, with fields:
%          * x0  - Sampled indices along 1st dim, as a 2D grid
%          * y0  - Sampled indices along 2nd dim, as a 2D grid
%          * z0  - Sampled indices along 3rd dimension
%          * o   - ...
%
% Subsample a lattice and return useful info to map between the full 
% lattice and the subsampled lattice.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

sk        = max([1 1 1],round(samp*[1 1 1]./vs)); 
if dm(3) == 1
    % Images are 2D
    sk(3) = 1;
end

% For multiplying and dividing displacements to map from the subsampled 
% voxel indices and the actual image voxel indices.
sk4 = reshape(sk,[1 1 1 3]);

% Mapping from indices of subsampled voxels to indices of voxels in image(s).
MT = [sk(1) 0     0     (1 - sk(1)); ...
      0     sk(2) 0     (1 - sk(2)); ...
      0     0     sk(3) (1 - sk(3)); ...
      0     0     0     1];

[x0,y0,o] = ndgrid(1:sk(1):dm(1),1:sk(2):dm(2),1);
z0        = 1:sk(3):dm(3);

% Grid struct
grd.x0 = x0;
grd.y0 = y0;
grd.z0 = z0;    
grd.o  = o;

% Info struct
subsmp.sk  = sk;
subsmp.sk4 = sk4;
subsmp.MT  = MT;
subsmp.dm  = [size(grd.x0) length(grd.z0)]; % Sub-sampled dimensions
%==========================================================================