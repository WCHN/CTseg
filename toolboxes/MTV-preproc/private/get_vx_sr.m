function vx = get_vx_sr(Nii)
% Get smallest voxel size and return as size of super-resolved images
% voxels.
%
% Nii - [1 x C cell array] Array with observed data
% vx  - Super-resolved images voxel size
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

mnvx = Inf;
C    = numel(Nii);
for c=1:C
    N = numel(Nii{c});
    for n=1:N
        mat = Nii{c}(n).mat; 
        vx  = sqrt(sum(mat(1:3,1:3).^2));  
        if min(vx) < mnvx
            mnvx = min(vx);
        end
    end
end
vx = round(mnvx,2)*ones(1,3);
%==========================================================================