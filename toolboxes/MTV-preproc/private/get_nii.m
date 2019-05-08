function img = get_nii(nii,z,dim)
% Read data
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 2, z   = 0; end
if nargin < 3, dim = 3; end

if iscell(nii)
    nii = nii{1};
    N   = numel(nii);
    img = cell(1,N);
    for n=1:N
        if z == 0
            img{n} = single(nii(n).dat());
        else
            img{n} = single(select_slices(nii(n).dat,dim,z));    
        end
    end
else
    if z == 0
        img = single(nii.dat());
    else
        img = single(select_slices(nii.dat,dim,z));    
    end   
end
%==========================================================================