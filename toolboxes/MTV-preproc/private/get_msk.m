function msk = get_msk(f,mu)
% Get mask defining what voxels to regard as missing data
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

msk = isfinite(f) & f ~= 0;
if nargin > 1
    msk = msk & isfinite(mu);
end
msk = msk(:);
%==========================================================================