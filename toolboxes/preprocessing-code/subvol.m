function VO = subvol(V,bb,prefix,deg,constrain_mx)
% Extract a subvolume
% FORMAT VO = subvol(V,bb,prefix)
% V      - SPM volume object
% bb     - bounding box
% prefix - file prefix (if empty -> overwrites)
% VO     - resized image
%
% Example:
%     V = spm_vol(spm_select(1,'image'));
%     subvol(V,[32 64 ; 1 64 ; 1 48]');
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin < 3, prefix       = 'sv_'; end
if nargin < 4, deg          = 0;     end
if nargin < 5, constrain_mx = true;     end

bb      = round(bb);
bb      = sort(bb);
bb(1,:) = max(bb(1,:),[1 1 1]);
if constrain_mx
    bb(2,:) = min(bb(2,:),V.dim(1:3));
end

VO            = V;
[pth,nam,ext] = fileparts(V.fname);
VO.fname      = fullfile(pth,[prefix nam ext]);
VO.dim(1:3)   = diff(bb)+1;
VO.mat        = V.mat*spm_matrix((bb(1,:)-1));

VO = spm_create_vol(VO);
for z=1:VO.dim(3)
    M   = V.mat\VO.mat*spm_matrix([0 0 z]);
    img = spm_slice_vol(V,M,VO.dim(1:2),deg);
    VO  = spm_write_plane(VO,img,z);
end
%==========================================================================