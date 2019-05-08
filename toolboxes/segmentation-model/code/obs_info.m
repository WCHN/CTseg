function [dm,mat,vs,C,nam,V,fnames,chn_names] = obs_info(dat)
% FORMAT [dm,mat,vs,C,nam,V,fnames,chn_names] = obs_info(dat)
% dat       - Subject's data structure (one subject)
% dm        - Observed image dimension
% mat       - Observed image orientation matrix
% vs        - Observed image voxel size
% C         - Number of channels
% nam       - Modality name
% V         - spm_vol handle of observed images
% fnames    - Filenames of observed images
% chn_names - Channel names
%
% Return useful information about the observed image.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if isfield(dat.modality{1},'channel')    
    C         = numel(dat.modality{1}.channel);
    fnames    = cell(1,C);
    chn_names = cell(1,C);
    for c=1:C
       fnames{c}    = dat.modality{1}.channel{c}.nii.dat.fname;
       chn_names{c} = dat.modality{1}.channel{c}.name;
    end    
    dm  = dat.modality{1}.channel{1}.nii.dat.dim;
    mat = dat.modality{1}.channel{1}.nii.mat;    
else
    fnames    = {dat.modality{1}.nii.dat.fname};
    C         = 1;
    chn_names = {dat.modality{1}.name};
    dm        = dat.modality{1}.nii.dat.dim;
    mat       = dat.modality{1}.nii.mat;        
end
nam = cell(1,C);
for c=1:C
   [~,nam{c}] = fileparts(fnames{c});
end
fnames  = char(fnames);
V       = spm_vol(fnames);
vs      = spm_misc('vxsize',mat);

if numel(dm)==2, dm(3) = 1; end

if any(dm == 1) && find(dm == 1) ~= 3
    error('find(dm == 1) ~= 3')
end
%==========================================================================