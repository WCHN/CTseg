function p_out = spm_CTseg_warp(direction, p_in, p_y, p_atlas, prefix, odir, interp)
% Warp an atlas using CTseg deformations. 
% FORMAT p_out = spm_CTseg_warp(direction, p_in, p_y, p_atlas, prefix, odir, interp)
%
% OBS: The atlas needs to be in alignment with the SPM atlas!
%
% ARGS:
% --------------
% direction (str):
%   'atlas': Warp image to atlas space.
%   'image': Warp atlas to image space.
% p_in (str): Path to CT image, or its native space segmentation.
% p_y (str): Path to deformation (obtained from fitting CTseg to p_in).
% p_atlas (str): Path to atlas.
% prefix (str): Prefix of warped images. Default is 'w'.
% odir (str): Output directory. Default is same as input.
% interp (int): Interpolation order. Default is 1 (trilinear).
%
% RETURNS:
% --------------
% p_warped (str): Path to warped image/atlas.
%
%_______________________________________________________________________
if nargin < 5, prefix = 'w'; end
if nargin < 6, odir   = ''; end
if nargin < 7, interp = 1; end

% get parameters
% image
n_img = nifti(p_in);
% CTseg template
dir_ctseg = fileparts(mfilename('fullpath'));
p_mu      = fullfile(dir_ctseg,'mu_CTseg.nii');
n_mu      = nifti(p_mu);
M_mu      = n_mu.mat;
% atlas
n_atlas  = nifti(p_atlas);
M_atlas  = n_atlas.mat;
dm_atlas = n_atlas.dat.dim(1:3);
% affine from mni to mu
load(fullfile(dir_ctseg, 'Mmni.mat'), 'Mmni');

% define outputs
[pth,nam,ext] = fileparts(n_img.dat.fname);
if isempty(odir)
    odir = pth;
end
if ~exist(odir,'dir')
    mkdir(odir);
end
p_out = fullfile(odir,[prefix nam ext]);

% make affine warp
M_kl       = Mmni\M_atlas;
M          = M_mu*M_kl;
y          = spm_CTseg_util('affine',dm_atlas,M);
y          = reshape(y, [dm_atlas(1:3), 1, 3]);
p_y_mni2mu = fullfile(odir,'y_mni2mu.nii');
spm_CTseg_util('write_nii',p_y_mni2mu,single(y),M_atlas,'mu-to-mni affine');

% define warping
if strcmpi(direction,'atlas')
    defs.comp{1}.inv.comp{1}.def = {p_y};
    defs.comp{1}.inv.space = {p_mu};
    defs.comp{2}.def = {p_y_mni2mu};
    defs.out{1}.savedef.ofname = 'temp.nii';
    defs.out{1}.savedef.savedir.saveusr = {odir};
    defs.out{2}.pull.fnames = {p_in};
elseif strcmpi(direction,'image')
    defs.comp{1}.inv.comp{1}.def = {p_y_mni2mu};
    defs.comp{1}.inv.space = {p_atlas};
    defs.comp{2}.def = {p_y};
    defs.comp{3}.id.space = {p_in};
    defs.out{1}.savedef.ofname = 'temp.nii';
    defs.out{1}.savedef.savedir.saveusr = {odir};
    defs.out{2}.pull.fnames = {p_atlas};
else
    error(['Undefined direction: ' direction])
end
defs.out{2}.pull.savedir.saveusr = {odir};
defs.out{2}.pull.interp = interp;
defs.out{2}.pull.mask = 1;
defs.out{2}.pull.fwhm = [0 0 0];
defs.out{2}.pull.prefix = prefix;

% do warping
matlabbatch{1}.spm.util.defs = defs;
spm_jobman('run',matlabbatch);

% clean up
spm_unlink(fullfile(odir, 'y_temp.nii'));
spm_unlink(p_y_mni2mu);
%==========================================================================