function Nii = RunPreproc(Nii,opt)

addpath('/home/mbrud/dev/mbrud/code/matlab/MTV-preproc')

if nargin < 2, opt = struct; end
% Do
if ~isfield(opt,'do'),          opt.do          = struct; end
if ~isfield(opt.do,'res_orig'), opt.do.res_orig = false; end
if ~isfield(opt.do,'real_mni'), opt.do.real_mni = true; end
if ~isfield(opt.do,'crop'),     opt.do.crop     = true; end
if ~isfield(opt.do,'coreg'),    opt.do.coreg    = true; end
if ~isfield(opt.do,'denoise'),  opt.do.denoise  = true; end
if ~isfield(opt.do,'reslice'),  opt.do.reslice  = true; end
% Output directory
if ~isfield(opt,'dir_out'),     opt.dir_out     = 'output'; end
% Reslice options
if ~isfield(opt,'reslice'),     opt.reslice     = struct; end
if ~isfield(opt.reslice,'ref'), opt.reslice.ref = 1; end

if ~iscell(Nii)
    Nii = {Nii};
end

ResliceRef = 2; 

% Copy so to not overwrite originals
Nii = make_copies(Nii,opt.dir_out);

if opt.do.res_orig
    % Reset origin (important for CT)
    Nii = reset_origin(Nii,VoxSize);
end

if opt.do.real_mni
    % Realing to MNI space
    Nii = realign2mni(Nii);
end

if opt.do.crop
    % Remove uneccesary data
    Nii = crop(Nii);
end

if opt.do.coreg
    % Coreg
    Nii = coreg(Nii);
end

if opt.do.denoise
    % Denoise
    Nii = denoise(Nii);

    % Coreg (one more time after denoising)
    Nii = coreg(Nii);
end

if opt.do.reslice
    % Make same size
    Nii = reslice(Nii,opt.reslice.ref);
end
%==========================================================================