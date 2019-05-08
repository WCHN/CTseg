function [dat,opt] = init_reg(dat,opt)
% FORMAT [dat,opt] = init_reg(dat,opt)
% dat   - Subjects data structure
% opt   - Options structure
%
% Init registration related variables:
% * opt.reg.B:              Lie basis
% * opt.reg.aff_reg_ICO:    Affine precision matrix
% * dat.reg.r:              Affine coordinates in Lie algebra (all zero)
% * dat.reg.v:              Create velocity nifti file on disk
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

S0   = numel(dat);
dm   = obs_info(dat{1});
type = opt.reg.aff_type;
if dm(3) > 1
    flat = '3d';
else
    flat = '2d';
end
dt = [spm_type('float32') spm_platform('bigend')];

% Affine basis function
B                   = spm_misc('affine_basis', type, flat);
Nr                  = size(B,3);
opt.reg.B           = B;
opt.reg.aff_reg_ICO = opt.reg.aff_reg*eye(Nr);

% Adjust ICO so that we do not regularise rigid part
if dm(3) > 1
    % 3D: Six translation and rotation parameters
    M = diag([zeros(1,6) ones(1,Nr - 6)]);
else
    % 2D: Three translation and rotation parameters
    M = diag([zeros(1,3) ones(1,Nr - 3)]);
end
opt.reg.aff_reg_ICO = M.*opt.reg.aff_reg_ICO;

% Loop over all subjects
for s=1:S0
    % Subject image info
    [dm,mat,vs] = obs_info(dat{s});
    
    % If subsampling (samp > 0), adjust dimensions of velocities
    % accordingly
    subsamp = get_subsampling_grid(dm,vs,opt.seg.samp);
    dm      = subsamp.dm;
    clear subsamp
    
    % Create velocities
    v = zeros([dm(1:3),3],'single');
    if opt.template.do
        pth_v        = fullfile(opt.dir_vel,['v' num2str(s) '.nii']);
        spm_misc('create_nii',pth_v,v,mat,dt,'vel');
        dat{s}.reg.v = nifti(pth_v);
    else
        if S0 == 1
            dat{s}.reg.v = v;
        else
            pth_v        = fullfile(dat{s}.dir.vel,'v.nii');
            spm_misc('create_nii',pth_v,v,mat,dt,'vel');
            dat{s}.reg.v = nifti(pth_v);
        end
    end
    
    % Create affine parameters
    dat{s}.reg.r = zeros([Nr,1]);
end
%==========================================================================