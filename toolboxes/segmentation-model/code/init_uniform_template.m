function model = init_uniform_template(dat,opt)
% FORMAT model = init_uniform_template(dat,opt)
% dat   - Subjects data structure
% opt   - Options structure
% model - Model structure
%
% Initialise the template as "uniform"
% -> Same value everywhere
% -> Same probability for all classes
% Create template files on disk
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

K         = opt.template.K;
dir_model = opt.dir_model;
S0        = numel(dat);
vx        = opt.template.vs;
keep_neck = opt.template.keep_neck;
do_resiz  = opt.template.resize;
    
% Collect orientation matrices and dimensions from all input images
%--------------------------------------------------------------------------
mats = [];
dms  = [];   
for s=1:S0
    if isfield(dat{s}.modality{1},'channel')
        C = numel(dat{s}.modality{1}.channel);
        for c=1:C
            mat  = dat{s}.modality{1}.channel{c}.nii.mat;           
            mats = cat(3,mats,mat);

            dm  = dat{s}.modality{1}.channel{c}.nii.dat.dim;                          
            if numel(dm)==2, dm(3) = 1; end
            dms = cat(1,dms,dm);
        end
    else
        mat  = dat{s}.modality{1}.nii.mat;           
        mats = cat(3,mats,mat);

        dm  = dat{s}.modality{1}.nii.dat.dim;                          
        if numel(dm)==2, dm(3) = 1; end
        dms = cat(1,dms,dm);
    end
end

% Compute average orientation matrix and dimensions
%--------------------------------------------------------------------------
[M,d] = spm_misc('compute_avg_mat',mats,dms);

if keep_neck
    % If we model spine, remove a little bit of the bottom of the template
    V.dim = d;
    V.mat = M;
    bb    = world_bb(V);

    mrg = 25; % mm to remove from bottom (we are in world space, and vx = 1)
    if d(3) == 1
        bb(1,2) = bb(1,2) + mrg;
    else
        bb(1,3) = bb(1,3) + mrg;
    end
    
    mn = bb(1,:);
    mx = bb(2,:);
    M  = spm_matrix([mn 0 0 0 ones(1,3)])*spm_matrix([-1 -1 -1]);    
    d  = ceil(M \ [mx 1]' - 0.1)';
    d  = d(1:3);
end

% Write template to disk
%--------------------------------------------------------------------------
pth_template  = fullfile(dir_model,'template.nii');   
[pth,nam,ext] = fileparts(pth_template);

vols = cell(K,1);
for k=1:K    
    
    img = zeros(d,'single');
    
    vols{k} = fullfile(pth,[nam num2str(k) ext]);
    spm_misc('create_nii',vols{k},img,M,[spm_type('float32') spm_platform('bigend')],'template');
end
clear img

matlabbatch                       = cell(1,1);
matlabbatch{1}.spm.util.cat.vols  = vols;
matlabbatch{1}.spm.util.cat.name  = pth_template;    
matlabbatch{1}.spm.util.cat.dtype = 0;
spm_jobman('run',matlabbatch);

for k=1:K
    delete(vols{k});
end

% Set template
model              = struct;
model.template.nii = nifti(pth_template);

% Change voxel size and/or crop template
if ~isempty(vx) && d(3) == 1
    dm0  = d;
    mat0 = M;
    vx0  = sqrt(sum(mat(1:3,1:3).^2));
       
    % New orientation matrix and dimensions
    ds    = vx0./vx;
    D     = diag([ds 1]);
    mat   = mat0/D;
    dm    = floor(D(1:3,1:3)*dm0')';    
    dm(3) = 1;

    % Update template
    model.template.nii = spm_misc('create_nii',pth_template,zeros([dm K],'single'),mat,[spm_type('float32') spm_platform('bigend')],'template');        
    d                  = dm;
elseif (opt.template.resize || ~isempty(vx)) && d(3) > 1    
    model     = resize_template(model,do_resiz,vx,keep_neck);
    d         = model.template.nii.dat.dim;
end

if ~isempty(opt.lesion.hemi)
    for k=1:K          
        % Assumes there are two lesion classes, one for each hemisphere
        if k == opt.lesion.hemi{1}
            img(1:floor(d(1)/2 - 0.05*d(1)),:,:)  = log(1e-2);
            img = spm_imbasics('smooth_img',img,12);
        elseif k == opt.lesion.hemi{2}
            img(ceil(d(1)/2 + 0.05*d(1)):end,:,:) = log(1e-2);
            img = spm_imbasics('smooth_img',img,12);
        end

        model.template.nii.dat(:,:,:,k) = img;
    end
end

% Get values to be used for FOV voxels when warping 
model = init_template_bg(model,opt);

% Init template objval
model.template.objval.post  = 0;
model.template.objval.likel = 0;
model.template.objval.pr    = 0;

alpha           = ones(1,K)*opt.prop.reg;
PropPrior.alpha = alpha;
PropPrior.norm  = 0;
model.PropPrior = PropPrior;

% Save PropPrior
fname = fullfile(opt.dir_model,'PropPrior.mat');
save(fname,'PropPrior');
%==========================================================================