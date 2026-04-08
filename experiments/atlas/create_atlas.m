function create_atlas(mri_type, vox, pth_out, inu_mult, v_mult)
% Create a CTseg atlas nonlinearly aligned to MNI space.
% Runs Multi-Brain registration on MNI image(s) with uninformative
% priors, then warps mu_CTseg.nii to MNI space using the resulting
% deformation.
%
% Creates (in CTseg root directory, or at pth_out if specified):
%   mu_CTseg_spm.nii  — CTseg atlas warped to MNI space (4D, K-1 classes)
%
% Usage:
%   create_atlas                          % avg152T1, TPM voxel size (~1.5mm)
%   create_atlas('single')                % single_subj_T1, TPM voxel size
%   create_atlas('mni')                   % MNI152_T1_1mm, TPM voxel size
%       Source: https://github.com/Jfortin1/MNITemplate/blob/master/inst/extdata/MNI152_T1_1mm.nii.gz
%   create_atlas('icbm_asym')             % ICBM152 2009c asym T1+T2+PD, TPM voxel size
%   create_atlas('icbm_sym')              % ICBM152 2009c sym T1+T2+PD, TPM voxel size
%   create_atlas('icbm_asym', 1)          % ICBM152 2009c asym T1+T2+PD, 1mm isotropic
%   create_atlas('icbm_sym', 1.5)         % ICBM152 2009c sym T1+T2+PD, 1.5mm isotropic
%   create_atlas('single', 1)             % single_subj_T1, 1mm isotropic
%   create_atlas('avg', 1)                % avg152T1, 1mm isotropic
%   create_atlas('mni', 1.5, p, 1, 4)     % custom output path (p), inu_reg mult (im) and v_settings mult (vm)
%
% Optional args:
%   pth_out  (char):   Output file path (default: mu_CTseg_spm.nii in CTseg dir)
%   inu_mult (double): Multiplier for inu_reg (default 1 -> inu_reg=1e4)
%   v_mult   (double): Multiplier for v_settings (default 4)
%
% Requires:
%   - SPM12 on MATLAB path (with Shoot, Longitudinal, MB toolboxes)

    if nargin < 1, mri_type = 'avg'; end
    if nargin < 2, vox      = []; end
    if nargin < 3, pth_out  = ''; end
    if nargin < 4, inu_mult = 1; end
    if nargin < 5, v_mult   = 4; end

    spm('defaults', 'fmri');
    spm_jobman('initcfg');

    % Paths
    dir_ctseg = fileparts(which('spm_CTseg'));
    pth_mu = fullfile(dir_ctseg, 'mu_CTseg.nii');

    % Select MNI image(s)
    spm_dir = spm('Dir');
    exp_dir = fileparts(fileparts(mfilename('fullpath')));  % experiments/
    switch lower(mri_type)
        case 'avg'
            pth_imgs = {fullfile(spm_dir, 'canonical', 'avg152T1.nii')};
        case 'single'
            pth_imgs = {fullfile(spm_dir, 'canonical', 'single_subj_T1.nii')};
        case 'mni'
            pth_imgs = {fullfile(spm_dir, 'canonical', 'MNI152_T1_1mm.nii')};
        case 'icbm_asym'
            icbm_dir = fullfile(exp_dir, 'mni_icbm152_nlin_asym_09c_nifti', ...
                'mni_icbm152_nlin_asym_09c');
            pth_imgs = {
                fullfile(icbm_dir, 'mni_icbm152_t1_tal_nlin_asym_09c.nii')
                fullfile(icbm_dir, 'mni_icbm152_t2_tal_nlin_asym_09c.nii')
                fullfile(icbm_dir, 'mni_icbm152_pd_tal_nlin_asym_09c.nii')
            };
        case 'icbm_sym'
            icbm_dir = fullfile(exp_dir, 'mni_icbm152_nlin_sym_09c_nifti', ...
                'mni_icbm152_nlin_sym_09c');
            pth_imgs = {
                fullfile(icbm_dir, 'mni_icbm152_t1_tal_nlin_sym_09c.nii')
                fullfile(icbm_dir, 'mni_icbm152_t2_tal_nlin_sym_09c.nii')
                fullfile(icbm_dir, 'mni_icbm152_pd_tal_nlin_sym_09c.nii')
            };
        otherwise
            error('Unknown mri_type: %s. Use ''avg'', ''single'', ''mni'', ''icbm_asym'', or ''icbm_sym''.', mri_type);
    end
    for i = 1:numel(pth_imgs)
        if ~exist(pth_imgs{i}, 'file')
            error('Image not found: %s', pth_imgs{i});
        end
    end
    fprintf('Using %d MRI channel(s):\n', numel(pth_imgs));
    for i = 1:numel(pth_imgs)
        fprintf('  Channel %d: %s\n', i, pth_imgs{i});
    end

    % Output directory (temporary)
    odir = fullfile(dir_ctseg, 'temp_spm_template');
    if ~exist(odir, 'dir'), mkdir(odir); end

    % Run MB to register MNI image(s) to CTseg template
    fprintf('Running MB registration...\n');
    run              = struct;
    run.mu.exist     = {pth_mu};
    run.onam         = 'spm_template';
    run.odir         = {odir};
    run.v_settings   = [0.00001 0 0.4 0.1 0.4] * v_mult;
    run.aff          = 'Aff(3)';
    run.gmm.pr.file          = {''};  % empty = uninformative prior (computed from image)
    run.gmm.pr.hyperpriors   = [];
    for c = 1:numel(pth_imgs)
        run.gmm.chan(c).images      = pth_imgs(c);
        run.gmm.chan(c).modality    = 1;  % MRI
        run.gmm.chan(c).inu.inu_reg = 1e4 * inu_mult;
    end

    out        = struct;
    out.result = {fullfile(odir, ['mb_fit_' run.onam '.mat'])};

    jobs{1}.spm.tools.mb.run = run;
    jobs{2}.spm.tools.mb.out = out;
    spm_jobman('run', jobs);

    % Load results to get deformation field
    res = load(out.result{1});
    dat = res.dat;
    pth_y = dat(1).psi.dat.fname;
    fprintf('Deformation field: %s\n', pth_y);

    % Copy deformation to CTseg directory
    pth_y_spm = fullfile(dir_ctseg, 'y_CTseg_spm.nii');
    copyfile(pth_y, pth_y_spm);
    fprintf('Saved deformation: %s\n', pth_y_spm);

    % Warp mu_CTseg.nii (4D, K-1 channels) to MNI space
    fprintf('Warping CTseg template to MNI space...\n');
    Nii_mu = nifti(pth_mu);
    K1     = Nii_mu.dat.dim(4);  % number of channels (K-1)

    % Warp each channel using spm.util.defs (pull through forward deformation)
    warped_channels = cell(1, K1);
    Mout = [];
    for k = 1:K1
        % Write single channel to temp file
        pth_chan = fullfile(odir, sprintf('mu_chan%d.nii', k));
        dat_chan = single(Nii_mu.dat(:,:,:,k));
        spm_CTseg_util('write_nii', pth_chan, dat_chan, Nii_mu.mat, ...
            sprintf('CTseg template channel %d', k), 'float32');

        % Warp to MNI space
        pth_tpm = fullfile(spm_dir, 'tpm', 'TPM.nii');
        if ~isempty(vox)
            % Create reference at requested voxel size (derived from TPM bounding box)
            pth_ref = fullfile(odir, 'ref_space.nii');
            if ~exist(pth_ref, 'file')
                create_ref_space(pth_tpm, vox, pth_ref);
            end
        else
            pth_ref = pth_tpm;
        end
        matlabbatch = {};
        matlabbatch{1}.spm.util.defs.comp{1}.def     = {pth_y_spm};
        matlabbatch{1}.spm.util.defs.comp{2}.id.space = {pth_ref};
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {pth_chan};
        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {odir};
        matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
        matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = 'w';
        spm_jobman('run', matlabbatch);

        % Read warped channel
        pth_wchan = fullfile(odir, sprintf('wmu_chan%d.nii', k));
        Nii_w = nifti(pth_wchan);
        warped_channels{k} = single(Nii_w.dat(:,:,:));
        Mout = Nii_w.mat;

        % Cleanup temp channel files
        delete(pth_chan);
        delete(pth_wchan);
    end

    % Assemble 4D warped template
    dm_out = size(warped_channels{1});
    mu_spm = zeros([dm_out K1], 'single');
    for k = 1:K1
        mu_spm(:,:,:,k) = warped_channels{k};
    end

    % Write output template
    if ~isempty(pth_out)
        pth_mu_spm = pth_out;
        if ~endsWith(pth_mu_spm, '.nii')
            pth_mu_spm = [pth_mu_spm '.nii'];
        end
    else
        pth_mu_spm = fullfile(dir_ctseg, 'mu_CTseg_spm.nii');
    end
    Nii_out         = nifti;
    Nii_out.dat     = file_array(pth_mu_spm, [dm_out K1], 'float32', 0);
    Nii_out.mat     = Mout;
    Nii_out.mat0    = Mout;
    Nii_out.descrip = 'CTseg template in MNI space';
    create(Nii_out);
    Nii_out.dat(:,:,:,:) = mu_spm;

    fprintf('Saved warped template: %s\n', pth_mu_spm);

    % Cleanup
    if exist(pth_y_spm, 'file'), delete(pth_y_spm); end
    if exist(odir, 'dir'), rmdir(odir, 's'); end

    fprintf('Done.\n');
end


function create_ref_space(pth_tpm, vox, pth_out)
% Create a reference NIfTI at specified voxel size, covering the TPM FOV.
    if numel(vox) == 1, vox = vox * [1 1 1]; end
    V_tpm  = spm_vol(pth_tpm);
    vx_tpm = sqrt(sum(V_tpm(1).mat(1:3,1:3).^2));
    D      = diag([vx_tpm ./ vox, 1]);
    mat    = V_tpm(1).mat / D;
    dim    = floor(D(1:3,1:3) * V_tpm(1).dim(1:3)')';
    Nii         = nifti;
    Nii.dat     = file_array(pth_out, dim, 'float32', 0);
    Nii.mat     = mat;
    Nii.mat0    = mat;
    Nii.descrip = sprintf('Reference space at %gmm', vox(1));
    create(Nii);
    Nii.dat(:,:,:) = zeros(dim, 'single');
end
