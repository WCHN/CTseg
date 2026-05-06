function create_atlas(vox, pth_out, v_mult, spm_fov, fwhm)
% Create a CTseg atlas aligned to SPM tissue probability (TPM) space.
%
% The CTseg template (mu_CTseg.nii) is registered directly to SPM's
% TPM.nii as a categorical (tissue-probability) input, using Multi-Brain.
% Because both source and target are tissue probability fields, no
% intensity modelling is needed. The resulting deformation is then used
% to warp mu_CTseg.nii into SPM/TPM space at the requested voxel size.
%
% Output filename defaults to mu_CTseg_spm15.nii (1.5 mm) or
% mu_CTseg_spm10.nii (1.0 mm).
%
% Usage:
%   create_atlas                                 % 1.5 mm isotropic
%   create_atlas(1.0)                            % 1.0 mm isotropic
%   create_atlas(1.5, 'my_atlas.nii')            % custom output path
%   create_atlas(1.5, '', 4, true, 4)            % v_settings multiplier, SPM FOV, smoothing FWHM
%
% Args:
%   vox     (double): Output voxel size in mm (default 1.5).
%   pth_out (char):   Output file path (default: mu_CTseg_spm{10,15}.nii
%                     in CTseg root).
%   v_mult  (double): Multiplier for v_settings (default 4).
%   spm_fov (bool):   Output FOV same as SPM atlas (default true).
%   fwhm    (double): FWHM in mm for Gaussian smoothing of TPM channels
%                     before registration (default 0 = no smoothing).
%
% Requires:
%   - SPM12 on MATLAB path (with Shoot, Longitudinal, MB toolboxes)

    if nargin < 1 || isempty(vox),     vox     = 1.5; end
    if nargin < 2,                     pth_out = ''; end
    if nargin < 3 || isempty(v_mult),  v_mult  = 4; end
    if nargin < 4 || isempty(spm_fov), spm_fov = true; end
    if nargin < 5 || isempty(fwhm),    fwhm    = 0; end

    spm('defaults', 'fmri');
    spm_jobman('initcfg');

    % Paths
    dir_ctseg = fileparts(which('spm_CTseg'));
    pth_mu    = fullfile(dir_ctseg, 'models', 'mu_CTseg.nii');
    spm_dir   = spm('Dir');
    pth_tpm   = fullfile(spm_dir, 'tpm', 'TPM.nii');

    if ~exist(pth_mu, 'file')
        error('CTseg template not found: %s', pth_mu);
    end
    if ~exist(pth_tpm, 'file')
        error('SPM TPM not found: %s', pth_tpm);
    end

    % Output directory (temporary)
    odir = fullfile(dir_ctseg, 'temp_spm_template');
    if ~exist(odir, 'dir'), mkdir(odir); end

    % -----------------------------------------------------------------
    % Categorical registration: TPM tissue priors -> CTseg template
    % -----------------------------------------------------------------
    fprintf('Running MB registration (categorical: TPM -> mu_CTseg)...\n');

    % Load TPM and extract first K channels (drop implicit background)
    Nii_tpm  = nifti(pth_tpm);
    K        = size(nifti(pth_mu).dat, 4);
    tpm_data = single(Nii_tpm.dat(:,:,:,1:K));
    fprintf('  TPM: using %d of %d channels (K=%d)\n', ...
        K, size(Nii_tpm.dat, 4), K);

    % Optionally smooth each TPM channel
    if fwhm > 0
        fprintf('  Smoothing TPM channels with FWHM=%gmm...\n', fwhm);
        pth_tmp = fullfile(odir, 'smooth_tmp.nii');
        for k = 1:K
            spm_CTseg_util('write_nii', pth_tmp, tpm_data(:,:,:,k), ...
                Nii_tpm.mat, 'tmp', 'float32');
            spm_smooth(pth_tmp, pth_tmp, [fwhm fwhm fwhm]);
            Nii_tmp = nifti(pth_tmp);
            tpm_data(:,:,:,k) = single(Nii_tmp.dat(:,:,:));
        end
        delete(pth_tmp);
    end

    % Write K-channel TPM to temp 4D file
    pth_tpm_input = fullfile(odir, 'tpm_input.nii');
    dm              = size(tpm_data);
    Nii_cat         = nifti;
    Nii_cat.dat     = file_array(pth_tpm_input, dm, 'float32', 0);
    Nii_cat.mat     = Nii_tpm.mat;
    Nii_cat.mat0    = Nii_tpm.mat;
    Nii_cat.descrip = 'TPM input for categorical MB registration';
    create(Nii_cat);
    Nii_cat.dat(:,:,:,:) = tpm_data;
    clear tpm_data

    % Build config and call spm_mb_init/spm_mb_fit directly
    mb_cfg              = struct;
    mb_cfg.mu.exist     = {pth_mu};
    mb_cfg.onam         = 'spm_template';
    mb_cfg.odir         = {odir};
    mb_cfg.v_settings   = [0.00001 0 0.4 0.1 0.4] * v_mult;
    mb_cfg.aff          = 'Aff(3)';
    mb_cfg.del_settings = Inf;
    mb_cfg.accel        = 0.8;
    mb_cfg.min_dim      = 8;
    mb_cfg.tol          = 0.001;
    mb_cfg.sampdens     = 2;
    mb_cfg.save         = true;
    mb_cfg.nworker      = 0;
    mb_cfg.cat          = {{pth_tpm_input}};
    mb_cfg.gmm          = [];

    [dat_mb, sett_mb] = spm_mb_init(mb_cfg);
    [dat_mb, ~]       = spm_mb_fit(dat_mb, sett_mb);

    pth_y = dat_mb(1).psi.dat.fname;
    fprintf('Deformation field: %s\n', pth_y);

    % Copy deformation to CTseg directory
    pth_y_spm = fullfile(dir_ctseg, 'y_CTseg_spm.nii');
    copyfile(pth_y, pth_y_spm);
    fprintf('Saved deformation: %s\n', pth_y_spm);

    % Warp mu_CTseg.nii (4D, K-1 channels) to SPM space
    fprintf('Warping CTseg template to SPM space...\n');
    Nii_mu = nifti(pth_mu);
    K1     = Nii_mu.dat.dim(4);  % number of channels (K-1)

    % Reference space: SPM TPM (optionally) at requested voxel size
    if spm_fov
        pth_ref_base = pth_tpm;
    else
        pth_ref_base = pth_mu;
    end
    if ~isempty(vox)
        pth_ref = fullfile(odir, 'ref_space.nii');
        create_ref_space(pth_ref_base, vox, pth_ref);
    else
        pth_ref = pth_ref_base;
    end

    warped_channels = cell(1, K1);
    Mout = [];
    for k = 1:K1
        % Write single channel to temp file
        pth_chan = fullfile(odir, sprintf('mu_chan%d.nii', k));
        dat_chan = single(Nii_mu.dat(:,:,:,k));
        spm_CTseg_util('write_nii', pth_chan, dat_chan, Nii_mu.mat, ...
            sprintf('CTseg template channel %d', k), 'float32');

        matlabbatch = {};
        matlabbatch{1}.spm.util.defs.comp{1}.def      = {pth_y_spm};
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

        delete(pth_chan);
        delete(pth_wchan);
    end

    % Assemble 4D warped template
    dm_out = size(warped_channels{1});
    mu_spm = zeros([dm_out K1], 'single');
    for k = 1:K1
        mu_spm(:,:,:,k) = warped_channels{k};
    end

    % Output path: derive from vox if not explicitly set
    if ~isempty(pth_out)
        pth_mu_spm = pth_out;
        if ~endsWith(pth_mu_spm, '.nii')
            pth_mu_spm = [pth_mu_spm '.nii'];
        end
    else
        if abs(vox - 1.0) < 1e-6
            pth_mu_spm = fullfile(dir_ctseg, 'mu_CTseg_spm10.nii');
        elseif abs(vox - 1.5) < 1e-6
            pth_mu_spm = fullfile(dir_ctseg, 'mu_CTseg_spm15.nii');
        else
            pth_mu_spm = fullfile(dir_ctseg, sprintf('mu_CTseg_spm%gmm.nii', vox));
        end
    end
    Nii_out         = nifti;
    Nii_out.dat     = file_array(pth_mu_spm, [dm_out K1], 'float32', 0);
    Nii_out.mat     = Mout;
    Nii_out.mat0    = Mout;
    Nii_out.descrip = 'CTseg template in SPM space';
    create(Nii_out);
    Nii_out.dat(:,:,:,:) = mu_spm;

    fprintf('Saved warped template: %s\n', pth_mu_spm);

    % Cleanup
    if exist(pth_y_spm, 'file'), delete(pth_y_spm); end
    if exist(odir, 'dir'), rmdir(odir, 's'); end

    fprintf('Done.\n');
end


function create_ref_space(pth, vox, pth_out)
% Create a reference NIfTI at specified voxel size, covering the reference FOV.
    if numel(vox) == 1, vox = vox * [1 1 1]; end
    V_ref  = spm_vol(pth);
    vx_ref = sqrt(sum(V_ref(1).mat(1:3,1:3).^2));
    D      = diag([vx_ref ./ vox, 1]);
    mat    = V_ref(1).mat / D;
    dim    = floor(D(1:3,1:3) * V_ref(1).dim(1:3)')';
    Nii         = nifti;
    Nii.dat     = file_array(pth_out, dim, 'float32', 0);
    Nii.mat     = mat;
    Nii.mat0    = mat;
    Nii.descrip = sprintf('Reference space at %gmm', vox(1));
    create(Nii);
    Nii.dat(:,:,:) = zeros(dim, 'single');
end
