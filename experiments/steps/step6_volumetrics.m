function step6_volumetrics(cfg)
% Step 6: Compute TBV and TIV from native-space tissue maps for all three
% methods (MR-SPM, SPM-CT, CTseg). Tissue maps are masked by an
% intracranial mask derived from the CTseg atlas (GM+WM+CSF > 0.5),
% warped to native space via the SPM-MR deformation field.
% Also collects runtimes. Saves to cfg.out_dir/volumetrics.mat.

    fprintf('=== Step 6: Volumetrics (TBV/TIV) ===\n');
    subjects = get_subjects(cfg);
    n = numel(subjects);

    % Create intracranial mask in template space from CTseg atlas
    icv_mask_template = create_icv_mask(cfg.pth_mu);
    fprintf('  ICV mask: %d / %d template voxels (%.1f%%)\n', ...
        sum(icv_mask_template(:)), numel(icv_mask_template), ...
        100 * sum(icv_mask_template(:)) / numel(icv_mask_template));

    tbv_mr    = nan(n, 1);  tiv_mr    = nan(n, 1);
    tbv_spm   = nan(n, 1);  tiv_spm   = nan(n, 1);
    tbv_ctseg = nan(n, 1);  tiv_ctseg = nan(n, 1);
    rt_spm_mr  = nan(n, 1);
    rt_spm_ct  = nan(n, 1);
    rt_ctseg   = nan(n, 1);
    valid = false(n, 1);

    for i = 1:n
        subj_dir = fullfile(cfg.data_dir, subjects(i).name);
        mr_files    = get_tissue_files(subj_dir, 'mr_spm');
        ct_files    = get_tissue_files(subj_dir, 'ct_spm');
        ctseg_files = get_tissue_files(subj_dir, 'ct_ctseg');

        % Require all three methods
        if ~exist(mr_files{1}, 'file') || ...
           ~exist(ct_files{1}, 'file') || ...
           ~exist(ctseg_files{1}, 'file')
            continue;
        end

        % Get SPM-MR deformation to warp template mask to native space
        [~, nam_mr] = fileparts(mr_files{1});
        nam_mr = regexprep(nam_mr, '^c[123]', '');  % c1pp_mr -> pp_mr
        pth_y = fullfile(subj_dir, ['y_' nam_mr '.nii']);
        if ~exist(pth_y, 'file')
            fprintf('  [%d/%d] %s — skipping (no MR deformation)\n', i, n, subjects(i).name);
            continue;
        end
        valid(i) = true;

        % Warp ICV mask to native space
        V_tissue = spm_vol(mr_files{1});
        icv_mask_native = warp_mask_to_native(icv_mask_template, cfg.pth_mu, pth_y, V_tissue);

        % All methods: compute TBV/TIV with ICV mask
        [tbv_mr(i),    tiv_mr(i)]    = compute_volumes(mr_files,    icv_mask_native);
        [tbv_spm(i),   tiv_spm(i)]   = compute_volumes(ct_files,    icv_mask_native);
        [tbv_ctseg(i), tiv_ctseg(i)] = compute_volumes(ctseg_files, icv_mask_native);

        % Runtimes
        rt_file = fullfile(subj_dir, 'runtime_spm_mr.mat');
        if exist(rt_file, 'file')
            R = load(rt_file); rt_spm_mr(i) = R.runtime;
        end
        rt_file = fullfile(subj_dir, 'runtime_spm_ct.mat');
        if exist(rt_file, 'file')
            R = load(rt_file); rt_spm_ct(i) = R.runtime;
        end
        vol_file = fullfile(subj_dir, 'vol_CTseg.mat');
        if exist(vol_file, 'file')
            S = load(vol_file);
            if isfield(S, 'runtime'), rt_ctseg(i) = S.runtime; end
        end

        fprintf('  [%d/%d] %s  MR: TBV=%.0f TIV=%.0f  SPM-CT: TBV=%.0f TIV=%.0f  CTseg: TBV=%.0f TIV=%.0f ml\n', ...
            i, n, subjects(i).name, tbv_mr(i), tiv_mr(i), tbv_spm(i), tiv_spm(i), tbv_ctseg(i), tiv_ctseg(i));
    end

    % Save
    volumetrics = struct();
    volumetrics.subject_names = {subjects(valid).name};
    volumetrics.tbv_mr        = tbv_mr(valid);
    volumetrics.tiv_mr        = tiv_mr(valid);
    volumetrics.tbv_spm       = tbv_spm(valid);
    volumetrics.tiv_spm       = tiv_spm(valid);
    volumetrics.tbv_ctseg     = tbv_ctseg(valid);
    volumetrics.tiv_ctseg     = tiv_ctseg(valid);
    volumetrics.rt_spm_mr     = rt_spm_mr(valid);
    volumetrics.rt_spm_ct     = rt_spm_ct(valid);
    volumetrics.rt_ctseg      = rt_ctseg(valid);

    % Print runtime summary
    fprintf('\n  Runtimes (segmentation only):\n');
    for method = {'spm_mr', 'spm_ct', 'ctseg'}
        rt = volumetrics.(['rt_' method{1}]);
        rt = rt(~isnan(rt));
        if ~isempty(rt)
            fprintf('    %-8s  mean=%.1fs  median=%.1fs  std=%.1fs  range=[%.1f, %.1f]s  (n=%d)\n', ...
                method{1}, mean(rt), median(rt), std(rt), min(rt), max(rt), numel(rt));
        end
    end

    save(fullfile(cfg.out_dir, 'volumetrics.mat'), 'volumetrics');
    fprintf('  Saved to %s\n', fullfile(cfg.out_dir, 'volumetrics.mat'));
end

function icv_mask = create_icv_mask(pth_mu)
% Create binary intracranial mask from CTseg atlas.
% The atlas stores K-1 log-probability channels. Softmax (matching Multi-Brain's
% softmax0/LSE0) with an implicit background class gives tissue probabilities.
% ICV = GM+WM+CSF (classes 1-3).
    Nii = nifti(pth_mu);
    mu = single(Nii.dat(:,:,:,:));  % dim x K-1 (log-space)
    dm = size(mu);

    % Softmax (Multi-Brain's softmax0: numerically stable, implicit background)
    mu_prob = softmax0(mu);
    % Append implicit background class
    mu_prob = cat(4, mu_prob, max(1 - sum(mu_prob, 4), 0));

    % Sum intracranial classes (GM + WM + CSF = channels 1-3)
    icv_prob = sum(mu_prob(:,:,:,1:3), 4);

    % Threshold
    icv_mask = icv_prob > 0.5;
    fprintf('  Created ICV mask from atlas (%dx%dx%d, threshold=0.5)\n', dm(1), dm(2), dm(3));
end

function P = softmax0(mu, ax)
% Safe softmax function (matches LSE0 from Multi-Brain).
% Handles implicit background class (mu=0).
    if nargin < 2, ax = 4; end
    mx  = max(mu, [], ax);
    E   = exp(bsxfun(@minus, mu, mx));
    den = sum(E, ax) + exp(-mx);
    P   = bsxfun(@rdivide, E, den);
end

function mask_native = warp_mask_to_native(icv_mask, pth_mu, pth_y, V_target)
% Warp template-space ICV mask to native space using SPM deformation.
% The SPM deformation y maps template voxels to native world coordinates.
% We use it to sample the template mask at each native voxel's corresponding
% template location.
%
% Strategy: for each native voxel, find its template-space coordinate by
% inverting the deformation lookup. Since the deformation y is defined on
% the template grid (mapping template→native), we need to find for each
% native voxel which template voxel maps closest to it.
%
% Simpler approach: use SPM's spm_diffeo to pull the mask through the
% inverse deformation. SPM stores the deformation as template→native,
% so we write the mask as a NIfTI in template space, then use
% spm_deformations to pull it to native space.

    dm_native = V_target.dim;

    % Write ICV mask as temporary NIfTI in template space
    V_mu = spm_vol(pth_mu);
    V_mu1 = V_mu(1);
    tmp_mask = [tempname '.nii'];
    Vo = V_mu1;
    Vo.fname = tmp_mask;
    Vo.dt = [spm_type('float32') 0];
    Vo.pinfo = [1; 0; 0];
    Vo.descrip = 'ICV mask';
    spm_write_vol(Vo, single(icv_mask));

    % Use SPM batch to warp mask to native space (pull through deformation)
    tmp_out = [tempname '_wmask.nii'];
    matlabbatch = {};
    matlabbatch{1}.spm.util.defs.comp{1}.def = {pth_y};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {tmp_mask};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fileparts(tmp_out)};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;  % nearest-neighbour for binary mask
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 0;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';
    spm_jobman('run', matlabbatch);

    % Read the warped mask — SPM writes it with 'w' prefix in the same dir
    [pth, nam, ext] = fileparts(tmp_mask);
    warped_mask_file = fullfile(pth, ['w' nam ext]);
    if exist(warped_mask_file, 'file')
        V_wm = spm_vol(warped_mask_file);
        mask_warped = spm_read_vols(V_wm);

        % Reslice to native tissue map grid if dimensions differ
        if isequal(V_wm.dim, dm_native) && norm(V_wm.mat - V_target.mat) < 1e-3
            mask_native = mask_warped > 0.5;
        else
            % Resample warped mask onto native tissue map grid
            mask_native = false(dm_native);
            M = V_wm.mat \ V_target.mat;  % native vox → warped mask vox
            [I, J, K] = ndgrid(single(1:dm_native(1)), single(1:dm_native(2)), single(1:dm_native(3)));
            coords = M(1:3,1:3) * [I(:) J(:) K(:)]' + M(1:3,4);
            % Nearest-neighbour sampling
            ci = round(coords(1,:)); cj = round(coords(2,:)); ck = round(coords(3,:));
            valid = ci >= 1 & ci <= V_wm.dim(1) & ...
                    cj >= 1 & cj <= V_wm.dim(2) & ...
                    ck >= 1 & ck <= V_wm.dim(3);
            idx = sub2ind(V_wm.dim, ci(valid), cj(valid), ck(valid));
            vals = mask_warped(idx) > 0.5;
            mask_native_lin = false(1, prod(dm_native));
            mask_native_lin(valid) = vals;
            mask_native = reshape(mask_native_lin, dm_native);
        end

        delete(warped_mask_file);
    else
        warning('Warped mask not found: %s — using all voxels', warped_mask_file);
        mask_native = true(dm_native);
    end

    % Clean up temp file
    if exist(tmp_mask, 'file'), delete(tmp_mask); end
end

function [tbv, tiv] = compute_volumes(tissue_files, mask)
% Compute TBV and TIV from native-space tissue probability maps, masked
% by an intracranial mask.
%   tissue_files: cell array {GM, WM, CSF} of file paths
%   mask:         3D logical array (same dims as tissue maps)
%   TBV = GM + WM, TIV = GM + WM + CSF (in ml)
    V = spm_vol(tissue_files{1});
    vx_vol = abs(det(V.mat(1:3,1:3)));  % voxel volume in mm^3

    tbv = 0;
    tiv = 0;
    for t = 1:3
        Y = spm_read_vols(spm_vol(tissue_files{t}));
        Y(~mask) = 0;
        tissue_ml = sum(Y(:)) * vx_vol / 1000;  % mm^3 -> ml
        if t <= 2
            tbv = tbv + tissue_ml;
        end
        tiv = tiv + tissue_ml;
    end
end
