function step4b_normalisation_metrics(cfg)
% Step 4b: Evaluate normalisation quality by warping the SAME native MR
% tissue maps (silver standard) to SPM/MNI space using three deformations:
%   (1) SPM-MR deformation (y_pp_mr.nii)        — reference
%   (2) SPM-CT deformation (y_pp_ct.nii)         — test
%   (3) CTseg-CT deformation (y_*_CTseg.nii)     — test
% and comparing (2) vs (1) and (3) vs (1) with Dice/HD95/ASSD.
%
% This isolates normalisation quality from segmentation quality:
% the tissue maps are identical, only the spatial mapping differs.
%
% All deformations produce output at 1.5mm in SPM space.
% CTseg deformation maps directly to MNI space when using an MNI-aligned atlas.

    fprintf('=== Step 4b: Normalisation metrics ===\n');
    pth_mu  = resolve_atlas(cfg.mu);
    pth_tpm = fullfile(spm('Dir'), 'tpm', 'TPM.nii');
    subjects = get_subjects(cfg);
    n = numel(subjects);
    thresh = cfg.dice_thresh;


    % Preallocate: SPM-CT vs SPM-MR
    dice_spm   = nan(n, cfg.n_tissues);
    hd95_spm   = nan(n, cfg.n_tissues);
    assd_spm   = nan(n, cfg.n_tissues);
    % CTseg vs SPM-MR
    dice_ctseg = nan(n, cfg.n_tissues);
    hd95_ctseg = nan(n, cfg.n_tissues);
    assd_ctseg = nan(n, cfg.n_tissues);
    valid = false(n, 1);

    for i = 1:n
        subj_dir = fullfile(cfg.data_dir, subjects(i).name);
        mr_files = get_tissue_files(subj_dir, 'mr_spm');

        % Check MR tissue maps exist
        if ~exist(mr_files{1}, 'file'), continue; end

        % Find SPM-MR deformation
        [~, nam_mr] = fileparts(find_nii(subj_dir, 'pp_mr'));
        pth_y_mr = fullfile(subj_dir, ['y_' nam_mr '.nii']);
        if ~exist(pth_y_mr, 'file'), continue; end

        % Find SPM-CT deformation (prefer pp_ct_def, fall back to pp_ct)
        ct_nii = find_nii(subj_dir, 'pp_ct_def');
        if isempty(ct_nii), ct_nii = find_nii(subj_dir, 'pp_ct'); end
        [~, nam_ct] = fileparts(ct_nii);
        pth_y_ct = fullfile(subj_dir, ['y_' nam_ct '.nii']);
        if ~exist(pth_y_ct, 'file'), continue; end

        % Find CTseg deformation
        d_y = dir(fullfile(subj_dir, 'y_*_CTseg.nii'));
        if isempty(d_y), continue; end
        pth_y_ctseg = fullfile(subj_dir, d_y(1).name);

        valid(i) = true;
        fprintf('  [%d/%d] %s\n', i, n, subjects(i).name);

        for t = 1:cfg.n_tissues
            % Warp MR tissue map with all three deformations
            wc_mr    = warp_to_mni_spm(mr_files{t}, pth_y_mr, pth_tpm, subj_dir, 'wnorm_mr_');
            wc_spmct = warp_to_mni_spm(mr_files{t}, pth_y_ct, pth_tpm, subj_dir, 'wnorm_spmct_');
            wc_ctseg = warp_to_mni_ctseg(mr_files{t}, pth_y_ctseg, pth_tpm, subj_dir);

            % Read warped images
            V_mr    = spm_vol(wc_mr);
            V_spmct = spm_vol(wc_spmct);
            V_ctseg = spm_vol(wc_ctseg);
            Y_mr    = spm_read_vols(V_mr);
            Y_spmct = spm_read_vols(V_spmct);
            Y_ctseg = spm_read_vols(V_ctseg);

            % Binarise
            vx = sqrt(sum(V_mr.mat(1:3,1:3).^2));
            seg_mr    = Y_mr > thresh;
            seg_spmct = Y_spmct > thresh;
            seg_ctseg = Y_ctseg > thresh;

            % SPM-CT vs SPM-MR
            dice_spm(i, t) = compute_dice(seg_mr, seg_spmct);
            [hd95_spm(i, t), assd_spm(i, t)] = compute_surface_distances(seg_mr, seg_spmct, vx);

            % CTseg vs SPM-MR
            dice_ctseg(i, t) = compute_dice(seg_mr, seg_ctseg);
            [hd95_ctseg(i, t), assd_ctseg(i, t)] = compute_surface_distances(seg_mr, seg_ctseg, vx);

            % Clean up temporary warped files
            delete(wc_mr); delete(wc_spmct); delete(wc_ctseg);
        end

        fprintf('    Dice SPM-CT: [%.3f %.3f %.3f]  CTseg: [%.3f %.3f %.3f]\n', ...
            dice_spm(i,:), dice_ctseg(i,:));
    end

    % Save
    norm_metrics = struct();
    norm_metrics.subject_names = {subjects(valid).name};
    norm_metrics.tissue_names  = cfg.tissue_names;
    norm_metrics.dice_spm      = dice_spm(valid, :);
    norm_metrics.hd95_spm      = hd95_spm(valid, :);
    norm_metrics.assd_spm      = assd_spm(valid, :);
    norm_metrics.dice_ctseg    = dice_ctseg(valid, :);
    norm_metrics.hd95_ctseg    = hd95_ctseg(valid, :);
    norm_metrics.assd_ctseg    = assd_ctseg(valid, :);

    fprintf('  %d/%d subjects complete\n', sum(valid), n);

    % Print mean Dice per tissue class
    fprintf('\n  Mean Dice:\n');
    fprintf('  %-8s  %-10s  %-10s\n', 'Tissue', 'SPM-CT', 'CTseg');
    for t = 1:cfg.n_tissues
        fprintf('  %-8s  %.3f       %.3f\n', cfg.tissue_names{t}, ...
            mean(norm_metrics.dice_spm(:,t), 'omitnan'), ...
            mean(norm_metrics.dice_ctseg(:,t), 'omitnan'));
    end
    fprintf('  %-8s  %.3f       %.3f\n', 'Mean', ...
        mean(norm_metrics.dice_spm(:), 'omitnan'), ...
        mean(norm_metrics.dice_ctseg(:), 'omitnan'));

    save(fullfile(cfg.out_dir, 'norm_metrics.mat'), 'norm_metrics');
    fprintf('\n  Saved to %s\n', fullfile(cfg.out_dir, 'norm_metrics.mat'));
end


function pth_out = warp_to_mni_spm(tissue_file, def_file, pth_mu, subj_dir, prefix)
% Warp a native tissue map to MNI using SPM normalise-write.
% SPM deformations are defined on the template grid → need pull.
    matlabbatch = {};
    matlabbatch{1}.spm.util.defs.comp{1}.comp{1}.def         = {def_file};
    matlabbatch{1}.spm.util.defs.comp{1}.space               = {pth_mu};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {tissue_file};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {subj_dir};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = prefix;
    spm_jobman('run',matlabbatch);

    [~, nam, ext] = fileparts(tissue_file);
    pth_out = fullfile(subj_dir, [prefix nam ext]);
end


function pth_out = warp_to_mni_ctseg(tissue_file, def_file, pth_mu, subj_dir)
% Warp a native tissue map to SPM space using CTseg deformation.
% CTseg deformations are defined on the native grid → need push.
    matlabbatch = {};
    matlabbatch{1}.spm.util.defs.comp{1}.def                 = {def_file};
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames          = {tissue_file};
    matlabbatch{1}.spm.util.defs.out{1}.push.weight          = {''};
    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {subj_dir};
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.file        = {pth_mu};
    matlabbatch{1}.spm.util.defs.out{1}.push.prefix          = 'wnorm_ctseg_';
    spm_jobman('run',matlabbatch);

    [~, nam, ext] = fileparts(tissue_file);
    pth_out = fullfile(subj_dir, ['wnorm_ctseg_' nam ext]);
end