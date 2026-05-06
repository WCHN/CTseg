function step4_compute_metrics(cfg)
% Step 4: Compute segmentation metrics.
% For each subject, compares CT-SPM and CT-CTseg tissue maps against
% MR-SPM silver standard. Computes Dice, Hausdorff95, and ASSD.
% Only includes subjects where ALL three methods produced output.
% Saves results to cfg.out_dir/metrics.mat.

    fprintf('=== Step 4: Computing segmentation metrics ===\n');
    subjects = get_subjects(cfg);
    n = numel(subjects);
    thresh = cfg.dice_thresh;

    % Preallocate: subjects x tissues
    dice_spm   = nan(n, cfg.n_tissues);
    dice_ctseg = nan(n, cfg.n_tissues);
    hd95_spm   = nan(n, cfg.n_tissues);
    hd95_ctseg = nan(n, cfg.n_tissues);
    assd_spm   = nan(n, cfg.n_tissues);
    assd_ctseg = nan(n, cfg.n_tissues);
    valid = false(n, 1);
    n_skipped = 0;

    for i = 1:n
        subj_dir = fullfile(cfg.data_dir, subjects(i).name);
        mr_files    = get_tissue_files(subj_dir, 'mr_spm');
        ct_spm_f    = get_tissue_files(subj_dir, 'ct_spm');
        ct_ctseg_f  = get_tissue_files(subj_dir, 'ct_ctseg');

        % Require ALL three methods to have completed
        all_exist = exist(mr_files{1}, 'file') && ...
                    exist(ct_spm_f{1}, 'file') && ...
                    exist(ct_ctseg_f{1}, 'file');
        if ~all_exist
            n_skipped = n_skipped + 1;
            fprintf('  [%d/%d] %s - SKIP (incomplete: MR=%d SPM-CT=%d CTseg=%d)\n', ...
                i, n, subjects(i).name, ...
                exist(mr_files{1}, 'file') == 2, ...
                exist(ct_spm_f{1}, 'file') == 2, ...
                exist(ct_ctseg_f{1}, 'file') == 2);
            continue;
        end
        valid(i) = true;

        % Get voxel size from MR reference
        V_mr = spm_vol(mr_files{1});
        vx = sqrt(sum(V_mr.mat(1:3,1:3).^2));

        for t = 1:cfg.n_tissues
            mr_vol = spm_read_vols(spm_vol(mr_files{t}));
            mr_seg = mr_vol > thresh;

            % CT-SPM
            ct_vol = spm_read_vols(spm_vol(ct_spm_f{t}));
            ct_seg = ct_vol > thresh;
            dice_spm(i, t) = compute_dice(mr_seg, ct_seg);
            [hd95_spm(i, t), assd_spm(i, t)] = compute_surface_distances(mr_seg, ct_seg, vx);

            % CT-CTseg
            ct_vol = spm_read_vols(spm_vol(ct_ctseg_f{t}));
            ct_seg = ct_vol > thresh;
            dice_ctseg(i, t) = compute_dice(mr_seg, ct_seg);
            [hd95_ctseg(i, t), assd_ctseg(i, t)] = compute_surface_distances(mr_seg, ct_seg, vx);
        end

        fprintf('  [%d/%d] %s  Dice SPM: [%.3f %.3f %.3f]  CTseg: [%.3f %.3f %.3f]\n', ...
            i, n, subjects(i).name, dice_spm(i,:), dice_ctseg(i,:));
    end

    % Save only complete subjects
    metrics = struct();
    metrics.subject_names = {subjects(valid).name};
    metrics.tissue_names  = cfg.tissue_names;
    metrics.dice_spm      = dice_spm(valid, :);
    metrics.dice_ctseg    = dice_ctseg(valid, :);
    metrics.hd95_spm      = hd95_spm(valid, :);
    metrics.hd95_ctseg    = hd95_ctseg(valid, :);
    metrics.assd_spm      = assd_spm(valid, :);
    metrics.assd_ctseg    = assd_ctseg(valid, :);

    fprintf('  %d/%d subjects complete, %d skipped (incomplete)\n', ...
        sum(valid), n, n_skipped);
    save(fullfile(cfg.out_dir, 'metrics.mat'), 'metrics');
    fprintf('  Saved to %s\n', fullfile(cfg.out_dir, 'metrics.mat'));
end
