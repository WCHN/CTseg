function step5_average_normalised(cfg)
% Step 5: Compute average normalised CT images and voxelwise CoV.
% Used for visual and quantitative comparison of normalisation quality.
% A sharper average and lower CoV indicate more consistent normalisation.

    fprintf('=== Step 5: Computing average normalised CTs ===\n');
    subjects = get_subjects(cfg);

    % Accumulators for mean and sum-of-squares (for CoV)
    sum_mr    = [];  sum2_mr    = [];  V_mr    = [];  n_mr    = 0;
    sum_spm   = [];  sum2_spm   = [];  V_spm   = [];  n_spm   = 0;
    sum_ctseg = [];  sum2_ctseg = [];  V_ctseg = [];  n_ctseg = 0;

    for i = 1:numel(subjects)
        subj_dir = fullfile(cfg.data_dir, subjects(i).name);

        % SPM-MR normalised CT (CT warped with MR deformation)
        wct_mr = fullfile(subj_dir, 'wmr_pp_ct_def.nii');
        if ~exist(wct_mr, 'file')
            wct_mr = fullfile(subj_dir, 'wmr_pp_ct.nii');
        end
        if exist(wct_mr, 'file')
            V = spm_vol(wct_mr);
            Y = single(spm_read_vols(V));
            if isempty(sum_mr)
                sum_mr  = zeros(size(Y), 'single');
                sum2_mr = zeros(size(Y), 'single');
                V_mr    = V;
            end
            sum_mr  = sum_mr  + Y;
            sum2_mr = sum2_mr + Y.^2;
            n_mr    = n_mr + 1;
        end

        % SPM-CT normalised CT (prefer deformed CT output)
        wct_spm = fullfile(subj_dir, 'wpp_ct_def.nii');
        if ~exist(wct_spm, 'file')
            wct_spm = fullfile(subj_dir, 'wpp_ct.nii');
        end
        if exist(wct_spm, 'file')
            V = spm_vol(wct_spm);
            Y = single(spm_read_vols(V));
            if isempty(sum_spm)
                sum_spm  = zeros(size(Y), 'single');
                sum2_spm = zeros(size(Y), 'single');
                V_spm    = V;
            end
            sum_spm  = sum_spm  + Y;
            sum2_spm = sum2_spm + Y.^2;
            n_spm    = n_spm + 1;
        end

        % CTseg normalised CT (prefer deformed CT output)
        wct_ctseg = fullfile(subj_dir, 'w_ctseg_pp_ct_def.nii');
        if ~exist(wct_ctseg, 'file')
            wct_ctseg = fullfile(subj_dir, 'w_ctseg_pp_ct.nii');
        end
        if exist(wct_ctseg, 'file')
            V = spm_vol(wct_ctseg);
            Y = single(spm_read_vols(V));
            if isempty(sum_ctseg)
                sum_ctseg  = zeros(size(Y), 'single');
                sum2_ctseg = zeros(size(Y), 'single');
                V_ctseg    = V;
            end
            sum_ctseg  = sum_ctseg  + Y;
            sum2_ctseg = sum2_ctseg + Y.^2;
            n_ctseg    = n_ctseg + 1;
        end
    end

    % Compute and write averages + CoV maps
    norm_results = struct();

    if n_mr > 0
        [avg, cov_map, cov_brain] = compute_avg_and_cov(sum_mr, sum2_mr, n_mr);

        V_mr.dt    = [spm_type('float32') 0];
        V_mr.fname = fullfile(cfg.out_dir, 'avg_normalised_ct_mr.nii');
        spm_write_vol(V_mr, avg);

        V_mr.fname = fullfile(cfg.out_dir, 'cov_normalised_ct_mr.nii');
        spm_write_vol(V_mr, cov_map);

        norm_results.mr_n         = n_mr;
        norm_results.mr_mean_cov  = cov_brain;

        fprintf('  SPM-MR: n=%d  mean brain CoV=%.4f\n', n_mr, cov_brain);
    end

    if n_spm > 0
        [avg, cov_map, cov_brain] = compute_avg_and_cov(sum_spm, sum2_spm, n_spm);

        V_spm.dt    = [spm_type('float32') 0];
        V_spm.fname = fullfile(cfg.out_dir, 'avg_normalised_ct_spm.nii');
        spm_write_vol(V_spm, avg);

        V_spm.fname = fullfile(cfg.out_dir, 'cov_normalised_ct_spm.nii');
        spm_write_vol(V_spm, cov_map);

        norm_results.spm_n         = n_spm;
        norm_results.spm_mean_cov  = cov_brain;

        fprintf('  SPM:   n=%d  mean brain CoV=%.4f\n', n_spm, cov_brain);
    end

    if n_ctseg > 0
        [avg, cov_map, cov_brain] = compute_avg_and_cov(sum_ctseg, sum2_ctseg, n_ctseg);

        V_ctseg.dt    = [spm_type('float32') 0];
        V_ctseg.fname = fullfile(cfg.out_dir, 'avg_normalised_ct_ctseg.nii');
        spm_write_vol(V_ctseg, avg);

        V_ctseg.fname = fullfile(cfg.out_dir, 'cov_normalised_ct_ctseg.nii');
        spm_write_vol(V_ctseg, cov_map);

        norm_results.ctseg_n         = n_ctseg;
        norm_results.ctseg_mean_cov  = cov_brain;

        fprintf('  CTseg: n=%d  mean brain CoV=%.4f\n', n_ctseg, cov_brain);
    end

    save(fullfile(cfg.out_dir, 'normalisation.mat'), 'norm_results');
    fprintf('  Saved to %s\n', fullfile(cfg.out_dir, 'normalisation.mat'));
end


function [avg, cov_map, cov_brain] = compute_avg_and_cov(S, S2, n)
% Compute mean, CoV map, and mean brain CoV from running sums.
    avg = S / n;

    % Sample variance via Bessel-corrected formula, clamped to avoid negative from numerics
    variance = max((S2 - n * avg.^2) / (n - 1), 0);
    sd = sqrt(variance);

    % CoV = sd / |mean|, avoid division by zero
    cov_map = zeros(size(avg), 'single');
    brain_mask = avg > 20 & avg < 100;  % voxels with mean HU in brain parenchyma range
    cov_map(brain_mask) = sd(brain_mask) ./ abs(avg(brain_mask));

    % Mean CoV within brain
    cov_brain = mean(cov_map(brain_mask));
end
