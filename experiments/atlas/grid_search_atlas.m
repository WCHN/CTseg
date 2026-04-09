function grid_search_atlas()
% Grid search over create_atlas parameters to find optimal
% inu_reg and v_settings multipliers for the CTseg MNI-space atlas.
%
% Tests all combinations on a single test subject and reports
% native-space and normalised-space Dice scores.
%
% Usage:
%   grid_search_atlas

    exp_dir = fileparts(fileparts(mfilename('fullpath')));  % experiments/
    addpath(exp_dir);
    addpath(fullfile(exp_dir, 'atlas'));
    addpath(fullfile(exp_dir, 'utils'));
    addpath(fullfile(exp_dir, 'steps'));
    addpath(fullfile(exp_dir, 'analysis'));
    addpath(fullfile(exp_dir, '..'));  % CTseg root

    cfg = config();
    cfg.test_mode = true;
    cfg.test_n_subjects = 1;

    spm('defaults', 'fmri');
    spm_jobman('initcfg');

    % Output to separate folder
    grid_dir = fullfile(cfg.exp_dir, 'results_grid_search');
    if ~exist(grid_dir, 'dir'), mkdir(grid_dir); end
    cfg.out_dir = grid_dir;

    % Grid parameters
    inu_mults = [0.01, 0.1, 1, 10, 100];
    v_mults   = [1, 2, 3, 4];

    n_inu = numel(inu_mults);
    n_v   = numel(v_mults);
    n_t   = cfg.n_tissues;

    % Preallocate results
    results = struct();
    results.inu_mults    = inu_mults;
    results.v_mults      = v_mults;
    results.dice_native  = nan(n_inu, n_v, n_t);
    results.dice_norm    = nan(n_inu, n_v, n_t);

    dir_ctseg = fileparts(which('spm_CTseg'));

    % Get test subject name
    subjects = get_subjects(cfg);
    fprintf('Test subject: %s\n', subjects(1).name);
    fprintf('Grid: %d inu_mults x %d v_mults = %d combinations\n\n', n_inu, n_v, n_inu*n_v);

    % Ensure SPM-MR and SPM-CT outputs exist (needed by step4/4b)
    step1_segment_mr(cfg);
    step2_segment_ct_spm(cfg);

    for ii = 1:n_inu
        for vi = 1:n_v
            fprintf('\n========================================\n');
            fprintf('  inu_mult=%.2f (inu_reg=%.0f), v_mult=%d\n', ...
                inu_mults(ii), 1e4*inu_mults(ii), v_mults(vi));
            fprintf('========================================\n');

            % 1. Create atlas with these params
            pth_mu = fullfile(grid_dir, sprintf('mu_grid_inu%.2f_v%d.nii', ...
                inu_mults(ii), v_mults(vi)));

            try
                create_atlas('mni', 1.5, pth_mu, inu_mults(ii), v_mults(vi));
            catch ME
                fprintf('  ** Atlas creation FAILED: %s\n', ME.message);
                continue;
            end

            % 2. Update config to use this atlas
            cfg.pth_mu = pth_mu;

            % 3. Clean CTseg outputs for test subject
            clean_ctseg_outputs(cfg);

            % 4. Run step3, step4, step4b
            try
                step3_segment_ct_ctseg(cfg);
                step4_compute_metrics(cfg);
                step4b_normalisation_metrics(cfg);

                % 5. Load and store Dice scores
                M = load(fullfile(cfg.out_dir, 'metrics.mat'));
                results.dice_native(ii, vi, :) = M.metrics.dice_ctseg(1, :);

                NM = load(fullfile(cfg.out_dir, 'norm_metrics.mat'));
                results.dice_norm(ii, vi, :) = NM.norm_metrics.dice_ctseg(1, :);
            catch ME
                fprintf('  ** Pipeline FAILED: %s\n', ME.message);
            end

            % 6. Clean up atlas
            if exist(pth_mu, 'file'), delete(pth_mu); end
        end
    end

    % Final cleanup of CTseg outputs
    clean_ctseg_outputs(cfg);

    % Print results
    print_results(results, cfg);

    % Save
    save(fullfile(grid_dir, 'grid_search_results.mat'), 'results');
    fprintf('\nSaved to %s\n', fullfile(grid_dir, 'grid_search_results.mat'));
end


function clean_ctseg_outputs(cfg)
% Delete CTseg-related outputs for test subjects so step3/4/4b re-run.
% Keeps SPM MR/CT outputs (step1/2) intact.
    subjects = get_subjects(cfg);
    patterns = { ...
        'c01_*_CTseg.nii', 'c02_*_CTseg.nii', 'c03_*_CTseg.nii', ...
        'mwc01_*_CTseg.nii', 'mwc02_*_CTseg.nii', 'mwc03_*_CTseg.nii', ...
        'mwc1_ctseg_mni.nii', 'mwc2_ctseg_mni.nii', 'mwc3_ctseg_mni.nii', ...
        'w_ctseg_pp_ct.nii', 'w_ctseg_pp_ct_def.nii', ...
        'wnorm_ctseg_*.nii', 'wnorm_mr_*.nii', 'wnorm_spmct_*.nii', ...
        'y_*_CTseg.nii', 'v_*_CTseg.nii', ...
        'mb_fit_CTseg.mat', 'vol_CTseg.mat'};
    for i = 1:numel(subjects)
        subj_dir = fullfile(cfg.data_dir, subjects(i).name);
        for p = 1:numel(patterns)
            files = dir(fullfile(subj_dir, patterns{p}));
            for f = 1:numel(files)
                delete(fullfile(subj_dir, files(f).name));
            end
        end
    end
    % Also clean metrics files
    for f = {'metrics.mat', 'norm_metrics.mat'}
        pth = fullfile(cfg.out_dir, f{1});
        if exist(pth, 'file'), delete(pth); end
    end
end


function print_results(results, cfg)
% Print formatted grid search results.
    n_inu = numel(results.inu_mults);
    n_v   = numel(results.v_mults);
    n_t   = numel(cfg.tissue_names);

    fprintf('\n\n=== GRID SEARCH RESULTS ===\n\n');

    % Native-space Dice
    fprintf('--- Native-space Dice (CTseg) ---\n');
    fprintf('%-12s', 'inu\v_mult');
    for vi = 1:n_v
        fprintf('  v=%d      ', results.v_mults(vi));
    end
    fprintf('\n');
    best_native = -1; best_ii_n = 0; best_vi_n = 0;
    for ii = 1:n_inu
        fprintf('inu=%-8.2f', results.inu_mults(ii));
        for vi = 1:n_v
            d = squeeze(results.dice_native(ii, vi, :));
            m = mean(d, 'omitnan');
            if m > best_native
                best_native = m; best_ii_n = ii; best_vi_n = vi;
            end
            fprintf('  %.3f     ', m);
        end
        fprintf('\n');
    end
    fprintf('Best: inu_mult=%.2f, v_mult=%d (mean=%.3f)\n\n', ...
        results.inu_mults(best_ii_n), results.v_mults(best_vi_n), best_native);

    % Normalised-space Dice
    fprintf('--- Normalised-space Dice (CTseg) ---\n');
    fprintf('%-12s', 'inu\v_mult');
    for vi = 1:n_v
        fprintf('  v=%d      ', results.v_mults(vi));
    end
    fprintf('\n');
    best_norm = -1; best_ii_nm = 0; best_vi_nm = 0;
    for ii = 1:n_inu
        fprintf('inu=%-8.2f', results.inu_mults(ii));
        for vi = 1:n_v
            d = squeeze(results.dice_norm(ii, vi, :));
            m = mean(d, 'omitnan');
            if m > best_norm
                best_norm = m; best_ii_nm = ii; best_vi_nm = vi;
            end
            fprintf('  %.3f     ', m);
        end
        fprintf('\n');
    end
    fprintf('Best: inu_mult=%.2f, v_mult=%d (mean=%.3f)\n\n', ...
        results.inu_mults(best_ii_nm), results.v_mults(best_vi_nm), best_norm);

    % Per-tissue breakdown for best settings
    fprintf('--- Per-tissue Dice for best native setting (inu=%.2f, v=%d) ---\n', ...
        results.inu_mults(best_ii_n), results.v_mults(best_vi_n));
    for t = 1:n_t
        fprintf('  %-8s  native=%.3f  norm=%.3f\n', cfg.tissue_names{t}, ...
            results.dice_native(best_ii_n, best_vi_n, t), ...
            results.dice_norm(best_ii_n, best_vi_n, t));
    end

    fprintf('\n--- Per-tissue Dice for best normalised setting (inu=%.2f, v=%d) ---\n', ...
        results.inu_mults(best_ii_nm), results.v_mults(best_vi_nm));
    for t = 1:n_t
        fprintf('  %-8s  native=%.3f  norm=%.3f\n', cfg.tissue_names{t}, ...
            results.dice_native(best_ii_nm, best_vi_nm, t), ...
            results.dice_norm(best_ii_nm, best_vi_nm, t));
    end
end
