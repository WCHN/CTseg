function grid_search_vsettings()
% Grid search over CTseg v_settings multiplier to find optimal spatial
% regularisation for segmentation.
%
% Uses the best atlas (mu_CTseg_spm15.nii) and tests different
% v_settings multipliers on a single test subject.
%
% Usage:
%   grid_search_vsettings

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
    grid_dir = fullfile(cfg.exp_dir, 'results_grid_vsettings');
    if ~exist(grid_dir, 'dir'), mkdir(grid_dir); end
    cfg.out_dir = grid_dir;

    % Use best template from previous grid search
    dir_ctseg = fileparts(which('spm_CTseg'));
    cfg.pth_mu = fullfile(dir_ctseg, 'mu_CTseg_spm15.nii');
    if ~exist(cfg.pth_mu, 'file')
        error('Template not found: %s', cfg.pth_mu);
    end

    % Grid parameters
    v_mults = [0.5, 1, 2, 3, 4, 5, 6];
    n_v = numel(v_mults);
    n_t = cfg.n_tissues;

    % Preallocate results
    results = struct();
    results.v_mults      = v_mults;
    results.dice_native  = nan(n_v, n_t);
    results.dice_norm    = nan(n_v, n_t);

    % Get test subject name
    subjects = get_subjects(cfg);
    fprintf('Test subject: %s\n', subjects(1).name);
    fprintf('Template: %s\n', cfg.pth_mu);
    fprintf('Grid: %d v_settings multipliers\n\n', n_v);

    % Ensure SPM-MR and SPM-CT outputs exist (needed by step4/4b)
    step1_segment_mr(cfg);
    step2_segment_ct_spm(cfg);

    for vi = 1:n_v
        fprintf('\n========================================\n');
        fprintf('  v_mult=%.1f\n', v_mults(vi));
        fprintf('========================================\n');

        % 1. Set v_settings multiplier
        cfg.v_settings = v_mults(vi);

        % 2. Clean CTseg outputs for test subject
        clean_ctseg_outputs(cfg);

        % 3. Run step3, step4, step4b
        try
            step3_segment_ct_ctseg(cfg);
            step4_compute_metrics(cfg);
            step4b_normalisation_metrics(cfg);

            % 4. Load and store Dice scores
            M = load(fullfile(cfg.out_dir, 'metrics.mat'));
            results.dice_native(vi, :) = M.metrics.dice_ctseg(1, :);

            NM = load(fullfile(cfg.out_dir, 'norm_metrics.mat'));
            results.dice_norm(vi, :) = NM.norm_metrics.dice_ctseg(1, :);
        catch ME
            fprintf('  ** Pipeline FAILED: %s\n', ME.message);
        end
    end

    % Final cleanup
    clean_ctseg_outputs(cfg);

    % Print results
    print_results(results, cfg);

    % Save
    save(fullfile(grid_dir, 'grid_search_results.mat'), 'results');
    fprintf('\nSaved to %s\n', fullfile(grid_dir, 'grid_search_results.mat'));
end


function clean_ctseg_outputs(cfg)
% Delete CTseg-related outputs for test subjects so step3/4/4b re-run.
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
    for f = {'metrics.mat', 'norm_metrics.mat'}
        pth = fullfile(cfg.out_dir, f{1});
        if exist(pth, 'file'), delete(pth); end
    end
end


function print_results(results, cfg)
% Print formatted grid search results.
    n_v = numel(results.v_mults);
    n_t = numel(cfg.tissue_names);

    fprintf('\n\n=== GRID SEARCH RESULTS (v_settings) ===\n\n');

    % Header
    fprintf('%-10s', 'v_mult');
    for t = 1:n_t
        fprintf('  %s_nat  %s_nrm', cfg.tissue_names{t}, cfg.tissue_names{t});
    end
    fprintf('  Mean_nat  Mean_nrm\n');
    fprintf('%s\n', repmat('-', 1, 80));

    best_native = -1; best_vi_n = 0;
    best_norm   = -1; best_vi_nm = 0;

    for vi = 1:n_v
        fprintf('%-10.1f', results.v_mults(vi));
        d_nat = squeeze(results.dice_native(vi, :));
        d_nrm = squeeze(results.dice_norm(vi, :));
        for t = 1:n_t
            fprintf('  %.3f   %.3f ', d_nat(t), d_nrm(t));
        end
        m_nat = mean(d_nat, 'omitnan');
        m_nrm = mean(d_nrm, 'omitnan');
        fprintf('  %.3f     %.3f', m_nat, m_nrm);
        if m_nat > best_native
            best_native = m_nat; best_vi_n = vi;
        end
        if m_nrm > best_norm
            best_norm = m_nrm; best_vi_nm = vi;
        end
        fprintf('\n');
    end

    fprintf('\nBest native:     v_mult=%.1f (mean=%.3f)\n', ...
        results.v_mults(best_vi_n), best_native);
    fprintf('Best normalised: v_mult=%.1f (mean=%.3f)\n', ...
        results.v_mults(best_vi_nm), best_norm);
end
