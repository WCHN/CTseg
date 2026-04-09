function step7_grid_search(cfg)
% Grid search over machine type, FWHM, and GP optimisation iterations
% for the age/sex prediction experiment. All outputs are saved to a
% separate directory (cfg.out_dir/grid_search/) so that existing results
% are not overwritten.
%
% Grid:
%   Machine: 'gp', 'krr', 'rvr'
%   FWHM:    8, 12, 16
%   GP iters: 20 (default), 100  (only for 'gp')
%
% Usage:
%   cfg = config();
%   step7_grid_search(cfg);

    fprintf('=== Grid search: age/sex prediction ===\n');

    % --- Subject filtering (same as step7_prediction.m) ---
    demo = load_demographics(cfg);
    n_folds = 10;

    methods = {'mr_spm', 'ct_spm', 'ct_ctseg'};
    method_labels = {'SPM-MR', 'SPM-CT', 'CTseg-CT'};

    subjects = get_subjects(cfg);
    n = numel(subjects);

    has_data = false(n, 1);
    subj_demo_idx = zeros(n, 1);

    for i = 1:n
        subj_dir = fullfile(cfg.data_dir, subjects(i).name);
        all_valid = true;
        for mi = 1:numel(methods)
            mwc_f = get_mwc_files(subj_dir, methods{mi});
            for k = 1:3
                if ~exist(mwc_f{k}, 'file')
                    all_valid = false; break;
                end
                V = spm_vol(mwc_f{k});
                Y = spm_read_vols(V);
                if sum(Y(:)) == 0
                    all_valid = false; break;
                end
            end
            if ~all_valid, break; end
        end
        if ~all_valid, continue; end
        idx = find(strcmp(demo.ids, subjects(i).name));
        if isempty(idx), continue; end
        has_data(i) = true;
        subj_demo_idx(i) = idx(1);
    end

    valid_subj = find(has_data);
    n_valid = numel(valid_subj);
    fprintf('  %d subjects with complete mwc data\n', n_valid);
    if n_valid < 10
        fprintf('  Too few subjects — skipping.\n'); return;
    end

    ages  = demo.age(subj_demo_idx(valid_subj));
    sexes = demo.sex(subj_demo_idx(valid_subj));
    has_age = ~isnan(ages);
    has_sex = ~isnan(sexes);
    has_both = has_age & has_sex;
    n_both = sum(has_both);
    idx_both = find(has_both);
    n_folds = min(n_folds, n_both);
    fprintf('  %d subjects with both age and sex\n', n_both);

    % --- Define grid ---
    machines  = {'gp', 'krr', 'rvr'};
    fwhms     = [8, 12, 16];
    gp_iters  = [20, 100];

    % Build list of configurations
    configs = {};
    for m = 1:numel(machines)
        for f = 1:numel(fwhms)
            if strcmp(machines{m}, 'gp')
                for gi = 1:numel(gp_iters)
                    configs{end+1} = struct('machine', machines{m}, ...
                        'fwhm', fwhms(f), 'gp_f', gp_iters(gi)); %#ok<AGROW>
                end
            else
                configs{end+1} = struct('machine', machines{m}, ...
                    'fwhm', fwhms(f), 'gp_f', NaN); %#ok<AGROW>
            end
        end
    end
    n_configs = numel(configs);
    fprintf('  %d configurations to test\n\n', n_configs);

    % --- Output directory ---
    grid_dir = fullfile(cfg.out_dir, 'grid_search');
    if ~exist(grid_dir, 'dir'), mkdir(grid_dir); end

    % --- Results storage ---
    results = struct();
    results.configs = configs;
    results.methods = method_labels;
    % Per config × method: age_rmse, age_r, sex_auc
    results.age_rmse = NaN(n_configs, numel(methods));
    results.age_r    = NaN(n_configs, numel(methods));
    results.sex_auc  = NaN(n_configs, numel(methods));

    % --- Run grid ---
    for ci = 1:n_configs
        c = configs{ci};
        if isnan(c.gp_f)
            config_name = sprintf('%s_fwhm%d', c.machine, c.fwhm);
        else
            config_name = sprintf('%s_fwhm%d_f%d', c.machine, c.fwhm, c.gp_f);
        end
        fprintf('--- Config %d/%d: %s ---\n', ci, n_configs, config_name);

        for mi = 1:numel(methods)
            fprintf('  %s: ', method_labels{mi});

            % Create unique file copies per config (same as step7_prediction.m)
            link_dir = fullfile(grid_dir, sprintf('mwc_links_%s_%s', methods{mi}, config_name));
            if exist(link_dir, 'dir'), rmdir(link_dir, 's'); end; mkdir(link_dir);

            % For non-GP machines, only include regression (classification
            % is only implemented for GP in PredictPRoNTo)
            include_cls = strcmp(c.machine, 'gp');
            n_cols = 2 + include_cls;
            Data = cell(n_both, n_cols);
            for si = 1:n_both
                subj_dir = fullfile(cfg.data_dir, subjects(valid_subj(idx_both(si))).name);
                mwc_f = get_mwc_files(subj_dir, methods{mi});
                subj_id = subjects(valid_subj(idx_both(si))).name;
                mwc_unique = cell(1, 3);
                for k = 1:3
                    [~, ~, ext] = fileparts(mwc_f{k});
                    unique_name = fullfile(link_dir, sprintf('%s_mwc%d%s', subj_id, k, ext));
                    copyfile(mwc_f{k}, unique_name);
                    mwc_unique{k} = unique_name;
                end
                Data{si, 1} = mwc_unique;
                Data{si, 2} = double(ages(idx_both(si)));
                if include_cls
                    Data{si, 3} = logical(sexes(idx_both(si)));
                end
            end

            % PredictPRoNTo settings
            s = struct();
            s.DirRes       = fullfile(grid_dir, sprintf('pronto_%s_%s', methods{mi}, config_name));
            s.FWHM         = c.fwhm;
            s.CrsVal       = n_folds;
            s.Machine      = c.machine;
            s.ShowResults  = 0;
            s.DoProcess    = true;
            s.IncBg        = false;
            s.Msk          = true;
            s.CleanFeatureSet = true;
            if strcmp(c.machine, 'gp')
                s.GprArgs = sprintf('-l gauss -h -f %d', c.gp_f);
            end

            try
                PredictPRoNTo(Data, s);

                % Extract age results
                reg_file = fullfile(s.DirRes, 'Regression.mat');
                if exist(reg_file, 'file')
                    R = load(reg_file);
                    PRT = R.ResReg{1}.PRT;
                    y_true = []; y_pred = [];
                    for f = 1:numel(PRT.model.output.fold)
                        y_true = [y_true; PRT.model.output.fold(f).targets]; %#ok<AGROW>
                        y_pred = [y_pred; PRT.model.output.fold(f).predictions]; %#ok<AGROW>
                    end
                    results.age_rmse(ci, mi) = sqrt(mean((y_pred - y_true).^2));
                    results.age_r(ci, mi)    = corr(y_true, y_pred);
                end

                % Extract sex results (AUC from func_val)
                cls_file = fullfile(s.DirRes, 'Classification.mat');
                if exist(cls_file, 'file')
                    C = load(cls_file);
                    [~, ~, auc] = compute_roc_pronto(C.ResClass{1}.PRT, 1);
                    results.sex_auc(ci, mi) = auc;
                end

                fprintf('RMSE=%.2f, r=%.3f, AUC=%.3f\n', ...
                    results.age_rmse(ci, mi), results.age_r(ci, mi), results.sex_auc(ci, mi));
            catch ME
                fprintf('FAILED: %s\n', ME.message);
            end

            % Clean up temp copies
            if exist(link_dir, 'dir'), rmdir(link_dir, 's'); end
        end
        fprintf('\n');
    end

    % --- Print summary table ---
    fprintf('\n=== GRID SEARCH RESULTS ===\n\n');
    fprintf('%-20s', 'Config');
    for mi = 1:numel(methods)
        fprintf('  %s RMSE  %s r   %s AUC', ...
            method_labels{mi}, method_labels{mi}, method_labels{mi});
    end
    fprintf('\n');
    fprintf('%s\n', repmat('-', 1, 20 + numel(methods)*36));
    for ci = 1:n_configs
        c = configs{ci};
        if isnan(c.gp_f)
            config_name = sprintf('%s_fwhm%d', c.machine, c.fwhm);
        else
            config_name = sprintf('%s_fwhm%d_f%d', c.machine, c.fwhm, c.gp_f);
        end
        fprintf('%-20s', config_name);
        for mi = 1:numel(methods)
            fprintf('  %10.2f %6.3f %7.3f', ...
                results.age_rmse(ci, mi), results.age_r(ci, mi), results.sex_auc(ci, mi));
        end
        fprintf('\n');
    end

    % --- Save ---
    save(fullfile(grid_dir, 'grid_results.mat'), 'results');
    fprintf('\nResults saved to %s\n', fullfile(grid_dir, 'grid_results.mat'));
end
