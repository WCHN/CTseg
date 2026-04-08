function step7_prediction(cfg)
% Step 7: Age and sex prediction from normalised tissue maps.
% Follows the evaluation methodology of Brudfors et al. (arXiv:1909.01140,
% Sec 4.2): concatenated smoothed GM/WM/CSF maps -> GP regression/classification
% with 10-fold CV, using the PredictPRoNTo toolbox.
% Compares predictive performance across methods (SPM-MR, SPM-CT, CTseg-CT).
% Better normalisation + segmentation -> better prediction.
% Saves results to cfg.out_dir/prediction.mat.
%
% Requires: PredictPRoNTo (https://github.com/brudfors/PredictPRoNTo)
%           PRoNTo v2 (bundled in the PredictPRoNTo repo)

    fprintf('=== Step 7: Age/sex prediction ===\n');

    % Load demographics
    demo = load_demographics(cfg);

    % Settings (following Brudfors et al., arXiv:1909.01140, Sec 4.2)
    smooth_fwhm = 8;   % mm (optimised via grid search over 8/12/16 mm)
    n_folds = 10;      % (as in Brudfors 2019)

    methods = {'mr_spm', 'ct_spm', 'ct_ctseg'};
    method_labels = {'SPM-MR', 'SPM-CT', 'CTseg-CT'};

    % Find subjects with all three methods' mwc files and demographics
    subjects = get_subjects(cfg);
    n = numel(subjects);

    has_data = false(n, 1);
    subj_demo_idx = zeros(n, 1);

    for i = 1:n
        subj_dir = fullfile(cfg.data_dir, subjects(i).name);

        % Check all mwc files exist AND are non-empty for all methods
        all_valid = true;
        for mi = 1:numel(methods)
            mwc_f = get_mwc_files(subj_dir, methods{mi});
            for k = 1:3
                if ~exist(mwc_f{k}, 'file')
                    all_valid = false;
                    break;
                end
                % Check non-zero (SPM-CT can produce empty mwc files)
                V = spm_vol(mwc_f{k});
                Y = spm_read_vols(V);
                if sum(Y(:)) == 0
                    all_valid = false;
                    break;
                end
            end
            if ~all_valid, break; end
        end
        if ~all_valid, continue; end

        % Find in demographics
        idx = find(strcmp(demo.ids, subjects(i).name));
        if isempty(idx), continue; end

        has_data(i) = true;
        subj_demo_idx(i) = idx(1);
    end

    valid_subj = find(has_data);
    n_valid = numel(valid_subj);
    fprintf('  %d subjects with complete mwc data\n', n_valid);

    if n_valid < 10
        fprintf('  Too few subjects for prediction (need >= 10) — skipping.\n');
        return;
    end

    % Get ages and sexes for valid subjects
    ages = demo.age(subj_demo_idx(valid_subj));
    sexes = demo.sex(subj_demo_idx(valid_subj));

    has_age = ~isnan(ages);
    has_sex = ~isnan(sexes);
    fprintf('  With age: %d, with sex: %d\n', sum(has_age), sum(has_sex));

    % Adapt number of CV folds to sample size (based on subjects with BOTH
    % age and sex, since that is what gets sent to PredictPRoNTo)
    n_folds = min(n_folds, sum(has_age & has_sex));

    % Store results
    prediction = struct();
    prediction.methods = method_labels;
    prediction.n_subjects = n_valid;
    prediction.n_folds = n_folds;
    prediction.smooth_fwhm = smooth_fwhm;

    % Run PredictPRoNTo for each method
    for mi = 1:numel(methods)
        fprintf('\n  --- %s ---\n', method_labels{mi});

        % Build Data cell array for PredictPRoNTo
        % Column 1: cell of {mwc1, mwc2, mwc3} paths
        % Column 2: age (float) — for regression
        % Column 3: sex (logical) — for classification
        % Only include subjects with BOTH age and sex
        has_both = has_age & has_sex;
        n_both = sum(has_both);
        idx_both = find(has_both);

        if n_both < 2
            fprintf('  Too few subjects with both age and sex — skipping.\n');
            continue;
        end

        % Create symlinks/copies with unique per-subject filenames.
        % PredictPRoNTo's MakeFeatures uses the input filename as the output
        % filename — since all subjects' mwc files have the same basename
        % (e.g., mwc1pp_mr.nii), they overwrite each other in the features dir.
        link_dir = fullfile(cfg.out_dir, ['mwc_links_' methods{mi}]);
        if exist(link_dir, 'dir'), rmdir(link_dir, 's'); end; mkdir(link_dir);

        Data = cell(n_both, 3);
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
            Data{si, 2} = double(ages(idx_both(si)));     % float -> regression
            Data{si, 3} = logical(sexes(idx_both(si)));   % logical -> classification
        end

        % PredictPRoNTo settings
        s = struct();
        s.DirRes      = fullfile(cfg.out_dir, ['pronto_' methods{mi}]);
        s.FWHM        = smooth_fwhm;
        s.CrsVal      = n_folds;
        s.Machine     = 'gp';
        s.ShowResults  = 1;  % print to console only (no figure popups)
        s.DoProcess    = true;
        s.IncBg        = false;
        s.Msk          = true;  % SPM atlas mask
        s.CleanFeatureSet = true;

        % Run
        try
            PredictPRoNTo(Data, s);

            % Extract results from saved PRoNTo outputs
            reg_file = fullfile(s.DirRes, 'Regression.mat');
            cls_file = fullfile(s.DirRes, 'Classification.mat');

            if exist(reg_file, 'file')
                R = load(reg_file);
                PRT = R.ResReg{1}.PRT;
                % Collect predictions across folds
                y_true = []; y_pred = [];
                for f = 1:numel(PRT.model.output.fold)
                    y_true = [y_true; PRT.model.output.fold(f).targets]; %#ok<AGROW>
                    y_pred = [y_pred; PRT.model.output.fold(f).predictions]; %#ok<AGROW>
                end
                err = y_pred - y_true;
                prediction.(methods{mi}).age_true = y_true;
                prediction.(methods{mi}).age_pred = y_pred;
                prediction.(methods{mi}).age_rmse = sqrt(PRT.model.output.stats.mse);
                prediction.(methods{mi}).age_r2   = PRT.model.output.stats.r2;
                prediction.(methods{mi}).age_r    = corr(y_true, y_pred);
                prediction.(methods{mi}).age_mae  = mean(abs(err));
                prediction.(methods{mi}).age_sd   = std(err);
                prediction.(methods{mi}).age_bias = mean(err);
                fprintf('    Age: RMSE=%.2f, MAE=%.2f, r=%.3f, R2=%.3f (n=%d)\n', ...
                    prediction.(methods{mi}).age_rmse, prediction.(methods{mi}).age_mae, ...
                    prediction.(methods{mi}).age_r, prediction.(methods{mi}).age_r2, numel(y_true));
            end

            if exist(cls_file, 'file')
                C = load(cls_file);
                PRT = C.ResClass{1}.PRT;
                stats = PRT.model.output.stats;
                % Collect predictions across folds
                y_true = []; y_score = [];
                for f = 1:numel(PRT.model.output.fold)
                    y_true = [y_true; PRT.model.output.fold(f).targets]; %#ok<AGROW>
                    y_score = [y_score; PRT.model.output.fold(f).predictions]; %#ok<AGROW>
                end
                % PRoNTo targets: 1 for class 1, 2 for class 2
                % PRoNTo predictions: continuous decision values
                % Convert to binary: class 2 = positive (male=1)
                y_bin = double(y_true == 2);  % 0=female, 1=male
                auc = compute_auc_binary(y_bin, y_score);

                prediction.(methods{mi}).sex_true  = y_bin;
                prediction.(methods{mi}).sex_score = y_score;
                prediction.(methods{mi}).sex_pred  = double(y_score > 0);
                prediction.(methods{mi}).cls_file  = cls_file;
                % PRoNTo stats.acc: percentage if >1, proportion if <=1
                if stats.acc > 1
                    prediction.(methods{mi}).sex_acc = stats.acc / 100;
                else
                    prediction.(methods{mi}).sex_acc = stats.acc;
                end
                prediction.(methods{mi}).sex_auc   = auc;

                % Also store confidence interval from PRoNTo
                if isfield(stats, 'acc_lb')
                    if stats.acc_lb > 1
                        prediction.(methods{mi}).sex_acc_lb = stats.acc_lb / 100;
                        prediction.(methods{mi}).sex_acc_ub = stats.acc_ub / 100;
                    else
                        prediction.(methods{mi}).sex_acc_lb = stats.acc_lb;
                        prediction.(methods{mi}).sex_acc_ub = stats.acc_ub;
                    end
                end

                if isfield(stats, 'acc_lb')
                    fprintf('    Sex: AUC=%.3f, Acc=%.1f%% [%.1f, %.1f] (n=%d)\n', ...
                        auc, stats.acc, stats.acc_lb, stats.acc_ub, numel(y_true));
                else
                    fprintf('    Sex: AUC=%.3f, Acc=%.1f%% (n=%d)\n', ...
                        auc, stats.acc, numel(y_true));
                end
            end
        catch ME
            fprintf('  ** FAILED for %s: %s\n', method_labels{mi}, ME.message);
        end

        % Clean up temporary copies
        if exist(link_dir, 'dir'), rmdir(link_dir, 's'); end
    end

    save(fullfile(cfg.out_dir, 'prediction.mat'), 'prediction');
    fprintf('\n  Saved to %s\n', fullfile(cfg.out_dir, 'prediction.mat'));
end


function auc = compute_auc_binary(y_true, y_score)
% Compute AUC. y_true: +1/-1 or 1/0. y_score: continuous predictions.
    % Normalise labels to 0/1
    labels = y_true > 0;
    pos = y_score(labels == 1);
    neg = y_score(labels == 0);
    n_pos = numel(pos);
    n_neg = numel(neg);
    if n_pos == 0 || n_neg == 0
        auc = NaN;
        return;
    end
    auc = 0;
    for i = 1:n_pos
        auc = auc + sum(pos(i) > neg) + 0.5 * sum(pos(i) == neg);
    end
    auc = auc / (n_pos * n_neg);
end
