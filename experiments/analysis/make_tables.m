function make_tables(cfg)
% Generate LaTeX tables for the paper.
% Saves .tex files to cfg.fig_dir.

    fprintf('=== Generating LaTeX tables ===\n');

    S = load(fullfile(cfg.out_dir, 'stats.mat'));
    s = S.stats;

    %% Table 1: Segmentation metrics
    fid = fopen(fullfile(cfg.fig_dir, 'table_segmentation.tex'), 'w');
    fprintf(fid, '\\begin{table}[t]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{Segmentation accuracy compared to MRI silver standard. ');
    fprintf(fid, 'Values are median [IQR]. $p$-values from Wilcoxon signed-rank test (Bonferroni-corrected).}\n');
    fprintf(fid, '\\label{tab:segmentation}\n');
    fprintf(fid, '\\small\n');
    fprintf(fid, '\\begin{tabular}{llccc}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Metric & Method & GM & WM & CSF \\\\\n');
    fprintf(fid, '\\midrule\n');

    metric_names = {'dice', 'hd95', 'assd'};
    metric_labels = {'Dice', 'HD95 (mm)', 'ASSD (mm)'};

    for mi = 1:numel(metric_names)
        mn = metric_names{mi};

        % SPM row
        fprintf(fid, '%s & SPM', metric_labels{mi});
        for t = 1:cfg.n_tissues
            med = s.([mn '_spm']).median(t);
            iq  = s.([mn '_spm']).iqr(t,:);
            fprintf(fid, ' & %.3f [%.3f, %.3f]', med, iq(1), iq(2));
        end
        fprintf(fid, ' \\\\\n');

        % CTseg row
        fprintf(fid, ' & CTseg');
        for t = 1:cfg.n_tissues
            med = s.([mn '_ctseg']).median(t);
            iq  = s.([mn '_ctseg']).iqr(t,:);
            fprintf(fid, ' & %.3f [%.3f, %.3f]', med, iq(1), iq(2));
        end
        fprintf(fid, ' \\\\\n');

        % p-value row
        fprintf(fid, ' & $p$-value');
        for t = 1:cfg.n_tissues
            p = s.([mn '_p'])(t);
            if p < 0.001
                fprintf(fid, ' & $<$0.001');
            else
                fprintf(fid, ' & %.3f', p);
            end
        end
        fprintf(fid, ' \\\\\n');

        if mi < numel(metric_names)
            fprintf(fid, '\\midrule\n');
        end
    end

    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
    fprintf('  Saved table_segmentation.tex\n');

    %% Table 2: Volumetrics
    if isfield(s, 'tbv_spm_icc')
        fid = fopen(fullfile(cfg.fig_dir, 'table_volumetrics.tex'), 'w');
        fprintf(fid, '\\begin{table}[t]\n');
        fprintf(fid, '\\centering\n');
        fprintf(fid, '\\caption{Agreement of brain volume estimates (CT methods vs.\\ MRI reference).}\n');
        fprintf(fid, '\\label{tab:volumetrics}\n');
        fprintf(fid, '\\small\n');
        fprintf(fid, '\\begin{tabular}{llcccc}\n');
        fprintf(fid, '\\toprule\n');
        fprintf(fid, 'Measure & Method & ICC [95\\%% CI] & $r$ & Bias (ml) & LoA (ml) \\\\\n');
        fprintf(fid, '\\midrule\n');

        for vn = {'tbv', 'tiv'}
            var = vn{1};
            for mn = {'spm', 'ctseg'}
                method = mn{1};
                if strcmp(method, 'spm')
                    method_label = 'SPM';
                else
                    method_label = 'CTseg';
                end
                prefix = [var '_' method];
                fprintf(fid, '%s & %s & %.3f [%.3f, %.3f] & %.3f & %.1f & [%.1f, %.1f] \\\\\n', ...
                    upper(var), method_label, ...
                    s.([prefix '_icc']), s.([prefix '_icc_ci'])(1), s.([prefix '_icc_ci'])(2), ...
                    s.([prefix '_r']), ...
                    s.([prefix '_bias']), ...
                    s.([prefix '_loa'])(1), s.([prefix '_loa'])(2));
            end
            if strcmp(var, 'tbv')
                fprintf(fid, '\\midrule\n');
            end
        end

        fprintf(fid, '\\bottomrule\n');
        fprintf(fid, '\\end{tabular}\n');
        fprintf(fid, '\\end{table}\n');
        fclose(fid);
        fprintf('  Saved table_volumetrics.tex\n');
    end

    %% Table 3: Prediction (age regression + sex classification)
    pred_file = fullfile(cfg.out_dir, 'prediction.mat');
    if exist(pred_file, 'file')
        P = load(pred_file);
        p = P.prediction;

        methods = {'mr_spm', 'ct_spm', 'ct_ctseg'};
        method_labels = p.methods;

        fid = fopen(fullfile(cfg.fig_dir, 'table_prediction.tex'), 'w');
        fprintf(fid, '\\begin{table}[t]\n');
        fprintf(fid, '\\centering\n');
        fprintf(fid, '\\small\n');
        fprintf(fid, '\\caption{Age and sex prediction from normalised tissue maps ');
        fprintf(fid, '(10-fold CV, GP models, FWHM=%d\\,mm). ', p.smooth_fwhm);
        fprintf(fid, 'Following the evaluation paradigm of \\citet{brudfors2019superres}.}\n');
        fprintf(fid, '\\label{tab:prediction}\n');
        fprintf(fid, '\\small\n');
        fprintf(fid, '\\begin{tabular}{lcccccc}\n');
        fprintf(fid, '\\toprule\n');
        fprintf(fid, ' & \\multicolumn{4}{c}{Age regression} & \\multicolumn{2}{c}{Sex classification} \\\\\n');
        fprintf(fid, '\\cmidrule(lr){2-5} \\cmidrule(lr){6-7}\n');
        fprintf(fid, 'Method & RMSE (yr) & SD (yr) & Bias (yr) & $r$ & Acc (\\%%) & AUC \\\\\n');
        fprintf(fid, '\\midrule\n');

        for mi = 1:numel(methods)
            if ~isfield(p, methods{mi}), continue; end
            ms = p.(methods{mi});
            fprintf(fid, '%s', method_labels{mi});

            % Age metrics
            if isfield(ms, 'age_rmse')
                fprintf(fid, ' & %.2f & %.2f & %.2f & %.3f', ...
                    ms.age_rmse, ms.age_sd, ms.age_bias, ms.age_r);
            else
                fprintf(fid, ' & -- & -- & -- & --');
            end

            % Sex metrics
            if isfield(ms, 'sex_acc')
                % Convert to percentage; handle both proportion and
                % double-divided values from earlier step7 bug
                if ms.sex_acc < 0.01
                    acc_pct = ms.sex_acc * 10000;  % was divided by 100 twice
                elseif ms.sex_acc <= 1
                    acc_pct = ms.sex_acc * 100;    % proportion → percentage
                else
                    acc_pct = ms.sex_acc;           % already percentage
                end
                % Recompute AUC from func_val (continuous GP posteriors)
                % to match prt_plot_ROC / make_figures
                auc = ms.sex_auc;
                if isfield(ms, 'cls_file') && exist(ms.cls_file, 'file')
                    C = load(ms.cls_file);
                    [~, ~, auc] = compute_roc_pronto(C.ResClass{1}.PRT, 1);
                end
                fprintf(fid, ' & %.1f & %.3f', acc_pct, auc);
            else
                fprintf(fid, ' & -- & --');
            end

            fprintf(fid, ' \\\\\n');
        end

        fprintf(fid, '\\bottomrule\n');
        fprintf(fid, '\\end{tabular}\n');
        fprintf(fid, '\\end{table}\n');
        fclose(fid);
        fprintf('  Saved table_prediction.tex\n');
    end
end
