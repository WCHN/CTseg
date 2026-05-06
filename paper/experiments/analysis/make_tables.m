function make_tables(cfg)
% Generate LaTeX tables for the paper.
% Saves .tex files to cfg.fig_dir.
%
% Within each comparison column, the best value is bolded:
%   higher-is-better: Dice, ICC, r, Acc, AUC
%   lower-is-better:  HD95, ASSD, RMSE, SD, |Bias|, LoA range

    fprintf('=== Generating LaTeX tables ===\n');

    S = load(fullfile(cfg.out_dir, 'stats.mat'));
    s = S.stats;

    %% Table 1: Segmentation metrics
    fid = fopen(fullfile(cfg.fig_dir, 'table_segmentation.tex'), 'w');
    fprintf(fid, '\\begin{table}[t]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{Segmentation accuracy compared to MRI silver standard. ');
    fprintf(fid, 'Values are median [IQR]. $p$-values from Wilcoxon signed-rank test (Bonferroni-corrected). Best value per column in bold.}\n');
    fprintf(fid, '\\label{tab:segmentation}\n');
    fprintf(fid, '\\small\n');
    fprintf(fid, '\\begin{tabular}{llccc}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Metric & Method & GM & WM & CSF \\\\\n');
    fprintf(fid, '\\midrule\n');

    metric_names  = {'dice', 'hd95', 'assd'};
    metric_labels = {'Dice', 'HD95 (mm)', 'ASSD (mm)'};
    higher_better = [true, false, false];  % Dice up, HD95/ASSD down

    for mi = 1:numel(metric_names)
        mn = metric_names{mi};
        spm_med   = s.([mn '_spm']).median;
        ctseg_med = s.([mn '_ctseg']).median;

        % Per-tissue: which row is best?
        if higher_better(mi)
            ctseg_best = ctseg_med > spm_med;
        else
            ctseg_best = ctseg_med < spm_med;
        end

        % SPM row
        fprintf(fid, '%s & SPM', metric_labels{mi});
        for t = 1:cfg.n_tissues
            cell_str = sprintf('%.3f [%.3f, %.3f]', s.([mn '_spm']).median(t), ...
                s.([mn '_spm']).iqr(t,1), s.([mn '_spm']).iqr(t,2));
            fprintf(fid, ' & %s', bold_if(cell_str, ~ctseg_best(t)));
        end
        fprintf(fid, ' \\\\\n');

        % CTseg row
        fprintf(fid, ' & CTseg');
        for t = 1:cfg.n_tissues
            cell_str = sprintf('%.3f [%.3f, %.3f]', s.([mn '_ctseg']).median(t), ...
                s.([mn '_ctseg']).iqr(t,1), s.([mn '_ctseg']).iqr(t,2));
            fprintf(fid, ' & %s', bold_if(cell_str, ctseg_best(t)));
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
        fprintf(fid, '\\caption{Agreement of brain volume estimates (CT methods vs.\\ MRI reference). Best value per column (within each measure) in bold; for Bias, smaller magnitude is better, for LoA, tighter interval is better.}\n');
        fprintf(fid, '\\label{tab:volumetrics}\n');
        fprintf(fid, '\\small\n');
        fprintf(fid, '\\begin{tabular}{llcccc}\n');
        fprintf(fid, '\\toprule\n');
        fprintf(fid, 'Measure & Method & ICC [95\\%% CI] & $r$ & Bias (ml) & LoA (ml) \\\\\n');
        fprintf(fid, '\\midrule\n');

        vol_vars = {'tbv', 'tiv'};
        for vi = 1:numel(vol_vars)
            var = vol_vars{vi};
            spm_pref   = [var '_spm'];
            ctseg_pref = [var '_ctseg'];

            spm_icc   = s.([spm_pref   '_icc']);
            ctseg_icc = s.([ctseg_pref '_icc']);
            spm_r     = s.([spm_pref   '_r']);
            ctseg_r   = s.([ctseg_pref '_r']);
            spm_bias   = s.([spm_pref   '_bias']);
            ctseg_bias = s.([ctseg_pref '_bias']);
            spm_loa   = s.([spm_pref   '_loa']);
            ctseg_loa = s.([ctseg_pref '_loa']);

            % Best row per column
            ctseg_best_icc  = ctseg_icc > spm_icc;
            ctseg_best_r    = ctseg_r > spm_r;
            ctseg_best_bias = abs(ctseg_bias) < abs(spm_bias);
            spm_loa_range   = spm_loa(2)   - spm_loa(1);
            ctseg_loa_range = ctseg_loa(2) - ctseg_loa(1);
            ctseg_best_loa  = ctseg_loa_range < spm_loa_range;

            % SPM row
            spm_icc_str  = sprintf('%.3f [%.3f, %.3f]', spm_icc, s.([spm_pref '_icc_ci'])(1), s.([spm_pref '_icc_ci'])(2));
            spm_r_str    = sprintf('%.3f', spm_r);
            spm_bias_str = sprintf('%.1f',  spm_bias);
            spm_loa_str  = sprintf('[%.1f, %.1f]', spm_loa(1), spm_loa(2));
            fprintf(fid, '%s & SPM & %s & %s & %s & %s \\\\\n', ...
                upper(var), ...
                bold_if(spm_icc_str,  ~ctseg_best_icc), ...
                bold_if(spm_r_str,    ~ctseg_best_r), ...
                bold_if(spm_bias_str, ~ctseg_best_bias), ...
                bold_if(spm_loa_str,  ~ctseg_best_loa));

            % CTseg row
            ctseg_icc_str  = sprintf('%.3f [%.3f, %.3f]', ctseg_icc, s.([ctseg_pref '_icc_ci'])(1), s.([ctseg_pref '_icc_ci'])(2));
            ctseg_r_str    = sprintf('%.3f', ctseg_r);
            ctseg_bias_str = sprintf('%.1f',  ctseg_bias);
            ctseg_loa_str  = sprintf('[%.1f, %.1f]', ctseg_loa(1), ctseg_loa(2));
            fprintf(fid, '%s & CTseg & %s & %s & %s & %s \\\\\n', ...
                upper(var), ...
                bold_if(ctseg_icc_str,  ctseg_best_icc), ...
                bold_if(ctseg_r_str,    ctseg_best_r), ...
                bold_if(ctseg_bias_str, ctseg_best_bias), ...
                bold_if(ctseg_loa_str,  ctseg_best_loa));

            if vi < numel(vol_vars)
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

        % Collect column values for best-finding
        n_m = numel(methods);
        rmse_v = nan(1,n_m); sd_v = nan(1,n_m); bias_v = nan(1,n_m);
        rcorr_v = nan(1,n_m); acc_v = nan(1,n_m); auc_v = nan(1,n_m);
        for mi = 1:n_m
            if ~isfield(p, methods{mi}), continue; end
            ms = p.(methods{mi});
            if isfield(ms, 'age_rmse')
                rmse_v(mi)  = ms.age_rmse;
                sd_v(mi)    = ms.age_sd;
                bias_v(mi)  = ms.age_bias;
                rcorr_v(mi) = ms.age_r;
            end
            if isfield(ms, 'sex_acc')
                if ms.sex_acc < 0.01
                    acc_v(mi) = ms.sex_acc * 10000;
                elseif ms.sex_acc <= 1
                    acc_v(mi) = ms.sex_acc * 100;
                else
                    acc_v(mi) = ms.sex_acc;
                end
                auc = ms.sex_auc;
                if isfield(ms, 'cls_file') && exist(ms.cls_file, 'file')
                    C = load(ms.cls_file);
                    [~, ~, auc] = compute_roc_pronto(C.ResClass{1}.PRT, 1);
                end
                auc_v(mi) = auc;
            end
        end
        [~, best_rmse] = min(rmse_v);
        [~, best_sd]   = min(sd_v);
        [~, best_bias] = min(abs(bias_v));
        [~, best_r]    = max(rcorr_v);
        [~, best_acc]  = max(acc_v);
        [~, best_auc]  = max(auc_v);

        fid = fopen(fullfile(cfg.fig_dir, 'table_prediction.tex'), 'w');
        fprintf(fid, '\\begin{table}[t]\n');
        fprintf(fid, '\\centering\n');
        fprintf(fid, '\\small\n');
        fprintf(fid, '\\caption{Age and sex prediction from normalised tissue maps ');
        fprintf(fid, '(10-fold CV, GP models, FWHM=%d\\,mm). ', p.smooth_fwhm);
        fprintf(fid, 'Following the evaluation paradigm of \\citet{brudfors2019superres}. Best value per column in bold; for Bias, smaller magnitude is better.}\n');
        fprintf(fid, '\\label{tab:prediction}\n');
        fprintf(fid, '\\small\n');
        fprintf(fid, '\\begin{tabular}{lcccccc}\n');
        fprintf(fid, '\\toprule\n');
        fprintf(fid, ' & \\multicolumn{4}{c}{Age regression} & \\multicolumn{2}{c}{Sex classification} \\\\\n');
        fprintf(fid, '\\cmidrule(lr){2-5} \\cmidrule(lr){6-7}\n');
        fprintf(fid, 'Method & RMSE (yr) & SD (yr) & Bias (yr) & $r$ & Acc (\\%%) & AUC \\\\\n');
        fprintf(fid, '\\midrule\n');

        for mi = 1:n_m
            if ~isfield(p, methods{mi}), continue; end
            fprintf(fid, '%s', method_labels{mi});

            % Age metrics
            if ~isnan(rmse_v(mi))
                fprintf(fid, ' & %s & %s & %s & %s', ...
                    bold_if(sprintf('%.2f', rmse_v(mi)),  mi == best_rmse), ...
                    bold_if(sprintf('%.2f', sd_v(mi)),    mi == best_sd), ...
                    bold_if(sprintf('%.2f', bias_v(mi)),  mi == best_bias), ...
                    bold_if(sprintf('%.3f', rcorr_v(mi)), mi == best_r));
            else
                fprintf(fid, ' & -- & -- & -- & --');
            end

            % Sex metrics
            if ~isnan(acc_v(mi))
                fprintf(fid, ' & %s & %s', ...
                    bold_if(sprintf('%.1f', acc_v(mi)), mi == best_acc), ...
                    bold_if(sprintf('%.3f', auc_v(mi)), mi == best_auc));
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


function out = bold_if(str, is_best)
% Wrap str in \textbf{} if is_best is true.
    if is_best
        out = ['\textbf{' str '}'];
    else
        out = str;
    end
end
