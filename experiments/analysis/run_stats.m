function run_stats(cfg)
% Statistical tests for segmentation metrics and volumetrics.
% Saves results to cfg.out_dir/stats.mat and prints summary.

    fprintf('=== Statistical Analysis ===\n\n');
    alpha = cfg.alpha;

    stats = struct();
    stats.alpha = alpha;

    %% Processing summary: failures and completeness
    fprintf('Processing summary\n');
    fprintf('%s\n', repmat('-', 1, 70));

    subjects = get_subjects(cfg);
    n_total = numel(subjects);
    n_mr = 0; n_spm_ct = 0; n_ctseg = 0; n_complete = 0;
    for i = 1:n_total
        subj_dir = fullfile(cfg.data_dir, subjects(i).name);
        mr_ok    = exist(fullfile(subj_dir, 'c1pp_mr.nii'), 'file') == 2;
        spm_ok   = exist(fullfile(subj_dir, 'wpp_ct_def.nii'), 'file') == 2 || ...
                   exist(fullfile(subj_dir, 'wpp_ct.nii'), 'file') == 2;
        ctseg_ok = exist(fullfile(subj_dir, 'mwc1_ctseg_mni.nii'), 'file') == 2;
        n_mr     = n_mr + mr_ok;
        n_spm_ct = n_spm_ct + spm_ok;
        n_ctseg  = n_ctseg + ctseg_ok;
        n_complete = n_complete + (mr_ok & spm_ok & ctseg_ok);
    end

    stats.n_total    = n_total;
    stats.n_mr       = n_mr;
    stats.n_spm_ct   = n_spm_ct;
    stats.n_ctseg    = n_ctseg;
    stats.n_complete = n_complete;
    stats.n_failed_mr     = n_total - n_mr;
    stats.n_failed_spm_ct = n_total - n_spm_ct;
    stats.n_failed_ctseg  = n_total - n_ctseg;
    stats.n_excluded      = n_total - n_complete;

    fprintf('  Total subjects:     %d\n', n_total);
    fprintf('  SPM-MR succeeded:   %d (%d failed)\n', n_mr, n_total - n_mr);
    fprintf('  SPM-CT succeeded:   %d (%d failed)\n', n_spm_ct, n_total - n_spm_ct);
    fprintf('  CTseg succeeded:    %d (%d failed)\n', n_ctseg, n_total - n_ctseg);
    fprintf('  All three complete: %d (%d excluded)\n\n', n_complete, n_total - n_complete);

    % Load and report failure details if available
    fail_file = fullfile(cfg.out_dir, 'failures.mat');
    if exist(fail_file, 'file')
        F = load(fail_file);
        if ~isempty(F.all_failures)
            fprintf('  Failure details:\n');
            for i = 1:numel(F.all_failures)
                fprintf('    %s [%s]: %s\n', F.all_failures{i}.subject, ...
                    F.all_failures{i}.step, F.all_failures{i}.error);
            end
            fprintf('\n');
        end
    end

    %% Segmentation metrics
    M = load(fullfile(cfg.out_dir, 'metrics.mat'));
    m = M.metrics;
    n = size(m.dice_spm, 1);

    stats.n_subjects = n;

    fprintf('Segmentation metrics (n=%d)\n', n);
    fprintf('%-5s  %-8s  %-22s  %-22s  %-8s  %s\n', ...
        'Tissue', 'Metric', 'SPM (median [IQR])', 'CTseg (median [IQR])', 'p-value', 'Sig.');
    fprintf('%s\n', repmat('-', 1, 90));

    metric_names = {'dice', 'hd95', 'assd'};
    metric_labels = {'Dice', 'HD95', 'ASSD'};

    for mi = 1:numel(metric_names)
        mn = metric_names{mi};
        spm_data   = m.([mn '_spm']);
        ctseg_data = m.([mn '_ctseg']);

        stats.([mn '_p'])     = nan(1, cfg.n_tissues);
        stats.([mn '_spm'])   = struct();
        stats.([mn '_ctseg']) = struct();

        for t = 1:cfg.n_tissues
            x_all = spm_data(:, t);
            y_all = ctseg_data(:, t);

            % Paired subset: only subjects where both methods have data
            ok = ~isnan(x_all) & ~isnan(y_all);
            x_paired = x_all(ok);
            y_paired = y_all(ok);

            % Summary stats computed on paired subset (consistent with test)
            stats.([mn '_spm']).median(t)  = median(x_paired);
            stats.([mn '_spm']).iqr(t,:)   = prctile(x_paired, [25 75]);
            stats.([mn '_spm']).mean(t)    = mean(x_paired);
            stats.([mn '_spm']).std(t)     = std(x_paired);
            stats.([mn '_ctseg']).median(t)  = median(y_paired);
            stats.([mn '_ctseg']).iqr(t,:)   = prctile(y_paired, [25 75]);
            stats.([mn '_ctseg']).mean(t)    = mean(y_paired);
            stats.([mn '_ctseg']).std(t)     = std(y_paired);

            if numel(x_paired) > 1
                p = signrank(x_paired, y_paired);
            else
                p = NaN;
            end
            p_corr = min(p * cfg.n_tissues, 1);
            stats.([mn '_p'])(t) = p_corr;

            % Significance marker
            if strcmp(mn, 'dice')
                sig_str = ternary(p_corr < alpha & median(y_paired) > median(x_paired), '*', '');
            else
                sig_str = ternary(p_corr < alpha & median(y_paired) < median(x_paired), '*', '');
            end

            fprintf('%-5s  %-8s  %6.3f [%5.3f, %5.3f]    %6.3f [%5.3f, %5.3f]    %.4f   %s\n', ...
                cfg.tissue_names{t}, metric_labels{mi}, ...
                median(x_paired), prctile(x_paired, 25), prctile(x_paired, 75), ...
                median(y_paired), prctile(y_paired, 25), prctile(y_paired, 75), ...
                p_corr, sig_str);
        end
    end

    %% Volumetrics
    fprintf('\nVolumetrics\n');
    fprintf('%s\n', repmat('-', 1, 70));

    if exist(fullfile(cfg.out_dir, 'volumetrics.mat'), 'file')
        V = load(fullfile(cfg.out_dir, 'volumetrics.mat'));
        v = V.volumetrics;

        for mn = {'spm', 'ctseg'}
            method = mn{1};
            fprintf('\n  %s vs MR:\n', upper(method));
            for vn = {'tbv', 'tiv'}
                var = vn{1};
                x = v.([var '_mr']);
                y = v.([var '_' method]);
                ok = ~isnan(x) & ~isnan(y);
                x = x(ok); y = y(ok);

                % ICC(3,1)
                [icc_val, icc_lb, icc_ub] = compute_icc(x, y, alpha);

                % Bland-Altman stats
                diff = y - x;
                bias = mean(diff);
                loa  = 1.96 * std(diff);

                % Pearson correlation
                r = corr(x, y);

                prefix = [var '_' method];
                stats.([prefix '_icc'])    = icc_val;
                stats.([prefix '_icc_ci']) = [icc_lb, icc_ub];
                stats.([prefix '_bias'])   = bias;
                stats.([prefix '_loa'])    = [-loa, loa] + bias;
                stats.([prefix '_r'])      = r;

                fprintf('    %s:  ICC=%.3f [%.3f, %.3f]  r=%.3f  bias=%.1f ml  LoA=[%.1f, %.1f] ml  (n=%d)\n', ...
                    upper(var), icc_val, icc_lb, icc_ub, r, bias, bias-loa, bias+loa, numel(x));
            end
        end
        % Runtimes
        fprintf('\n  Runtimes (segmentation only):\n');
        for method = {'spm_mr', 'spm_ct', 'ctseg'}
            field = ['rt_' method{1}];
            if isfield(v, field)
                rt = v.(field);
                rt = rt(~isnan(rt));
                if ~isempty(rt)
                    stats.([field '_mean'])   = mean(rt);
                    stats.([field '_median']) = median(rt);
                    stats.([field '_std'])    = std(rt);
                    stats.([field '_range'])  = [min(rt), max(rt)];
                    fprintf('    %-8s  mean=%.1f +/- %.1f s, median=%.1f s, range=[%.1f, %.1f] s  (n=%d)\n', ...
                        method{1}, mean(rt), std(rt), median(rt), min(rt), max(rt), numel(rt));
                end
            end
        end
    else
        fprintf('  volumetrics.mat not found — skipping.\n');
    end

    %% Normalisation (CoV)
    fprintf('\nSpatial normalisation (CoV)\n');
    fprintf('%s\n', repmat('-', 1, 70));

    if exist(fullfile(cfg.out_dir, 'normalisation.mat'), 'file')
        N = load(fullfile(cfg.out_dir, 'normalisation.mat'));
        nr = N.norm_results;

        if isfield(nr, 'mr_mean_cov')
            stats.norm_mr_cov = nr.mr_mean_cov;
            fprintf('  SPM-MR: mean brain CoV = %.4f  (n=%d)\n', nr.mr_mean_cov, nr.mr_n);
        end
        if isfield(nr, 'spm_mean_cov')
            stats.norm_spm_cov = nr.spm_mean_cov;
            fprintf('  SPM-CT: mean brain CoV = %.4f  (n=%d)\n', nr.spm_mean_cov, nr.spm_n);
        end
        if isfield(nr, 'ctseg_mean_cov')
            stats.norm_ctseg_cov = nr.ctseg_mean_cov;
            fprintf('  CTseg:  mean brain CoV = %.4f  (n=%d)\n', nr.ctseg_mean_cov, nr.ctseg_n);
        end
    else
        fprintf('  normalisation.mat not found — skipping.\n');
    end

    %% Normalisation metrics (step4b)
    fprintf('\nNormalisation metrics (warped MR overlap)\n');
    fprintf('%s\n', repmat('-', 1, 70));

    if exist(fullfile(cfg.out_dir, 'norm_metrics.mat'), 'file')
        NM = load(fullfile(cfg.out_dir, 'norm_metrics.mat'));
        nm = NM.norm_metrics;
        n_nm = size(nm.dice_spm, 1);

        fprintf('  n=%d subjects\n', n_nm);
        fprintf('  %-5s  %-8s  %-22s  %-22s  %-8s  %s\n', ...
            'Tissue', 'Metric', 'SPM-CT (median [IQR])', 'CTseg (median [IQR])', 'p-value', 'Sig.');

        metric_names_nm  = {'dice', 'hd95', 'assd'};
        metric_labels_nm = {'Dice', 'HD95', 'ASSD'};
        for mi = 1:numel(metric_names_nm)
            mn = metric_names_nm{mi};
            spm_data   = nm.([mn '_spm']);
            ctseg_data = nm.([mn '_ctseg']);

            stats.(['norm_' mn '_p']) = nan(1, cfg.n_tissues);

            for t = 1:cfg.n_tissues
                ok = ~isnan(spm_data(:,t)) & ~isnan(ctseg_data(:,t));
                x_paired = spm_data(ok, t);
                y_paired = ctseg_data(ok, t);

                med_x = median(x_paired);  iq_x = prctile(x_paired, [25 75]);
                med_y = median(y_paired);  iq_y = prctile(y_paired, [25 75]);

                stats.(['norm_' mn '_spm']).median(t)  = med_x;
                stats.(['norm_' mn '_spm']).iqr(t,:)   = iq_x;
                stats.(['norm_' mn '_ctseg']).median(t) = med_y;
                stats.(['norm_' mn '_ctseg']).iqr(t,:)  = iq_y;

                if numel(x_paired) > 1
                    p = signrank(x_paired, y_paired);
                else
                    p = NaN;
                end
                p_corr = min(p * cfg.n_tissues, 1);
                stats.(['norm_' mn '_p'])(t) = p_corr;

                if strcmp(mn, 'dice')
                    sig_str = ternary(p_corr < alpha & med_y > med_x, '*', '');
                else
                    sig_str = ternary(p_corr < alpha & med_y < med_x, '*', '');
                end

                fprintf('  %-5s  %-8s  %6.3f [%5.3f, %5.3f]    %6.3f [%5.3f, %5.3f]    %.4f   %s\n', ...
                    cfg.tissue_names{t}, metric_labels_nm{mi}, ...
                    med_x, iq_x(1), iq_x(2), med_y, iq_y(1), iq_y(2), p_corr, sig_str);
            end
        end
    else
        fprintf('  norm_metrics.mat not found — skipping.\n');
    end

    save(fullfile(cfg.out_dir, 'stats.mat'), 'stats');
    fprintf('\nSaved to %s\n', fullfile(cfg.out_dir, 'stats.mat'));
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
