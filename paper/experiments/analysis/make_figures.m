function make_figures(cfg, only)
% Generate publication-quality figures.
% Saves .png and .pdf to cfg.fig_dir.
%
% Usage:
%   make_figures(cfg)        % generate all figures
%   make_figures(cfg, 0)     % generate only figure 0
%   make_figures(cfg, [4 5]) % generate figures 4 and 5
%
% Figure numbers:
%   0  — Example MR and CT data
%   0.5— CTseg template and intensity priors (use 0.5)
%   1  — Segmentation metrics (boxplots)
%   1.5— Normalisation metrics (use 1.5)
%   2  — Bland-Altman plots (volumetrics)
%   3  — Average normalised CTs
%   4  — Best and worst Dice examples
%   5  — Example tissue maps
%   6  — Demographics distributions
%   7  — Prediction results

    if nargin < 2, only = []; end
    do = @(n) isempty(only) || any(abs(only - n) < 0.01);

    fprintf('=== Generating figures ===\n');

    subjects = get_subjects(cfg);
    pth_mu   = resolve_atlas(cfg.mu);

    %% Figure 0: Example MR and CT data
    if do(0)
    % Show one subject with axial, coronal, and sagittal views for both MR and CT
    % Pick a subject from the middle of the dataset
    mid_idx = round(numel(subjects) / 2);
    subj_dir = fullfile(cfg.data_dir, subjects(mid_idx).name);

    mr_gz  = fullfile(subj_dir, 'pp_mr.nii.gz');
    mr_nii = fullfile(subj_dir, 'pp_mr.nii');
    ct_def = fullfile(subj_dir, 'pp_ct_def.nii');
    ct_gz   = fullfile(subj_dir, 'pp_ct.nii.gz');
    ct_nii  = fullfile(subj_dir, 'pp_ct.nii');
    if exist(ct_def, 'file'), ct_nii = ct_def; end
    cleanup_mr = false; cleanup_ct = false;
    if ~exist(mr_nii, 'file') && exist(mr_gz, 'file')
        gunzip(mr_gz, subj_dir); cleanup_mr = true;
    end
    if ~exist(ct_nii, 'file') && exist(ct_gz, 'file')
        gunzip(ct_gz, subj_dir); cleanup_ct = true;
    end

    if exist(mr_nii, 'file') && exist(ct_nii, 'file')
        mr_vol = spm_read_vols(spm_vol(mr_nii));
        ct_vol = spm_read_vols(spm_vol(ct_nii));
        dims = size(mr_vol);

        % Slice indices (roughly mid-brain)
        ax = round(0.45 * dims(3));  % axial (z)
        co = round(0.55 * dims(2));  % coronal (y)
        sa = round(0.50 * dims(1));  % sagittal (x)

        mr_slices = {rot90(squeeze(mr_vol(:,:,ax))), ...
                     rot90(squeeze(mr_vol(:,co,:))), ...
                     rot90(squeeze(mr_vol(sa,:,:)))};
        ct_slices = {rot90(squeeze(ct_vol(:,:,ax))), ...
                     rot90(squeeze(ct_vol(:,co,:))), ...
                     rot90(squeeze(ct_vol(sa,:,:)))};

        fig0 = figure('Name', 'Example Data', 'Position', [50 50 900 450], 'Color', 'w');

        for vi = 1:3
            % MR row
            ax_h = subplot_tight(2, 3, vi, [0.04 0.02]);
            imagesc(mr_slices{vi}); colormap(gray); axis off;
            set(ax_h, 'DataAspectRatioMode', 'auto');

            % CT row — window [0, 100] HU to emphasise soft tissue
            ax_h = subplot_tight(2, 3, 3 + vi, [0.04 0.02]);
            imagesc(ct_slices{vi}, [0 100]); colormap(gray); axis off;
            set(ax_h, 'DataAspectRatioMode', 'auto');
        end

        sgtitle(sprintf('Example paired MR and CT — %s', subjects(mid_idx).name), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
        save_figure(fig0, cfg.fig_dir, 'fig_example_data');
    end

    if cleanup_mr, delete(mr_nii); end
    if cleanup_ct, delete(ct_nii); end
    end  % do(0)

    %% Figure 0b: CTseg template and intensity priors
    if do(0.5)
    dir_ctseg = fileparts(which('spm_CTseg'));
    pth_prior = fullfile(dir_ctseg, 'models', 'prior_CTseg.mat');

    if exist(pth_mu, 'file')
        % Load template (multi-volume: K-1 tissue classes, last is implicit)
        Nii_mu = nifti(pth_mu);
        mu_vol = single(Nii_mu.dat(:,:,:,:));

        % Softmax to get tissue probabilities
        mu_soft = zeros(size(mu_vol,1), size(mu_vol,2), size(mu_vol,3), size(mu_vol,4)+1, 'single');
        for k = 1:size(mu_vol, 4)
            mu_soft(:,:,:,k) = exp(mu_vol(:,:,:,k));
        end
        mu_soft(:,:,:,end) = ones(size(mu_vol,1), size(mu_vol,2), size(mu_vol,3), 'single');
        mu_soft = bsxfun(@rdivide, mu_soft, sum(mu_soft, 4));

        dims = size(mu_soft);
        ax = round(0.50 * dims(3));
        co = round(0.55 * dims(2));
        sa = round(0.50 * dims(1));

        K = size(mu_soft, 4);
        tissue_labels = {'GM', 'WM', 'CSF', 'Bone', 'Soft tissue', 'Background'};
        tissue_colors = [0.8 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; ...
                         0.9 0.9 0.2; 0.9 0.5 0.2; 0.5 0.5 0.5];

        % Show all 6 tissue classes: 3 rows (axial/coronal/sagittal) x 6 columns (tissues)
        n_classes = min(K, 6);
        view_names = {'Axial', 'Coronal', 'Sagittal'};
        fig0b = figure('Name', 'CTseg Template', 'Position', [50 50 1200 550], 'Color', 'k');

        % Build montage: tile all slices into a single image
        gap = 2;  % pixels between tiles
        % Get slice sizes for each view
        slices = cell(3, n_classes);
        for vi = 1:3
            for ki = 1:n_classes
                if vi == 1
                    slices{vi,ki} = rot90(squeeze(mu_soft(:,:,ax,ki)));
                elseif vi == 2
                    slices{vi,ki} = rot90(squeeze(mu_soft(:,co,:,ki)));
                else
                    slices{vi,ki} = rot90(squeeze(mu_soft(sa,:,:,ki)));
                end
            end
        end

        % Compute row heights and column widths (max per row/col)
        row_h = zeros(3, 1); col_w = zeros(1, n_classes);
        for vi = 1:3
            for ki = 1:n_classes
                [h, w] = size(slices{vi,ki});
                row_h(vi) = max(row_h(vi), h);
                col_w(ki) = max(col_w(ki), w);
            end
        end

        % Title row height
        title_h = 20;
        total_h = title_h + sum(row_h) + gap * 2;
        total_w = sum(col_w) + gap * (n_classes - 1);

        montage = zeros(total_h, total_w, 'single');
        y0 = title_h;
        for vi = 1:3
            x0 = 0;
            for ki = 1:n_classes
                [h, w] = size(slices{vi,ki});
                % Centre within cell
                dy = round((row_h(vi) - h) / 2);
                dx = round((col_w(ki) - w) / 2);
                montage(y0+dy+1:y0+dy+h, x0+dx+1:x0+dx+w) = slices{vi,ki};
                x0 = x0 + col_w(ki) + gap;
            end
            y0 = y0 + row_h(vi) + gap;
        end

        fig0b = figure('Name', 'CTseg Template', 'Color', 'k');
        imshow(montage, [0 1], 'Border', 'tight');
        set(fig0b, 'Color', 'k');

        % Add tissue labels as text annotations
        x0 = 0;
        for ki = 1:n_classes
            tx = x0 + col_w(ki) / 2;
            text(tx, title_h * 0.6, tissue_labels{min(ki, numel(tissue_labels))}, ...
                 'Color', 'w', 'FontSize', 9, 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold');
            x0 = x0 + col_w(ki) + gap;
        end

        save_figure(fig0b, cfg.fig_dir, 'fig_ctseg_template', true);
    end

    % Plot original (groupwise) CTseg template for appendix
    pth_mu_orig = fullfile(dir_ctseg, 'models', 'mu_CTseg.nii');
    if exist(pth_mu_orig, 'file')
        Nii_mu_orig = nifti(pth_mu_orig);
        mu_vol_orig = single(Nii_mu_orig.dat(:,:,:,:));

        mu_soft_orig = zeros(size(mu_vol_orig,1), size(mu_vol_orig,2), size(mu_vol_orig,3), size(mu_vol_orig,4)+1, 'single');
        for k = 1:size(mu_vol_orig, 4)
            mu_soft_orig(:,:,:,k) = exp(mu_vol_orig(:,:,:,k));
        end
        mu_soft_orig(:,:,:,end) = ones(size(mu_vol_orig,1), size(mu_vol_orig,2), size(mu_vol_orig,3), 'single');
        mu_soft_orig = bsxfun(@rdivide, mu_soft_orig, sum(mu_soft_orig, 4));

        dims_orig = size(mu_soft_orig);
        ax_orig = round(0.50 * dims_orig(3));
        co_orig = round(0.55 * dims_orig(2));
        sa_orig = round(0.50 * dims_orig(1));

        K_orig = size(mu_soft_orig, 4);
        tissue_labels = {'GM', 'WM', 'CSF', 'Bone', 'Soft tissue', 'Background'};
        n_classes_orig = min(K_orig, 6);

        fig_orig = figure('Name', 'CTseg Original Template', 'Position', [50 50 1200 550], 'Color', 'k');
        gap = 2;
        slices_orig = cell(3, n_classes_orig);
        for vi = 1:3
            for ki = 1:n_classes_orig
                if vi == 1
                    slices_orig{vi,ki} = rot90(squeeze(mu_soft_orig(:,:,ax_orig,ki)));
                elseif vi == 2
                    slices_orig{vi,ki} = rot90(squeeze(mu_soft_orig(:,co_orig,:,ki)));
                else
                    slices_orig{vi,ki} = rot90(squeeze(mu_soft_orig(sa_orig,:,:,ki)));
                end
            end
        end

        row_h_orig = zeros(3, 1); col_w_orig = zeros(1, n_classes_orig);
        for vi = 1:3
            for ki = 1:n_classes_orig
                [h, w] = size(slices_orig{vi,ki});
                row_h_orig(vi) = max(row_h_orig(vi), h);
                col_w_orig(ki) = max(col_w_orig(ki), w);
            end
        end

        title_h = 20;
        total_h_orig = title_h + sum(row_h_orig) + gap * 2;
        total_w_orig = sum(col_w_orig) + gap * (n_classes_orig - 1);

        montage_orig = zeros(total_h_orig, total_w_orig, 'single');
        y0 = title_h;
        for vi = 1:3
            x0 = 0;
            for ki = 1:n_classes_orig
                [h, w] = size(slices_orig{vi,ki});
                dy = round((row_h_orig(vi) - h) / 2);
                dx = round((col_w_orig(ki) - w) / 2);
                montage_orig(y0+dy+1:y0+dy+h, x0+dx+1:x0+dx+w) = slices_orig{vi,ki};
                x0 = x0 + col_w_orig(ki) + gap;
            end
            y0 = y0 + row_h_orig(vi) + gap;
        end

        imshow(montage_orig, [0 1], 'Border', 'tight');
        set(fig_orig, 'Color', 'k');

        x0 = 0;
        for ki = 1:n_classes_orig
            tx = x0 + col_w_orig(ki) / 2;
            text(tx, title_h * 0.6, tissue_labels{min(ki, numel(tissue_labels))}, ...
                 'Color', 'w', 'FontSize', 9, 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold');
            x0 = x0 + col_w_orig(ki) + gap;
        end

        save_figure(fig_orig, cfg.fig_dir, 'fig_ctseg_template_original', true);
    end

    % Plot GMM intensity priors
    if exist(pth_prior, 'file')
        S = load(pth_prior);

        tissue_labels = {'GM', 'WM', 'CSF', 'Bone', 'Soft tissue', 'Background'};
        tissue_colors = [0.8 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; ...
                         0.9 0.9 0.2; 0.9 0.5 0.2; 0.5 0.5 0.5];

        % prior_CTseg.mat format:
        %   pr:    {1×6} cell — {means(1×G), b(1×G), W(1×1×G), n(1×G), ...}
        %   mg_ix: (1×G) — maps each Gaussian to a tissue class (1..K)
        % Means are in HU+1024 space (CT data is stored with 1024 offset).
        pr_m  = S.pr{1};       % (1×G) means
        pr_b  = S.pr{2};       % (1×G) precision scaling
        pr_W  = S.pr{3};       % (1×1×G) Wishart scale
        pr_n  = S.pr{4};       % (1×G) Wishart df
        mg_ix = S.mg_ix;       % (1×G) tissue class index per Gaussian

        G = numel(pr_m);
        fprintf('  Prior: %d Gaussians across %d tissue classes\n', G, max(mg_ix(:)));

        % Convert means from internal (HU+1024) to HU
        hu_offset = 1024;
        mu_hu = pr_m(:)' - hu_offset;
        sd_hu = zeros(1, G);
        for g = 1:G
            Wg = pr_W(1, 1, g);
            sd_hu(g) = sqrt(abs(1 / (pr_n(g) * Wg)));
        end

        for g = 1:G
            fprintf('    Gaussian %2d: class=%d (%s), mu=%.0f HU, sd=%.0f\n', ...
                g, mg_ix(g), tissue_labels{min(mg_ix(g), numel(tissue_labels))}, ...
                mu_hu(g), sd_hu(g));
        end

        % X-axis range from data
        x_lo = min(mu_hu - 4*sd_hu);
        x_hi = max(mu_hu + 4*sd_hu);
        x = linspace(x_lo, x_hi, 1000);

        fig0c = figure('Name', 'CTseg Intensity Priors', ...
                       'Position', [50 50 800 400], 'Color', 'w');
        hold on;
        plotted_classes = [];
        for g = 1:G
            ki = mg_ix(g);
            ci = min(ki, size(tissue_colors, 1));
            y = normpdf(x, mu_hu(g), sd_hu(g));
            if ~ismember(ki, plotted_classes)
                plot(x, y, '-', 'Color', tissue_colors(ci,:), 'LineWidth', 2, ...
                     'DisplayName', tissue_labels{min(ki, numel(tissue_labels))});
                plotted_classes = [plotted_classes, ki]; %#ok<AGROW>
            else
                plot(x, y, '-', 'Color', tissue_colors(ci,:), 'LineWidth', 2, ...
                     'HandleVisibility', 'off');
            end
        end
        hold off;
        xlabel('Intensity (HU)', 'FontSize', 12);
        ylabel('Probability density', 'FontSize', 12);
        title('CTseg learned intensity priors', 'FontSize', 12);
        legend('Location', 'northeast', 'FontSize', 10);
        grid on; box on;
        save_figure(fig0c, cfg.fig_dir, 'fig_ctseg_intensity_priors');
    end
    end  % do(0.5)

    %% Load metrics for subsequent figures
    M = load(fullfile(cfg.out_dir, 'metrics.mat'));
    m = M.metrics;

    %% Figure 1: Segmentation metrics (boxplots)
    if do(1)
    fig1 = figure('Name', 'Segmentation Metrics', 'Position', [50 50 1200 800], ...
                  'Color', 'w');

    metric_fields = {'dice', 'hd95', 'assd'};
    metric_titles = {'Dice Coefficient', '95th Percentile Hausdorff Distance (mm)', ...
                     'Average Symmetric Surface Distance (mm)'};

    for mi = 1:3
        subplot_tight(1, 3, mi, [0.08 0.05]);
        spm_data   = m.([metric_fields{mi} '_spm']);
        ctseg_data = m.([metric_fields{mi} '_ctseg']);

        % Combine for grouped boxplot
        n = size(spm_data, 1);
        data   = [];
        groups = {};
        tissues = {};
        for t = 1:cfg.n_tissues
            data = [data; spm_data(:,t); ctseg_data(:,t)]; %#ok<AGROW>
            groups = [groups; repmat({'SPM'}, n, 1); repmat({'CTseg'}, n, 1)]; %#ok<AGROW>
            tissues = [tissues; repmat(cfg.tissue_names(t), 2*n, 1)]; %#ok<AGROW>
        end

        % Use positions to create grouped layout
        positions = [];
        xtick_pos = [];
        xtick_lab = {};
        colors = [0.4 0.6 0.9; 0.9 0.4 0.3];  % blue=SPM, red=CTseg
        hold on;
        for t = 1:cfg.n_tissues
            idx_spm   = strcmp(tissues, cfg.tissue_names{t}) & strcmp(groups, 'SPM');
            idx_ctseg = strcmp(tissues, cfg.tissue_names{t}) & strcmp(groups, 'CTseg');

            pos_spm   = (t-1)*3 + 1;
            pos_ctseg = (t-1)*3 + 2;

            bp1 = boxplot_single(data(idx_spm),   pos_spm,   colors(1,:));
            bp2 = boxplot_single(data(idx_ctseg), pos_ctseg, colors(2,:));

            xtick_pos = [xtick_pos, pos_spm + 0.5]; %#ok<AGROW>
            xtick_lab = [xtick_lab, cfg.tissue_names(t)]; %#ok<AGROW>
        end
        hold off;

        set(gca, 'XTick', xtick_pos, 'XTickLabel', xtick_lab, 'FontSize', 11);
        ylabel(metric_titles{mi}, 'FontSize', 12);
        title(metric_titles{mi}, 'FontSize', 13);
        grid on; box on;

        % Legend (only on first subplot)
        if mi == 1
            h1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'k');
            h2 = patch(NaN, NaN, colors(2,:), 'EdgeColor', 'k');
            legend([h1 h2], {'SPM', 'CTseg'}, 'Location', 'southwest', 'FontSize', 10);
        end
    end

    % sgtitle('Segmentation metrics', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
    save_figure(fig1, cfg.fig_dir, 'fig_segmentation_metrics');
    end  % do(1)

    %% Figure 1b: Normalisation metrics (step4b)
    if do(1.5)
    if exist(fullfile(cfg.out_dir, 'norm_metrics.mat'), 'file')
        NM = load(fullfile(cfg.out_dir, 'norm_metrics.mat'));
        nm = NM.norm_metrics;

        fig1b = figure('Name', 'Normalisation Metrics', 'Position', [50 50 900 800], ...
                        'Color', 'w');

        % HD95 omitted: at 1.5mm template resolution, bwdist quantises
        % surface distances so HD95 collapses to a few discrete values.
        norm_fields = {'dice', 'assd'};
        norm_titles = {'Dice Coefficient', 'Average Symmetric Surface Distance (mm)'};
        colors = [0.4 0.6 0.9; 0.9 0.4 0.3];  % blue=SPM-CT, red=CTseg

        for mi = 1:numel(norm_fields)
            subplot_tight(1, 2, mi, [0.08 0.05]);
            hold on;
            for t = 1:cfg.n_tissues
                data_spm   = nm.([norm_fields{mi} '_spm'])(:, t);
                data_ctseg = nm.([norm_fields{mi} '_ctseg'])(:, t);

                pos_spm   = (t-1)*3 + 1;
                pos_ctseg = (t-1)*3 + 2;
                boxplot_single(data_spm,   pos_spm,   colors(1,:));
                boxplot_single(data_ctseg, pos_ctseg, colors(2,:));
            end
            hold off;
            xtick_pos = 1.5:3:(cfg.n_tissues*3);
            set(gca, 'XTick', xtick_pos, 'XTickLabel', cfg.tissue_names, 'FontSize', 11);
            ylabel(norm_titles{mi}, 'FontSize', 12);
            title(norm_titles{mi}, 'FontSize', 13);
            grid on; box on;

            if mi == 1
                h1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'k');
                h2 = patch(NaN, NaN, colors(2,:), 'EdgeColor', 'k');
                legend([h1 h2], {'SPM-CT', 'CTseg'}, 'Location', 'southwest', 'FontSize', 10);
            end
        end

        % sgtitle('Normalisation metrics (warped MR tissue maps vs SPM-MR reference)', ...
        %         'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
        save_figure(fig1b, cfg.fig_dir, 'fig_normalisation_metrics');

        % --- Example: warped MR GM for one subject (3 rows x 3 views) ---
        % Pick a subject near the median CTseg normalisation Dice
        mean_norm_dice = mean(nm.dice_ctseg, 2, 'omitnan');
        [~, sorted_idx] = sort(mean_norm_dice);
        mid_idx = sorted_idx(round(numel(sorted_idx)/2));
        subj_name = nm.subject_names{mid_idx};
        subj_dir  = fullfile(cfg.data_dir, subj_name);

        % Find deformation fields and MR GM tissue map
        mr_files = get_tissue_files(subj_dir, 'mr_spm');
        [~, nam_mr] = fileparts(find_nii(subj_dir, 'pp_mr'));
        pth_y_mr = fullfile(subj_dir, ['y_' nam_mr '.nii']);
        ct_nii = find_nii(subj_dir, 'pp_ct_def');
        if isempty(ct_nii), ct_nii = find_nii(subj_dir, 'pp_ct'); end
        [~, nam_ct] = fileparts(ct_nii);
        pth_y_ct = fullfile(subj_dir, ['y_' nam_ct '.nii']);
        d_y = dir(fullfile(subj_dir, 'y_*_CTseg.nii'));

        if exist(mr_files{1}, 'file') && exist(pth_y_mr, 'file') && ...
           exist(pth_y_ct, 'file') && ~isempty(d_y)
            pth_y_ctseg = fullfile(subj_dir, d_y(1).name);

            % Warp GM with all three deformations
            wc_mr    = warp_for_figure(mr_files{1}, pth_y_mr,    pth_mu, subj_dir, 'wnorm_mr_');
            wc_spmct = warp_for_figure(mr_files{1}, pth_y_ct,    pth_mu, subj_dir, 'wnorm_spmct_');
            wc_ctseg = warp_ctseg_for_figure(mr_files{1}, pth_y_ctseg, pth_mu, subj_dir);

            V_mr    = spm_vol(wc_mr);
            V_spmct = spm_vol(wc_spmct);
            V_ctseg = spm_vol(wc_ctseg);
            Y_mr    = spm_read_vols(V_mr);
            Y_spmct = spm_read_vols(V_spmct);
            Y_ctseg = spm_read_vols(V_ctseg);

            % Reslice CTseg to SPM grid if needed
            if ~isequal(V_mr.dim, V_ctseg.dim) || norm(V_mr.mat - V_ctseg.mat) > 1e-4
                Y_ctseg = reslice_to_target(V_ctseg, V_mr);
            end

            % Clean up temporary files
            delete(wc_mr); delete(wc_spmct); delete(wc_ctseg);

            % Pick slices
            dims = size(Y_mr);
            vox_c = round(V_mr.mat \ [0; -20; 20; 1]);
            ax = max(1, min(round(vox_c(3)), dims(3)));
            co = max(1, min(round(vox_c(2)), dims(2)));
            sa = max(1, min(round(vox_c(1)), dims(1)));

            fig1c = figure('Name', 'Normalisation Example', ...
                           'Position', [50 50 900 600], 'Color', 'w');

            view_titles = {'Axial', 'Coronal', 'Sagittal'};
            row_labels  = {'SPM-MR def.', 'SPM-CT def.', 'CTseg def.'};
            mr_slices    = {rot90(squeeze(Y_mr(:,:,ax))), ...
                            rot90(squeeze(Y_mr(:,co,:))), ...
                            rot90(squeeze(Y_mr(sa,:,:)))};
            spmct_slices = {rot90(squeeze(Y_spmct(:,:,ax))), ...
                            rot90(squeeze(Y_spmct(:,co,:))), ...
                            rot90(squeeze(Y_spmct(sa,:,:)))};
            ctseg_slices = {rot90(squeeze(Y_ctseg(:,:,ax))), ...
                            rot90(squeeze(Y_ctseg(:,co,:))), ...
                            rot90(squeeze(Y_ctseg(sa,:,:)))};

            all_slices = {mr_slices, spmct_slices, ctseg_slices};
            for ri = 1:3
                for vi = 1:3
                    subplot_tight(3, 3, (ri-1)*3 + vi, [0.04 0.02]);
                    imagesc(all_slices{ri}{vi}, [0 1]); colormap(gray);
                    axis image; set(gca, 'XTick', [], 'YTick', []);
                    if ri == 1, title(view_titles{vi}, 'FontSize', 11); end
                    if vi == 1
                        ylabel(row_labels{ri}, 'FontSize', 12, 'FontWeight', 'bold');
                    end
                end
            end

            sgtitle(sprintf('Warped MR GM — %s', subj_name), ...
                    'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
            save_figure(fig1c, cfg.fig_dir, 'fig_normalisation_example');
        end
    end
    end  % do(1.5)

    %% Figure 2: Bland-Altman plots (volumetrics)
    if do(2)
    % 2x2 grid: rows = SPM-CT, CTseg; columns = TBV, TIV
    if exist(fullfile(cfg.out_dir, 'volumetrics.mat'), 'file')
        V = load(fullfile(cfg.out_dir, 'volumetrics.mat'));
        v = V.volumetrics;

        fig2 = figure('Name', 'Volumetrics', 'Position', [50 50 1000 700], 'Color', 'w');

        subplot_tight(2, 2, 1, [0.10 0.08]);
        bland_altman_plot(v.tbv_mr, v.tbv_spm, 'TBV (ml)', 'SPM');

        subplot_tight(2, 2, 2, [0.10 0.08]);
        bland_altman_plot(v.tiv_mr, v.tiv_spm, 'TIV (ml)', 'SPM');

        subplot_tight(2, 2, 3, [0.10 0.08]);
        bland_altman_plot(v.tbv_mr, v.tbv_ctseg, 'TBV (ml)', 'CTseg');

        subplot_tight(2, 2, 4, [0.10 0.08]);
        bland_altman_plot(v.tiv_mr, v.tiv_ctseg, 'TIV (ml)', 'CTseg');

        % sgtitle('Bland–Altman: CT vs MR volume estimates', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
        save_figure(fig2, cfg.fig_dir, 'fig_volumetrics_bland_altman');
    end
    end  % do(2)

    %% Figure 3: Average normalised CTs (SPM-MR vs SPM-CT vs CTseg)
    if do(3)
    % 1 row x 3 columns: axial average normalised CT for SPM-MR, SPM-CT, CTseg
    avg_mr_f    = fullfile(cfg.out_dir, 'avg_normalised_ct_mr.nii');
    avg_spm_f   = fullfile(cfg.out_dir, 'avg_normalised_ct_spm.nii');
    avg_ctseg_f = fullfile(cfg.out_dir, 'avg_normalised_ct_ctseg.nii');
    if exist(avg_spm_f, 'file') && exist(avg_ctseg_f, 'file')
        V_spm   = spm_vol(avg_spm_f);
        Y_spm   = spm_read_vols(V_spm);

        V_ctseg = spm_vol(avg_ctseg_f);
        Y_ctseg = reslice_to_target(V_ctseg, V_spm);

        has_mr = exist(avg_mr_f, 'file');
        if has_mr
            V_mr = spm_vol(avg_mr_f);
            Y_mr = reslice_to_target(V_mr, V_spm);
        end

        dims = size(Y_spm);

        % Pick axial slice in world coordinates
        world_coords = [0, -20, 20];
        vox = V_spm.mat \ [world_coords(1); world_coords(2); world_coords(3); 1];
        ax = max(1, min(round(vox(3)), dims(3)));

        clim = [0 100];

        if has_mr
            col_labels = {'SPM-MR', 'SPM-CT', 'CTseg'};
            slices = {rot90(squeeze(Y_mr(:,:,ax))), ...
                      rot90(squeeze(Y_spm(:,:,ax))), ...
                      rot90(squeeze(Y_ctseg(:,:,ax)))};
        else
            col_labels = {'SPM-CT', 'CTseg'};
            slices = {rot90(squeeze(Y_spm(:,:,ax))), ...
                      rot90(squeeze(Y_ctseg(:,:,ax)))};
        end

        n_cols = numel(slices);
        fig3 = figure('Name', 'Average Normalised CTs', ...
                      'Position', [50 50 300*n_cols 350], 'Color', 'w');

        for ci = 1:n_cols
            subplot_tight(1, n_cols, ci, [0.02 0.01]);
            imagesc(slices{ci}, clim); colormap(gray);
            axis image; set(gca, 'XTick', [], 'YTick', []);
            title(col_labels{ci}, 'FontSize', 12, 'FontWeight', 'bold');
        end

        save_figure(fig3, cfg.fig_dir, 'fig_average_normalised');
    end
    end  % do(3)

    %% Figure 4: Best and worst Dice examples
    if do(4)
    % 3 rows (MR ref, SPM-on-CT, CTseg-on-CT) x 4 columns (best/worst SPM, best/worst CTseg)
    if exist(fullfile(cfg.out_dir, 'metrics.mat'), 'file')
        mean_dice_spm   = mean(m.dice_spm, 2, 'omitnan');
        mean_dice_ctseg = mean(m.dice_ctseg, 2, 'omitnan');

        [~, best_spm]    = max(mean_dice_spm);
        [~, worst_spm]   = min(mean_dice_spm);
        [~, best_ctseg]  = max(mean_dice_ctseg);
        [~, worst_ctseg] = min(mean_dice_ctseg);

        col_idx   = [best_spm, worst_spm, best_ctseg, worst_ctseg];
        col_label = {sprintf('SPM best (%.3f)',   mean_dice_spm(best_spm)), ...
                     sprintf('SPM worst (%.3f)',  mean_dice_spm(worst_spm)), ...
                     sprintf('CTseg best (%.3f)', mean_dice_ctseg(best_ctseg)), ...
                     sprintf('CTseg worst (%.3f)',mean_dice_ctseg(worst_ctseg))};

        fig4 = figure('Name', 'Best/Worst Segmentations', ...
                      'Position', [50 50 1200 900], 'Color', 'w');
        colors = [1 0 0; 0 0 1; 0 1 0];  % GM=red, WM=blue, CSF=green
        row_labels = {sprintf('SPM-MR'), sprintf('SPM-CT'), sprintf('CTseg-CT')};
        methods = {'mr_spm', 'ct_spm', 'ct_ctseg'};

        % Preload background slices for each column subject
        bg = cell(4, 1); bg_z = zeros(4, 1); bg_cleanup = cell(4, 1);
        for ci = 1:4
            subj_dir = fullfile(cfg.data_dir, m.subject_names{col_idx(ci)});
            [bg{ci}, bg_z(ci), bg_cleanup{ci}] = load_bg_slice(subj_dir);
        end

        left_margin = 0.05;  % space for row labels
        ax_handles = gobjects(3, 4);
        for ri = 1:3
            for ci = 1:4
                subj_dir = fullfile(cfg.data_dir, m.subject_names{col_idx(ci)});
                tissue_f = get_tissue_files(subj_dir, methods{ri});
                ax_handles(ri,ci) = subplot_tight(3, 4, (ri-1)*4 + ci, [0.01 0.005]);
                % Shift right to make room for row labels
                pos = get(ax_handles(ri,ci), 'Position');
                pos(1) = pos(1) * (1 - left_margin) + left_margin;
                pos(3) = pos(3) * (1 - left_margin);
                set(ax_handles(ri,ci), 'Position', pos);
                if exist(tissue_f{1}, 'file')
                    show_seg_overlay(bg{ci}, tissue_f, bg_z(ci), colors, cfg.dice_thresh);
                else
                    imshow(repmat(bg{ci}, [1 1 3])); axis image off;
                end
                if ri == 1
                    title(sprintf('%s\n%s', col_label{ci}, m.subject_names{col_idx(ci)}), ...
                          'FontSize', 9);
                end
            end
            % Add row label as annotation
            mid_y = pos(2) + pos(4)/2;
            annotation(fig4, 'textbox', [0, mid_y - 0.05, left_margin, 0.10], ...
                       'String', row_labels{ri}, 'FontSize', 10, 'FontWeight', 'bold', ...
                       'Color', [0 0 0], 'EdgeColor', 'none', ...
                       'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                       'FitBoxToText', 'off');
        end

        % Cleanup temp files
        for ci = 1:4
            if ~isempty(bg_cleanup{ci}), delete(bg_cleanup{ci}); end
        end

        save_figure(fig4, cfg.fig_dir, 'fig_best_worst_segmentations');
    end
    end  % do(4)

    %% Figure 5: Example tissue maps for a random subject
    if do(5)
    % 3 rows (GM, WM, CSF) x 3 columns (MR-SPM, SPM-on-CT, CTseg-on-CT)
    if exist(fullfile(cfg.out_dir, 'metrics.mat'), 'file')
        % Pick a subject near the median Dice (representative, not extreme)
        mean_dice_all = mean([m.dice_spm, m.dice_ctseg], 2, 'omitnan');
        [~, sorted_idx] = sort(mean_dice_all);
        med_idx = sorted_idx(round(numel(sorted_idx)/2));

        subj_name = m.subject_names{med_idx};
        subj_dir  = fullfile(cfg.data_dir, subj_name);

        mr_files    = get_tissue_files(subj_dir, 'mr_spm');
        ct_spm_f    = get_tissue_files(subj_dir, 'ct_spm');
        ct_ctseg_f  = get_tissue_files(subj_dir, 'ct_ctseg');

        if exist(mr_files{1}, 'file')
            fig5 = figure('Name', 'Example Tissue Maps', ...
                          'Position', [50 50 1000 900], 'Color', 'w');

            % Get slice index from first tissue map
            V = spm_vol(mr_files{1});
            nz = V.dim(3);
            z = round(0.45 * nz);

            % Get per-tissue Dice for this subject
            dice_spm   = m.dice_spm(med_idx, :);    % [GM, WM, CSF]
            dice_ctseg = m.dice_ctseg(med_idx, :);

            all_files = {mr_files, ct_spm_f, ct_ctseg_f};
            col_titles_base = {'SPM on MR (ref)', 'SPM on CT', 'CTseg on CT'};
            tissue_labels = cfg.tissue_names;  % {'GM', 'WM', 'CSF'}

            for ti = 1:cfg.n_tissues
                for ci = 1:3
                    subplot_tight(cfg.n_tissues, 3, (ti-1)*3 + ci, [0.02 0.005]);
                    f = all_files{ci}{ti};
                    if exist(f, 'file')
                        Y = spm_read_vols(spm_vol(f));
                        slice = rot90(squeeze(Y(:,:,z)));
                        imagesc(slice, [0 1]);
                        colormap(gray); axis image off;
                    else
                        axis off;
                    end
                    if ci == 1
                        t_str = col_titles_base{ci};
                    elseif ci == 2
                        t_str = sprintf('%s (D=%.2f)', col_titles_base{ci}, dice_spm(ti));
                    else
                        t_str = sprintf('%s (D=%.2f)', col_titles_base{ci}, dice_ctseg(ti));
                    end
                    title(t_str, 'FontSize', 10);
                    if ci == 1, ylabel(tissue_labels{ti}, 'FontSize', 12); end
                end
            end

            save_figure(fig5, cfg.fig_dir, 'fig_example_tissue_maps');
        end
    end
    end  % do(5)

    %% Figure 6: Demographics distributions
    if do(6)
    % Show demographics only for the subjects used in prediction (after filtering)
    pred_file = fullfile(cfg.out_dir, 'prediction.mat');
    if exist(pred_file, 'file')
        P = load(pred_file);
        pp = P.prediction;
        % Use first method's data (all methods use the same subjects)
        pred_ages = pp.mr_spm.age_true;
        pred_sexes = pp.mr_spm.sex_true;  % 0=female, 1=male
    else
        % Fallback to full demographics if prediction not yet run
        demo = load_demographics(cfg);
        pred_ages = demo.age(~isnan(demo.age));
        pred_sexes = demo.sex(~isnan(demo.sex));
    end

    n_pred = numel(pred_ages);
    fig6 = figure('Name', 'Demographics', 'Position', [50 50 900 350], 'Color', 'w');

    subplot_tight(1, 2, 1, [0.10 0.08]);
    histogram(pred_ages, 15, 'FaceColor', [0.4 0.6 0.9], 'FaceAlpha', 1, 'EdgeColor', 'k');
    xlabel('Age (years)', 'FontSize', 12);
    ylabel('Count', 'FontSize', 12);
    title(sprintf('Age distribution (n=%d, mean=%.0f, std=%.0f)', ...
        numel(pred_ages), mean(pred_ages), std(pred_ages)), 'FontSize', 11);
    grid on; box on;

    subplot_tight(1, 2, 2, [0.10 0.08]);
    n_male = sum(pred_sexes == 1);
    n_female = sum(pred_sexes == 0);
    sex_labels = {'Male', 'Female'};
    sex_counts = [n_male, n_female];
    bar(sex_counts, 'FaceColor', [0.4 0.6 0.9], 'EdgeColor', 'k');
    set(gca, 'XTickLabel', sex_labels, 'FontSize', 11);
    ylabel('Count', 'FontSize', 12);
    title(sprintf('Sex distribution (n=%d)', n_male + n_female), 'FontSize', 11);
    grid on; box on;

    save_figure(fig6, cfg.fig_dir, 'fig_demographics');
    end  % do(6)

    %% Figure 7: Prediction results
    if do(7)
    pred_file = fullfile(cfg.out_dir, 'prediction.mat');
    if exist(pred_file, 'file')
        P = load(pred_file);
        p = P.prediction;

        methods = {'mr_spm', 'ct_spm', 'ct_ctseg'};
        method_labels = p.methods;
        colors = [0.2 0.6 0.2; 0.4 0.6 0.9; 0.9 0.4 0.3];  % green=MR, blue=SPM-CT, red=CTseg

        % --- Age: scatter plots (true vs predicted) ---
        has_age_data = isfield(p, 'mr_spm') && isfield(p.mr_spm, 'age_true');
        if has_age_data
            fig7a = figure('Name', 'Age Prediction', 'Position', [50 50 1200 350], 'Color', 'w');
            for mi = 1:3
                subplot_tight(1, 3, mi, [0.10 0.05]);
                s = p.(methods{mi});
                scatter(s.age_true, s.age_pred, 30, colors(mi,:), 'filled', 'MarkerEdgeColor', 'k');
                hold on;
                % Identity line
                lims = [min([s.age_true; s.age_pred]) - 5, max([s.age_true; s.age_pred]) + 5];
                plot(lims, lims, 'k--', 'LineWidth', 1);
                hold off;
                xlim(lims); ylim(lims);
                xlabel('True age (years)', 'FontSize', 11);
                ylabel('Predicted age (years)', 'FontSize', 11);
                title(sprintf('%s\nRMSE=%.1f, r=%.3f', method_labels{mi}, s.age_rmse, s.age_r), ...
                      'FontSize', 11);
                axis square; grid on; box on;
            end
            save_figure(fig7a, cfg.fig_dir, 'fig_age_prediction');
        end

        % --- Sex: ROC curves (all methods on one plot) ---
        % Uses func_val (continuous GP posteriors) from PRoNTo Classification.mat,
        % following the same logic as prt_plot_ROC.m
        has_sex_data = isfield(p, 'mr_spm') && isfield(p.mr_spm, 'cls_file');
        if has_sex_data
            fig7b = figure('Name', 'Sex Classification', 'Position', [50 50 500 450], 'Color', 'w');
            hold on;
            legend_entries = {};
            for mi = 1:3
                s = p.(methods{mi});
                % Load PRoNTo Classification.mat and extract func_val
                C = load(s.cls_file);
                PRT = C.ResClass{1}.PRT;
                [fpr, tpr, auc] = compute_roc_pronto(PRT, 1);
                plot(fpr, tpr, '-', 'Color', colors(mi,:), 'LineWidth', 2);
                legend_entries{end+1} = sprintf('%s (AUC=%.3f)', method_labels{mi}, auc); %#ok<AGROW>
            end
            plot([0 1], [0 1], 'k--', 'LineWidth', 1);  % chance line
            hold off;
            xlabel('False positive rate', 'FontSize', 12);
            ylabel('True positive rate', 'FontSize', 12);
            title('Sex classification ROC', 'FontSize', 12);
            legend(legend_entries, 'Location', 'southeast', 'FontSize', 10);
            axis square; grid on; box on;
            save_figure(fig7b, cfg.fig_dir, 'fig_sex_roc');
        end

        % --- Debug: example mwc input slices ---
        fig7c = figure('Name', 'Prediction Input (mwc)', 'Position', [50 50 1000 700], 'Color', 'w');
        % Show GM mwc for one subject across all 3 methods
        subjects = get_subjects(cfg);
        subj_dir = fullfile(cfg.data_dir, subjects(1).name);
        method_ids = {'mr_spm', 'ct_spm', 'ct_ctseg'};
        view_names = {'Axial', 'Coronal', 'Sagittal'};

        for mi = 1:3
            mwc_f = get_mwc_files(subj_dir, method_ids{mi});
            if ~exist(mwc_f{1}, 'file'), continue; end
            vol = spm_read_vols(spm_vol(mwc_f{1}));  % GM only
            dims = size(vol);
            ax = round(0.45 * dims(3));
            co = round(0.55 * dims(2));
            sa = round(0.50 * dims(1));
            slices = {rot90(squeeze(vol(:,:,ax))), ...
                      rot90(squeeze(vol(:,co,:))), ...
                      rot90(squeeze(vol(sa,:,:)))};
            for vi = 1:3
                subplot_tight(3, 3, (mi-1)*3 + vi, [0.03 0.02]);
                imagesc(slices{vi}, [0 0.5]); colormap(gray); axis image off;
                if mi == 1, title(view_names{vi}, 'FontSize', 11); end
                if vi == 1, ylabel(method_labels{mi}, 'FontSize', 11); end
            end
        end
        sgtitle(sprintf('GM modulated warped tissue maps — %s', subjects(1).name), 'FontSize', 13, 'FontWeight', 'bold', 'Color', 'k');
        save_figure(fig7c, cfg.fig_dir, 'fig_prediction_input_debug');
    end

    end  % do(7)

    fprintf('  Figures saved to %s\n', cfg.fig_dir);
end


%% ---- Local helpers ----

% compute_roc_pronto is in experiments/utils/compute_roc_pronto.m

function [ct_disp, z, cleanup_file] = load_bg_slice(subj_dir)
% Load a CT background slice for overlay figures.
% Returns normalised [0,1] slice, slice index, and file to clean up (or '').
    cleanup_file = '';
    ct_def = fullfile(subj_dir, 'pp_ct_def.nii');
    ct_nii = fullfile(subj_dir, 'pp_ct.nii');
    ct_gz  = fullfile(subj_dir, 'pp_ct.nii.gz');
    if exist(ct_def, 'file')
        ct_nii = ct_def;
    elseif ~exist(ct_nii, 'file') && exist(ct_gz, 'file')
        gunzip(ct_gz, subj_dir);
        cleanup_file = ct_nii;
    end
    if exist(ct_nii, 'file')
        ct_vol = spm_read_vols(spm_vol(ct_nii));
    else
        ct_vol = zeros(192, 192, 192);
    end
    nz = size(ct_vol, 3);
    z  = round(0.45 * nz);
    ct_slice = rot90(squeeze(ct_vol(:,:,z)));
    ct_disp  = max(0, min(1, (ct_slice + 50) / 150));  % window [-50, 100] HU
end

function show_seg_overlay(bg_slice, tissue_files, z, colors, thresh)
% Overlay binarised tissue maps on a grayscale background slice.
    rgb = repmat(bg_slice, [1 1 3]);
    alpha = 0.35;
    for t = 1:min(3, numel(tissue_files))
        if ~exist(tissue_files{t}, 'file'), continue; end
        vol = spm_read_vols(spm_vol(tissue_files{t}));
        seg = rot90(squeeze(vol(:,:,z))) > thresh;
        for ch = 1:3
            plane = rgb(:,:,ch);
            plane(seg) = (1-alpha)*plane(seg) + alpha*colors(t, ch);
            rgb(:,:,ch) = plane;
        end
    end
    imshow(rgb); axis image off;
end

function bp = boxplot_single(data, pos, col)
% Draw a single boxplot at a given position with a given colour.
    data = data(~isnan(data));
    if isempty(data), bp = []; return; end

    q = prctile(data, [25 50 75]);
    iqr_val = q(3) - q(1);
    lo = max(min(data), q(1) - 1.5*iqr_val);
    hi = min(max(data), q(3) + 1.5*iqr_val);
    w = 0.35;

    % Box
    bp = patch([pos-w pos+w pos+w pos-w], [q(1) q(1) q(3) q(3)], col, ...
               'EdgeColor', 'k', 'LineWidth', 1);
    % Median
    plot([pos-w pos+w], [q(2) q(2)], 'k-', 'LineWidth', 2);
    % Whiskers
    plot([pos pos], [lo q(1)], 'k-', 'LineWidth', 1);
    plot([pos pos], [q(3) hi], 'k-', 'LineWidth', 1);
    plot([pos-w/2 pos+w/2], [lo lo], 'k-', 'LineWidth', 1);
    plot([pos-w/2 pos+w/2], [hi hi], 'k-', 'LineWidth', 1);

    % Outliers
    outliers = data(data < lo | data > hi);
    if ~isempty(outliers)
        plot(pos*ones(size(outliers)), outliers, 'k.', 'MarkerSize', 6);
    end
end

function bland_altman_plot(x, y, label, method_name)
% Bland-Altman plot: x = reference (MR-SPM), y = test method.
%   method_name: string identifying the test method (e.g. 'SPM', 'CTseg').
    if nargin < 4, method_name = 'Test'; end
    ok = ~isnan(x) & ~isnan(y);
    x = x(ok); y = y(ok);
    mn   = (x + y) / 2;
    diff = y - x;
    mu   = mean(diff);
    sd   = std(diff);

    plot(mn, diff, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', [0.3 0.3 0.3]);
    hold on;
    xl = xlim;
    plot(xl, [mu mu], 'b-', 'LineWidth', 1.5);
    plot(xl, [mu+1.96*sd mu+1.96*sd], 'r--', 'LineWidth', 1);
    plot(xl, [mu-1.96*sd mu-1.96*sd], 'r--', 'LineWidth', 1);
    hold off;
    xlabel(['Mean ' label], 'FontSize', 11);
    ylabel(sprintf('Difference (%s - MR) %s', method_name, label), 'FontSize', 11);
    title(sprintf('%s %s  bias=%.1f  LoA=[%.1f, %.1f]', method_name, label, mu, mu-1.96*sd, mu+1.96*sd), ...
          'FontSize', 12);
    grid on; box on;
end

function Y_out = reslice_to_target(V_src, V_tgt)
% Reslice source volume to match target volume grid using trilinear interp.
%   V_src: spm_vol struct of source image
%   V_tgt: spm_vol struct of target image (defines output grid)
%   Y_out: resliced volume with dimensions of V_tgt
    Y_src = spm_read_vols(V_src);
    dims = V_tgt.dim(1:3);

    % Mapping from target voxels to source voxels
    M = V_src.mat \ V_tgt.mat;
    [x1, x2, x3] = ndgrid(single(1:dims(1)), single(1:dims(2)), single(1:dims(3)));
    y1 = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
    y2 = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
    y3 = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);

    Y_out = spm_bsplins(double(Y_src), double(y1), double(y2), double(y3), [1 1 1 0 0 0]);
    Y_out(~isfinite(Y_out)) = 0;
end

function pth_out = warp_for_figure(tissue_file, def_file, pth_mu_spm, subj_dir, prefix)
% Warp a native tissue map to SPM space using SPM deformation (pull).
% SPM deformations are on the template grid → pull directly.
    if nargin < 5, prefix = 'wnorm_spm_'; end
    matlabbatch = {};
    matlabbatch{1}.spm.util.defs.comp{1}.comp{1}.def         = {def_file};
    matlabbatch{1}.spm.util.defs.comp{1}.space               = {pth_mu_spm};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {tissue_file};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {subj_dir};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = prefix;
    spm_jobman('run', matlabbatch);
    [~, nam, ext] = fileparts(tissue_file);
    pth_out = fullfile(subj_dir, [prefix nam ext]);
end

function pth_out = warp_ctseg_for_figure(tissue_file, def_file, pth_mu_spm, subj_dir)
% Warp a native tissue map to SPM space using CTseg deformation (push).
% CTseg deformations are on the native grid → push to template.
    matlabbatch = {};
    matlabbatch{1}.spm.util.defs.comp{1}.def                 = {def_file};
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames          = {tissue_file};
    matlabbatch{1}.spm.util.defs.out{1}.push.weight          = {''};
    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {subj_dir};
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.file        = {pth_mu_spm};
    matlabbatch{1}.spm.util.defs.out{1}.push.prefix          = 'wnorm_ctseg_';
    spm_jobman('run', matlabbatch);
    [~, nam, ext] = fileparts(tissue_file);
    pth_out = fullfile(subj_dir, ['wnorm_ctseg_' nam ext]);
end

function save_figure(fig, out_dir, name, dark)
% Save figure as both PNG and PDF.
% By default, force white background. If dark=true, preserve figure colors.
    if nargin < 4, dark = false; end
    if ~dark
        set(fig, 'Color', 'w');
        all_ax = findall(fig, 'Type', 'axes');
        for i = 1:numel(all_ax)
            if strcmp(get(all_ax(i), 'Visible'), 'on')
                set(all_ax(i), 'Color', 'w', ...
                    'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
            end
            set(get(all_ax(i), 'Title'), 'Color', 'k');
        end
        all_text = findall(fig, 'Type', 'text');
        for i = 1:numel(all_text)
            % Skip text with 'overlay' or 'rowlabel' tag (intentionally colored)
            tag = get(all_text(i), 'Tag');
            if any(strcmp(tag, {'overlay', 'rowlabel'})), continue; end
            set(all_text(i), 'Color', 'k');
        end
        set(findall(fig, 'Type', 'legend'), 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
        bg = 'white';
    else
        bg = 'none';
    end
    exportgraphics(fig, fullfile(out_dir, [name '.png']), 'Resolution', 300, 'BackgroundColor', bg);
    exportgraphics(fig, fullfile(out_dir, [name '.pdf']), 'BackgroundColor', bg, 'ContentType', 'vector');
end

function h = subplot_tight(m, n, p, margins)
% subplot_tight — subplot with reduced margins.
%   margins = [vgap hgap] as fraction of figure size.
%   Reserves space at the top for sgtitle.
    if nargin < 4, margins = [0.04 0.04]; end
    vgap = margins(1);
    hgap = margins(2);
    top_margin = 0.06;  % reserve for sgtitle
    bot_margin = 0.02;

    % Convert linear index to row, col
    row = ceil(p / n);
    col = mod(p - 1, n) + 1;

    % Compute position
    usable_h = 1 - top_margin - bot_margin;
    usable_w = 1;
    plot_w = (usable_w - (n+1)*hgap) / n;
    plot_h = (usable_h - (m+1)*vgap) / m;
    x = hgap + (col - 1) * (plot_w + hgap);
    y = bot_margin + usable_h - vgap - row * plot_h - (row - 1) * vgap;

    h = axes('Position', [x, y, plot_w, plot_h]);
end
