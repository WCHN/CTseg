function review_segs(mode)
% Interactive viewer to review CTseg segmentation results.
% Click through subjects to inspect preprocessed images and tissue
% segmentations from all three methods (SPM-MR, SPM-CT, CTseg).
%
% Usage:
%   >> review_segs          % show segmentations (default)
%   >> review_segs('segs')  % show segmentations
%   >> review_segs('images') % show images only

    exp_dir = fileparts(fileparts(mfilename('fullpath')));  % experiments/
    addpath(exp_dir);
    addpath(fullfile(exp_dir, 'utils'));
    cfg = config();
    subjects = get_subjects(cfg);
    n = numel(subjects);

    % Load pre-computed metrics (from step4)
    metrics_file = fullfile(cfg.out_dir, 'metrics.mat');
    if exist(metrics_file, 'file')
        tmp = load(metrics_file, 'metrics');
        metrics = tmp.metrics;
    else
        metrics = [];
        warning('metrics.mat not found — run step4 first to see Dice scores.');
    end

    % Sort subjects by mean CTseg Dice (worst first) for easy failure inspection
    if ~isempty(metrics)
        subj_names = {subjects.name};
        mean_dice = nan(numel(subj_names), 1);
        for si = 1:numel(subj_names)
            mi = find(strcmp(metrics.subject_names, subj_names{si}), 1);
            if ~isempty(mi)
                mean_dice(si) = mean(metrics.dice_ctseg(mi,:));
            end
        end
        [~, sort_idx] = sort(mean_dice, 'ascend', 'MissingPlacement', 'last');
        subjects = subjects(sort_idx);
    end

    % State
    idx = 1;
    view_dim = 3;       % 1=sagittal, 2=coronal, 3=axial
    view_labels = {'Sagittal', 'Coronal', 'Axial'};
    if nargin < 1, mode = 'segs'; end
    show_segs = strcmpi(mode, 'segs');
    if show_segs && isempty(metrics)
        error('No metrics.mat found in %s. Run steps 1-4 first, or use viewer(''images'').', cfg.out_dir);
    end

    % Create figure
    fig = figure('Name', 'CTseg Viewer', 'NumberTitle', 'off', ...
                 'Units', 'normalized', 'Position', [0.05 0.05 0.85 0.85], ...
                 'Color', 'k', 'KeyPressFcn', @on_key);

    % Navigation buttons
    uicontrol(fig, 'Style', 'pushbutton', 'String', '< Prev', ...
              'Units', 'normalized', 'Position', [0.35 0.01 0.1 0.03], ...
              'Callback', @(~,~) go(-1));
    uicontrol(fig, 'Style', 'pushbutton', 'String', 'Next >', ...
              'Units', 'normalized', 'Position', [0.55 0.01 0.1 0.03], ...
              'Callback', @(~,~) go(1));
    hview = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Axial', ...
              'Units', 'normalized', 'Position', [0.46 0.01 0.08 0.03], ...
              'Callback', @(~,~) cycle_view());
    hmode = uicontrol(fig, 'Style', 'pushbutton', 'String', iff(show_segs, 'Images only', 'Segs'), ...
              'Units', 'normalized', 'Position', [0.15 0.01 0.1 0.03], ...
              'Callback', @(~,~) toggle_mode());

    % Title
    htitle = uicontrol(fig, 'Style', 'text', 'String', '', ...
                        'Units', 'normalized', 'Position', [0.2 0.96 0.6 0.03], ...
                        'FontSize', 14, 'FontWeight', 'bold', ...
                        'ForegroundColor', 'w', 'BackgroundColor', 'k', ...
                        'HorizontalAlignment', 'center');

    show_subject();

    function go(delta)
        idx = mod(idx - 1 + delta, n) + 1;
        show_subject();
    end

    function toggle_mode()
        show_segs = ~show_segs;
        if show_segs
            set(hmode, 'String', 'Images only');
        else
            set(hmode, 'String', 'Segs');
        end
        show_subject();
    end

    function cycle_view()
        view_dim = mod(view_dim, 3) + 1;  % 3->1->2->3
        set(hview, 'String', view_labels{view_dim});
        show_subject();
    end

    function on_key(~, evt)
        switch evt.Key
            case {'rightarrow', 'space'}
                go(1);
            case 'leftarrow'
                go(-1);
            case 'v'
                cycle_view();
            case 's'
                toggle_mode();
            case {'q', 'escape'}
                close(fig);
        end
    end

    function show_subject()
        subj = subjects(idx).name;
        subj_dir = fullfile(cfg.data_dir, subj);
        set(htitle, 'String', sprintf('%s  (%d / %d)', subj, idx, n));

        % Clear all existing axes
        delete(findobj(fig, 'Type', 'axes'));

        % Load preprocessed images (prefer deformed CT)
        pp_mr = load_nii(fullfile(subj_dir, 'pp_mr.nii.gz'));
        ct_def = fullfile(subj_dir, 'pp_ct_def.nii');
        if exist(ct_def, 'file')
            pp_ct = load_nii(ct_def);
        else
            pp_ct = load_nii(fullfile(subj_dir, 'pp_ct.nii.gz'));
        end

        sl_mr = round(size(pp_mr, view_dim) / 2);
        sl_ct = round(size(pp_ct, view_dim) / 2);

        if ~show_segs
            % --- Images only: MR and CT side by side ---
            ax1 = subplot(1, 2, 1, 'Parent', fig);
            show_slice(ax1, pp_mr, sl_mr, view_dim, 'gray', []);
            title(ax1, 'pp\_mr', 'Color', 'w', 'FontSize', 13);
            set(ax1, 'XTick', [], 'YTick', [], 'Color', 'k');
            axis(ax1, 'image');

            ax2 = subplot(1, 2, 2, 'Parent', fig);
            show_slice(ax2, pp_ct, sl_ct, view_dim, 'gray', [-10 110]);
            title(ax2, 'pp\_ct', 'Color', 'w', 'FontSize', 13);
            set(ax2, 'XTick', [], 'YTick', [], 'Color', 'k');
            axis(ax2, 'image');
        else
            % --- Full view: images + segmentations ---

            % Load SPM-MR segmentations (c1-3)
            mr = cell(1,3);
            for k = 1:3
                mr{k} = load_nii(fullfile(subj_dir, sprintf('c%dpp_mr.nii', k)));
            end

            % Load SPM-CT segmentations (c1-3)
            ct = cell(1,3);
            for k = 1:3
                ct{k} = load_nii(fullfile(subj_dir, sprintf('c%dpp_ct.nii', k)));
            end

            % Load CTseg segmentations (dynamic naming: c0K_*_CTseg.nii)
            ctseg = cell(1,3);
            for k = 1:3
                pat = sprintf('c%02d_*_CTseg.nii', k);
                f = dir(fullfile(subj_dir, pat));
                if ~isempty(f)
                    ctseg{k} = load_nii(fullfile(subj_dir, f(1).name));
                else
                    ctseg{k} = [];
                end
            end

            labels_row = {'', 'GM', 'WM', 'CSF'};
            labels_col = {'SPM-MR', 'SPM-CT', 'CTseg'};

            for r = 1:4
                for c = 1:3
                    ax = subplot(4, 3, (r-1)*3 + c, 'Parent', fig);
                    cla(ax); hold(ax, 'on');

                    if r == 1
                        if c == 1
                            show_slice(ax, pp_mr, sl_mr, view_dim, 'gray', []);
                            title(ax, 'pp\_mr', 'Color', 'w', 'FontSize', 11);
                        elseif c == 2
                            show_slice(ax, pp_ct, sl_ct, view_dim, 'gray', [-10 110]);
                            title(ax, 'pp\_ct', 'Color', 'w', 'FontSize', 11);
                        else
                            cla(ax);
                            set(ax, 'Color', 'k', 'XTick', [], 'YTick', []);
                            dice_str = get_dice_str(metrics, subj, cfg.tissue_names);
                            text(ax, 0.5, 0.5, dice_str, 'Color', 'w', ...
                                 'HorizontalAlignment', 'center', ...
                                 'VerticalAlignment', 'middle', ...
                                 'Units', 'normalized', 'FontSize', 10, ...
                                 'FontName', 'FixedWidth', 'Interpreter', 'tex');
                            title(ax, 'Dice (vs MR)', 'Color', 'w', 'FontSize', 11);
                        end
                    else
                        k = r - 1;
                        if c == 1
                            vol = mr{k};  sl = sl_mr;
                        elseif c == 2
                            vol = ct{k};  sl = sl_ct;
                        else
                            vol = ctseg{k};  sl = sl_ct;
                        end

                        if ~isempty(vol)
                            show_slice(ax, vol, round(size(vol,view_dim)/2), view_dim, 'gray', []);
                        else
                            set(ax, 'Color', 'k');
                            text(ax, 0.5, 0.5, 'missing', 'Color', 'r', ...
                                 'HorizontalAlignment', 'center', 'Units', 'normalized');
                        end

                        if c == 1
                            ylabel(ax, labels_row{r}, 'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold');
                        end
                        if r == 2
                            title(ax, labels_col{c}, 'Color', 'w', 'FontSize', 11);
                        end
                    end

                    set(ax, 'XTick', [], 'YTick', [], 'Color', 'k');
                    axis(ax, 'image');
                end
            end
        end

        drawnow;
    end
end

function vol = load_nii(fname)
% Load a NIfTI volume. Handles .nii and .nii.gz.
    if ~exist(fname, 'file')
        % Try without .gz
        fname2 = regexprep(fname, '\.gz$', '');
        if exist(fname2, 'file')
            fname = fname2;
        else
            vol = [];
            return;
        end
    end
    V = spm_vol(fname);
    vol = double(spm_read_vols(V));
end

function s = get_dice_str(metrics, subj, tissue_names)
% Format Dice scores for the current subject from pre-computed metrics.
    if isempty(metrics)
        s = 'No metrics.mat';
        return;
    end
    mi = find(strcmp(metrics.subject_names, subj), 1);
    if isempty(mi)
        s = 'Not in metrics';
        return;
    end
    lines = cell(1, numel(tissue_names) + 1);
    lines{1} = sprintf('%8s  SPM-CT  CTseg', '');
    for k = 1:numel(tissue_names)
        d_spm   = metrics.dice_spm(mi,k);
        d_ctseg = metrics.dice_ctseg(mi,k);
        if d_spm > d_ctseg
            s_spm   = sprintf('{\\bf%.3f}', d_spm);
            s_ctseg = sprintf('%.3f', d_ctseg);
        elseif d_ctseg > d_spm
            s_spm   = sprintf('%.3f', d_spm);
            s_ctseg = sprintf('{\\bf%.3f}', d_ctseg);
        else
            s_spm   = sprintf('%.3f', d_spm);
            s_ctseg = sprintf('%.3f', d_ctseg);
        end
        lines{k+1} = sprintf('%8s  %s   %s', tissue_names{k}, s_spm, s_ctseg);
    end
    s = strjoin(lines, newline);
end

function show_slice(ax, vol, sl, dim, cmap, clim)
% Display a single slice along the given dimension.
%   dim: 1=sagittal, 2=coronal, 3=axial
    if isempty(vol), return; end
    sl = max(1, min(sl, size(vol, dim)));
    switch dim
        case 1
            img = rot90(squeeze(vol(sl,:,:)), 3);
        case 2
            img = rot90(squeeze(vol(:,sl,:)), 3);
        case 3
            img = rot90(vol(:,:,sl), 3);
    end
    imagesc(ax, img);
    colormap(ax, cmap);
    if ~isempty(clim)
        caxis(ax, clim);
    end
    axis(ax, 'image', 'off');
end

function out = iff(cond, a, b)
    if cond, out = a; else, out = b; end
end
