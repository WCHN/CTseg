function set_origin()
% Interactively set the image origin for MR/CT pairs.
% Displays three orthogonal views of the MR. Click to place the crosshair
% at the desired origin (e.g. anterior commissure or centre of brain).
% Press Enter to confirm — the affine matrix of both MR and CT is updated
% so the clicked point becomes world coordinate [0,0,0].
% Press 'n' to skip a subject, 'q' to quit.
%
% Usage:
%   >> set_origin

    exp_dir = fileparts(fileparts(mfilename('fullpath')));  % experiments/
    addpath(exp_dir);
    addpath(fullfile(exp_dir, 'utils'));

    cfg = config();
    subjects = get_subjects(cfg);
    n = numel(subjects);

    fprintf('Set origin for %d subjects. Click to place origin, Enter to confirm, n to skip, q to quit.\n', n);

    for i = 1:n
        subj = subjects(i).name;
        subj_dir = fullfile(cfg.data_dir, subj);

        % Find MR and CT
        mr_file = find_raw_nii(subj_dir, 'mr');
        ct_file = find_raw_nii(subj_dir, 'ct');
        if isempty(mr_file)
            fprintf('  [%d/%d] %s - no mr.nii, SKIP\n', i, n, subj);
            continue;
        end

        fprintf('  [%d/%d] %s\n', i, n, subj);

        % Load MR volume — use SPM for the affine (consistent with write)
        V = spm_vol(mr_file);
        vol = double(spm_read_vols(V));
        M = V.mat;

        % Initial crosshair at volume centre
        dims = size(vol);
        pos = round(dims / 2);

        % Create figure
        fig = figure('Name', sprintf('Set Origin: %s (%d/%d)', subj, i, n), ...
                     'NumberTitle', 'off', 'Color', 'k', ...
                     'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

        done = false;
        action = '';

        while ~done
            draw_views(fig, vol, pos, subj, i, n, M);

            try
                [x, y, btn] = ginput(1);
            catch
                % Figure was closed
                return;
            end

            if isempty(btn)
                % Enter pressed — confirm
                action = 'confirm';
                done = true;
            elseif btn == 'n'
                action = 'skip';
                done = true;
            elseif btn == 'q'
                close(fig);
                fprintf('  Quit.\n');
                return;
            else
                % Mouse click — update crosshair based on which subplot
                ax = gca;
                ax_pos = get(ax, 'Position');
                % Determine which view was clicked
                tag = get(ax, 'Tag');
                x = round(x);
                y = round(y);
                switch tag
                    case 'axial'
                        pos(1) = clamp(x, 1, dims(1));
                        pos(2) = clamp(dims(2) - y + 1, 1, dims(2));
                    case 'coronal'
                        pos(1) = clamp(x, 1, dims(1));
                        pos(3) = clamp(dims(3) - y + 1, 1, dims(3));
                    case 'sagittal'
                        pos(2) = clamp(x, 1, dims(2));
                        pos(3) = clamp(dims(3) - y + 1, 1, dims(3));
                end
            end
        end

        close(fig);

        if strcmp(action, 'skip')
            fprintf('    Skipped.\n');
            continue;
        end

        % Compute world coordinate of selected voxel
        % SPM uses 1-based voxel indexing
        world = M * [pos(1); pos(2); pos(3); 1];

        % Translation to shift this point to [0,0,0]
        T = eye(4);
        T(1:3, 4) = -world(1:3);

        % Update MR affine
        M_new = T * M;
        update_affine(mr_file, M_new);
        fprintf('    MR origin set at voxel [%d %d %d], world [%.1f %.1f %.1f] -> [0 0 0]\n', ...
                pos, world(1:3));

        % Update CT affine with same translation
        if ~isempty(ct_file)
            V_ct = spm_vol(ct_file);
            M_ct_new = T * V_ct.mat;
            update_affine(ct_file, M_ct_new);
            fprintf('    CT origin updated with same translation.\n');
        end
    end

    fprintf('Done.\n');
end

function draw_views(fig, vol, pos, subj, idx, n, M)
    dims = size(vol);
    world = M * [pos(1); pos(2); pos(3); 1];

    % Auto-window: use percentiles of non-zero voxels
    nz = vol(vol ~= 0);
    if isempty(nz)
        clim = [0 1];
    else
        clim = [prctile(nz, 1) prctile(nz, 99)];
    end

    % Axial slice (view from above)
    ax1 = subplot(2, 2, 1, 'Parent', fig);
    img = flipud(squeeze(vol(:,:,pos(3)))');
    imagesc(ax1, img, clim);
    colormap(ax1, 'gray');
    hold(ax1, 'on');
    % Crosshair
    xline(ax1, pos(1), 'r', 'LineWidth', 0.5);
    yline(ax1, dims(2) - pos(2) + 1, 'r', 'LineWidth', 0.5);
    hold(ax1, 'off');
    title(ax1, sprintf('Axial (z=%d)', pos(3)), 'Color', 'w');
    set(ax1, 'Tag', 'axial', 'XTick', [], 'YTick', [], 'Color', 'k');
    axis(ax1, 'image');

    % Coronal slice
    ax2 = subplot(2, 2, 2, 'Parent', fig);
    img = flipud(squeeze(vol(:,pos(2),:))');
    imagesc(ax2, img, clim);
    colormap(ax2, 'gray');
    hold(ax2, 'on');
    xline(ax2, pos(1), 'r', 'LineWidth', 0.5);
    yline(ax2, dims(3) - pos(3) + 1, 'r', 'LineWidth', 0.5);
    hold(ax2, 'off');
    title(ax2, sprintf('Coronal (y=%d)', pos(2)), 'Color', 'w');
    set(ax2, 'Tag', 'coronal', 'XTick', [], 'YTick', [], 'Color', 'k');
    axis(ax2, 'image');

    % Sagittal slice
    ax3 = subplot(2, 2, 3, 'Parent', fig);
    img = flipud(squeeze(vol(pos(1),:,:))');
    imagesc(ax3, img, clim);
    colormap(ax3, 'gray');
    hold(ax3, 'on');
    xline(ax3, pos(2), 'r', 'LineWidth', 0.5);
    yline(ax3, dims(3) - pos(3) + 1, 'r', 'LineWidth', 0.5);
    hold(ax3, 'off');
    title(ax3, sprintf('Sagittal (x=%d)', pos(1)), 'Color', 'w');
    set(ax3, 'Tag', 'sagittal', 'XTick', [], 'YTick', [], 'Color', 'k');
    axis(ax3, 'image');

    % Info panel
    ax4 = subplot(2, 2, 4, 'Parent', fig);
    cla(ax4);
    set(ax4, 'Color', 'k', 'XTick', [], 'YTick', []);
    axis(ax4, 'off');
    info_str = { ...
        sprintf('Subject: %s (%d/%d)', subj, idx, n), ...
        '', ...
        sprintf('Voxel:  [%d, %d, %d]', pos), ...
        sprintf('World:  [%.1f, %.1f, %.1f]', world(1:3)), ...
        '', ...
        'Click to move crosshair', ...
        'Enter = confirm origin', ...
        'n = skip,  q = quit'};
    text(ax4, 0.1, 0.5, info_str, 'Color', 'w', 'FontSize', 11, ...
         'VerticalAlignment', 'middle', 'Interpreter', 'none');

    drawnow;
end

function update_affine(fname, M_new)
% Update the sform (qform untouched) of a NIfTI file using SPM.
    V = spm_vol(fname);
    spm_get_space(fname, M_new);
end

function fname = find_raw_nii(subj_dir, name)
% Find raw mr.nii or ct.nii (not pp_).
    fname = fullfile(subj_dir, [name '.nii']);
    if exist(fname, 'file'), return; end
    gz = fullfile(subj_dir, [name '.nii.gz']);
    if exist(gz, 'file')
        gunzip(gz, subj_dir);
        return;
    end
    fname = '';
end

function v = clamp(v, lo, hi)
    v = max(lo, min(hi, v));
end
