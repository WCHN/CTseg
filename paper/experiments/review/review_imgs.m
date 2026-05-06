function review_imgs(dir_data, file_prefix)
% Interactive viewer to review MR and CT images and flag subjects.
%
% Usage:
%   review_imgs('C:\Users\mbrudfors\Data\HN')              % raw data (mr.nii, ct.nii)
%   review_imgs('C:\Users\mbrudfors\Data\HN', 'pp_')       % preprocessed (pp_mr.nii, pp_ct.nii)
%
% Controls:
%   Right arrow / N  — next subject
%   Left arrow / P   — previous subject
%   F                — flag current subject (prints to command window)
%   Q                — quit and print summary
%
% Output: prints flagged subject names to command window for copy/paste.

if nargin < 1
    dir_data = 'C:\Users\mbrudfors\Data\HN';
end
if nargin < 2
    file_prefix = '';
end

% Find all subject directories
d = dir(dir_data);
subjects = d([d.isdir] & ~startsWith({d.name}, '.') & ~strcmp({d.name}, 'overviews') & ~strcmp({d.name}, 'results'));
n = numel(subjects);
fprintf('Found %d subjects in %s\n', n, dir_data);

flagged = false(n, 1);
idx = 1;

% Create figure
fig = figure('Name', 'Subject Review', 'Position', [50 50 1200 700], ...
             'Color', 'w', 'KeyPressFcn', @keypress);

show_subject();

    function keypress(~, event)
        switch event.Key
            case {'rightarrow', 'n'}
                idx = min(idx + 1, n);
                show_subject();
            case {'leftarrow', 'p'}
                idx = max(idx - 1, 1);
                show_subject();
            case 'f'
                flagged(idx) = ~flagged(idx);
                if flagged(idx)
                    fprintf('FLAGGED: %s\n', subjects(idx).name);
                else
                    fprintf('UNFLAGGED: %s\n', subjects(idx).name);
                end
                show_subject();
            case 'q'
                print_summary();
                close(fig);
        end
    end

    function show_subject()
        subj_dir = fullfile(dir_data, subjects(idx).name);

        % Find MR and CT images (prefer deformed CT when using pp_ prefix)
        mr_file = find_image(subj_dir, [file_prefix 'mr']);
        ct_file = '';
        if strcmp(file_prefix, 'pp_')
            ct_file = find_image(subj_dir, 'pp_ct_def');
        end
        if isempty(ct_file)
            ct_file = find_image(subj_dir, [file_prefix 'ct']);
        end

        clf(fig);

        % Row 1: MR (axial, coronal, sagittal)
        if ~isempty(mr_file)
            show_row(mr_file, 1, 'MR');
        else
            subplot(2, 3, 2);
            text(0.5, 0.5, 'MR not found', 'HorizontalAlignment', 'center', 'FontSize', 14);
            axis off;
        end

        % Row 2: CT (axial, coronal, sagittal)
        if ~isempty(ct_file)
            show_row(ct_file, 2, 'CT', [-50 100]);
        else
            subplot(2, 3, 5);
            text(0.5, 0.5, 'CT not found', 'HorizontalAlignment', 'center', 'FontSize', 14);
            axis off;
        end

        flag_str = '';
        if flagged(idx), flag_str = ' [FLAGGED]'; end
        sgtitle(sprintf('[%d/%d] %s%s  (F=flag, arrows=navigate, Q=quit)', ...
            idx, n, subjects(idx).name, flag_str), 'FontSize', 13);
    end

    function show_row(nii_file, row, label, clim)
        V = spm_vol(nii_file);
        Y = spm_read_vols(V);
        dims = size(Y);

        ax = round(0.45 * dims(3));
        co = round(0.55 * dims(2));
        sa = round(0.50 * dims(1));

        slices = {rot90(squeeze(Y(:,:,ax))), ...
                  rot90(squeeze(Y(:,co,:))), ...
                  rot90(squeeze(Y(sa,:,:)))};
        view_titles = {'Axial', 'Coronal', 'Sagittal'};

        for vi = 1:3
            subplot(2, 3, (row-1)*3 + vi);
            if nargin >= 4
                imagesc(slices{vi}, clim);
            else
                imagesc(slices{vi});
            end
            colormap(gray); axis image off;
            if row == 1, title(view_titles{vi}, 'FontSize', 11); end
            if vi == 1, ylabel(label, 'FontSize', 12); end
        end
    end

    function f = find_image(subj_dir, name)
        % Try .nii then .nii.gz
        f = fullfile(subj_dir, [name '.nii']);
        if exist(f, 'file'), return; end
        gz = fullfile(subj_dir, [name '.nii.gz']);
        if exist(gz, 'file')
            gunzip(gz, subj_dir);
            return;
        end
        f = '';
    end

    function print_summary()
        flagged_names = {subjects(flagged).name};
        fprintf('\n=== Flagged subjects (%d/%d) ===\n', sum(flagged), n);
        for i = 1:numel(flagged_names)
            fprintf('  %s\n', flagged_names{i});
        end
        if ~isempty(flagged_names)
            fprintf('\nMATLAB cell array:\n');
            fprintf('exclude = {');
            for i = 1:numel(flagged_names)
                if i > 1, fprintf(', '); end
                fprintf('''%s''', flagged_names{i});
            end
            fprintf('};\n');
        end
    end

end
