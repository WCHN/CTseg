function step1_segment_mr(cfg)
% Step 1: SPM unified segmentation on MR images.
% Produces native-space GM, WM, CSF tissue maps (c1/c2/c3pp_mr.nii)
% and modulated warped tissue maps (mwc1/mwc2/mwc3pp_mr.nii).
% Saves runtime to runtime_spm_mr.mat per subject.

    fprintf('=== Step 1: SPM segmentation on MR ===\n');
    subjects = get_subjects(cfg);
    failures = {};

    for i = 1:numel(subjects)
        subj_dir = fullfile(cfg.data_dir, subjects(i).name);

        % Find input file (.nii or .nii.gz)
        nii_file = find_nii(subj_dir, 'pp_mr');
        if isempty(nii_file), continue; end

        % Check for last file written
        done_file = fullfile(subj_dir, 'mwc1pp_mr.nii');
        if exist(done_file, 'file')
            fprintf('  [%d/%d] %s - SKIP\n', i, numel(subjects), subjects(i).name);
            continue;
        end

        fprintf('  [%d/%d] %s\n', i, numel(subjects), subjects(i).name);

        try
            t_start = tic;
            run_spm_segment(nii_file, true, false, true);
            runtime = toc(t_start);

            save(fullfile(subj_dir, 'runtime_spm_mr.mat'), 'runtime');
            fprintf('  SPM-MR runtime: %.1f s\n', runtime);

        catch ME
            fprintf('  ** FAILED: %s — %s\n', subjects(i).name, ME.message);
            failures{end+1} = struct('subject', subjects(i).name, ...
                'step', 'step1_segment_mr', 'error', ME.message); %#ok<AGROW>
        end

        cleanup_nii(nii_file);
    end

    save_failures(cfg, 'step1', failures);
end
