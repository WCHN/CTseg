function step2_segment_ct_spm(cfg)
% Step 2: SPM unified segmentation on CT images.
% Produces native-space GM, WM, CSF tissue maps, modulated warped tissue
% maps (mwc1-3pp_ct.nii), and normalised CT (wpp_ct.nii).
% Saves runtime (segmentation only, not normalise-write) to runtime_spm_ct.mat.

    fprintf('=== Step 2: SPM segmentation on CT ===\n');
    subjects = get_subjects(cfg);
    failures = {};
    pth_mu = resolve_atlas(cfg.mu);

    for i = 1:numel(subjects)
        subj_dir = fullfile(cfg.data_dir, subjects(i).name);

        % Find input file: prefer deformed CT, fall back to original
        nii_file = find_nii(subj_dir, 'pp_ct_def');
        if isempty(nii_file)
            nii_file = find_nii(subj_dir, 'pp_ct');
        end
        if isempty(nii_file), continue; end

        % Check for last file written (normalised CT)
        [~, nam] = fileparts(nii_file);
        done_file = fullfile(subj_dir, ['w' nam '.nii']);
        if exist(done_file, 'file')
            fprintf('  [%d/%d] %s - SKIP\n', i, numel(subjects), subjects(i).name);
            continue;
        end

        fprintf('  [%d/%d] %s\n', i, numel(subjects), subjects(i).name);

        try
            t_start = tic;
            run_spm_segment(nii_file, true, true, true);
            runtime = toc(t_start);

            save(fullfile(subj_dir, 'runtime_spm_ct.mat'), 'runtime');
            fprintf('  SPM-CT runtime: %.1f s\n', runtime);

            % Warp CT to MNI using SPM-CT deformation (pull, not timed)
            [~, nam] = fileparts(nii_file);
            def_file = fullfile(subj_dir, ['y_' nam '.nii']);
            warp_to_mni_spm(def_file, nii_file, pth_mu, subj_dir, 'w');

            % Keep deformation field (needed for step4b normalisation metrics)

            % Also warp CT to MNI using SPM-MR deformation (for step5 comparison)
            [~, nam_mr] = fileparts(find_nii(subj_dir, 'pp_mr'));
            y_mr_file = fullfile(subj_dir, ['y_' nam_mr '.nii']);
            if exist(y_mr_file, 'file')
                warp_to_mni_spm(y_mr_file, nii_file, pth_mu, subj_dir, 'wmr_');
            end
        catch ME
            fprintf('  ** FAILED: %s — %s\n', subjects(i).name, ME.message);
            failures{end+1} = struct('subject', subjects(i).name, ...
                'step', 'step2_segment_ct_spm', 'error', ME.message); %#ok<AGROW>
        end

        cleanup_nii(nii_file);
    end

    save_failures(cfg, 'step2', failures);
end


function warp_to_mni_spm(def_file, src_file, pth_mu_spm, subj_dir, prefix)
% Warp a native-space image to MNI using SPM deformation (pull).
    matlabbatch = {};
    matlabbatch{1}.spm.util.defs.comp{1}.comp{1}.def         = {def_file};
    matlabbatch{1}.spm.util.defs.comp{1}.space               = {pth_mu_spm};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {src_file};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {subj_dir};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 4;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = prefix;
    spm_jobman('run', matlabbatch);
end
