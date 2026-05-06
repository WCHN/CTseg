function step3_segment_ct_ctseg(cfg)
% Step 3: CTseg on CT images using the SPM-space template.
% Produces native-space GM, WM, CSF tissue maps and modulated warped
% tissue maps (at 1.5mm in SPM/MNI space, matching SPM mwc output).
% Saves CTseg volume estimates (TBV/TIV) and runtime.
% Warps CT to MNI space via the CTseg deformation (which maps directly
% to MNI space when using an MNI-aligned atlas).

    fprintf('=== Step 3: CTseg on CT ===\n');
    subjects = get_subjects(cfg);
    failures = {};

    % Resolve atlas shorthand to file path
    pth_mu = resolve_atlas(cfg.mu);

    for i = 1:numel(subjects)
        subj_dir = fullfile(cfg.data_dir, subjects(i).name);

        % Find input file: prefer deformed CT, fall back to original
        nii_file = find_nii(subj_dir, 'pp_ct_def');
        if isempty(nii_file)
            nii_file = find_nii(subj_dir, 'pp_ct');
        end
        if isempty(nii_file), continue; end

        % Check for last file written (renamed mwc)
        done_file = fullfile(subj_dir, 'mwc1_ctseg_mni.nii');
        if exist(done_file, 'file')
            fprintf('  [%d/%d] %s - SKIP\n', i, numel(subjects), subjects(i).name);
            continue;
        end

        fprintf('  [%d/%d] %s\n', i, numel(subjects), subjects(i).name);

        try
            % --- CTseg segmentation ---
            % Only run if native tissue maps don't exist yet
            d = dir(fullfile(subj_dir, 'c01_*_CTseg.nii'));
            if isempty(d)
                tc = false(6, 3);
                tc(1:3, 1) = true;  % native space for GM, WM, CSF
                tc(1:3, 3) = true;  % modulated warped for GM, WM, CSF

                t_start = tic;
                % Use SPM-space template; crop output to SPM TPM bounding box
                bb = spm_get_bbox(fullfile(spm('Dir'), 'tpm', 'TPM.nii'), 'old');
                [res, vol] = spm_CTseg(nii_file, subj_dir, tc, true, false, false, NaN, cfg.v_settings, [], cfg.mu, false, bb);
                runtime = toc(t_start);

                save(fullfile(subj_dir, 'vol_CTseg.mat'), 'vol', 'runtime');
                fprintf('  CTseg runtime: %.1f s\n', runtime);

                % Clean up mb_fit_CTseg.mat
                mb_fit = fullfile(subj_dir, 'mb_fit_CTseg.mat');
                if exist(mb_fit, 'file'), delete(mb_fit); end
            else
                fprintf('  CTseg native maps exist, skipping segmentation\n');
            end

            % --- Rename mwc files to standard names ---
            for k = 1:3
                d_mwc = dir(fullfile(subj_dir, sprintf('mwc%02d_*_CTseg.nii', k)));
                out_name = fullfile(subj_dir, sprintf('mwc%d_ctseg_mni.nii', k));
                if ~isempty(d_mwc) && ~exist(out_name, 'file')
                    movefile(fullfile(subj_dir, d_mwc(1).name), out_name);
                    fprintf('  Renamed %s -> mwc%d_ctseg_mni.nii\n', d_mwc(1).name, k);
                end
            end

            % --- Warp CT to SPM space ---
            % CTseg deformation is native→template, so use push
            warp_file = fullfile(subj_dir, 'w_ctseg_pp_ct.nii');
            if ~exist(warp_file, 'file')
                d_y = dir(fullfile(subj_dir, 'y_*_CTseg.nii'));
                if ~isempty(d_y)
                    pth_y = fullfile(subj_dir, d_y(1).name);
                    fprintf('  Warping CT to SPM space...\n');
                    matlabbatch = {};
                    matlabbatch{1}.spm.util.defs.comp{1}.def                 = {pth_y};
                    matlabbatch{1}.spm.util.defs.out{1}.push.fnames          = {nii_file};
                    matlabbatch{1}.spm.util.defs.out{1}.push.weight          = {''};
                    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {subj_dir};
                    matlabbatch{1}.spm.util.defs.out{1}.push.fov.file        = {pth_mu};
                    matlabbatch{1}.spm.util.defs.out{1}.push.prefix          = 'w_ctseg_';
                    spm_jobman('run', matlabbatch);
                end
            end
        catch ME
            fprintf('  ** FAILED: %s — %s\n', subjects(i).name, ME.message);
            failures{end+1} = struct('subject', subjects(i).name, ...
                'step', 'step3_segment_ct_ctseg', 'error', ME.message); %#ok<AGROW>
        end

        cleanup_nii(nii_file);
    end

    save_failures(cfg, 'step3', failures);
end
