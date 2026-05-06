function files = get_tissue_files(subj_dir, method)
% Return cell array of tissue map paths {GM, WM, CSF} for a given method.
%   method: 'mr_spm', 'ct_spm', or 'ct_ctseg'
    switch method
        case 'mr_spm'
            files = {fullfile(subj_dir, 'c1pp_mr.nii'), ...
                     fullfile(subj_dir, 'c2pp_mr.nii'), ...
                     fullfile(subj_dir, 'c3pp_mr.nii')};
        case 'ct_spm'
            % Prefer deformed CT outputs, fall back to original
            if exist(fullfile(subj_dir, 'c1pp_ct_def.nii'), 'file')
                files = {fullfile(subj_dir, 'c1pp_ct_def.nii'), ...
                         fullfile(subj_dir, 'c2pp_ct_def.nii'), ...
                         fullfile(subj_dir, 'c3pp_ct_def.nii')};
            else
                files = {fullfile(subj_dir, 'c1pp_ct.nii'), ...
                         fullfile(subj_dir, 'c2pp_ct.nii'), ...
                         fullfile(subj_dir, 'c3pp_ct.nii')};
            end
        case 'ct_ctseg'
            % CTseg names include MB subject index and temp prefix,
            % e.g. c01_1_00001_temp_pp_ct_CTseg.nii — find by pattern
            files = cell(1, 3);
            for k = 1:3
                pattern = fullfile(subj_dir, sprintf('c%02d_*_CTseg.nii', k));
                d = dir(pattern);
                if ~isempty(d)
                    files{k} = fullfile(subj_dir, d(1).name);
                else
                    files{k} = fullfile(subj_dir, sprintf('c%02d_CTseg.nii', k));
                end
            end
        otherwise
            error('Unknown method: %s', method);
    end
end
