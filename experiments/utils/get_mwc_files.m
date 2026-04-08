function files = get_mwc_files(subj_dir, method)
% Return cell array of modulated warped tissue map paths {GM, WM, CSF}.
%   method: 'mr_spm', 'ct_spm', or 'ct_ctseg'
    switch method
        case 'mr_spm'
            files = {fullfile(subj_dir, 'mwc1pp_mr.nii'), ...
                     fullfile(subj_dir, 'mwc2pp_mr.nii'), ...
                     fullfile(subj_dir, 'mwc3pp_mr.nii')};
        case 'ct_spm'
            % Prefer deformed CT outputs, fall back to original
            if exist(fullfile(subj_dir, 'mwc1pp_ct_def.nii'), 'file')
                files = {fullfile(subj_dir, 'mwc1pp_ct_def.nii'), ...
                         fullfile(subj_dir, 'mwc2pp_ct_def.nii'), ...
                         fullfile(subj_dir, 'mwc3pp_ct_def.nii')};
            else
                files = {fullfile(subj_dir, 'mwc1pp_ct.nii'), ...
                         fullfile(subj_dir, 'mwc2pp_ct.nii'), ...
                         fullfile(subj_dir, 'mwc3pp_ct.nii')};
            end
        case 'ct_ctseg'
            files = {fullfile(subj_dir, 'mwc1_ctseg_mni.nii'), ...
                     fullfile(subj_dir, 'mwc2_ctseg_mni.nii'), ...
                     fullfile(subj_dir, 'mwc3_ctseg_mni.nii')};
        otherwise
            error('Unknown method: %s', method);
    end
end
