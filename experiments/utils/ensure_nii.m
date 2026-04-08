function nii_file = ensure_nii(gz_file, subj_dir)
% Gunzip a .nii.gz file if the .nii doesn't already exist.
% Returns path to the .nii file.
    [~, nam] = fileparts(gz_file);  % 'pp_mr.nii'
    nii_file = fullfile(subj_dir, nam);
    if ~exist(nii_file, 'file')
        gunzip(gz_file, subj_dir);
    end
end
