function nii_file = find_nii(subj_dir, name)
% Find a NIfTI file, handling both .nii and .nii.gz.
% Returns path to .nii file (gunzips if needed), or '' if not found.
%
%   nii_file = find_nii(subj_dir, 'pp_mr')
%   % Returns pp_mr.nii if it exists, or gunzips pp_mr.nii.gz

    nii_file = fullfile(subj_dir, [name '.nii']);
    if exist(nii_file, 'file'), return; end

    gz_file = fullfile(subj_dir, [name '.nii.gz']);
    if exist(gz_file, 'file')
        gunzip(gz_file, subj_dir);
        return;
    end

    nii_file = '';  % not found
end
