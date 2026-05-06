function cleanup_nii(nii_file)
% Delete unzipped .nii if the .nii.gz source exists (to save disk space).
    gz_file = [nii_file '.gz'];
    if exist(gz_file, 'file') && exist(nii_file, 'file')
        delete(nii_file);
    end
end
