function run_spm_segment(nii_file, write_def, is_ct, write_mwc)
% Run SPM unified segmentation on a single NIfTI file.
%   nii_file:  path to .nii file
%   write_def: if true, write forward deformation field (y_*.nii)
%   is_ct:     if true, suppress bias field correction (not needed for CT)
%   write_mwc: if true, write modulated warped tissue maps (mwc1-3*.nii)
    if nargin < 2, write_def = false; end
    if nargin < 3, is_ct = false; end
    if nargin < 4, write_mwc = false; end

    tpm_file = fullfile(spm('Dir'), 'tpm', 'TPM.nii');

    % Bias field settings: heavy regularisation for CT (effectively off)
    if is_ct
        biasreg = 1;
    else
        biasreg = 0.001*(1/5);
    end
    biasfwhm = 60;

    matlabbatch = {};
    matlabbatch{1}.spm.spatial.preproc.channel.vols     = {[nii_file ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = biasreg;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = biasfwhm;
    matlabbatch{1}.spm.spatial.preproc.channel.write     = [0 0];

    ngaus_list = [1 1 2 3 4 2];
    for k = 1:6
        matlabbatch{1}.spm.spatial.preproc.tissue(k).tpm   = {[tpm_file ',' num2str(k)]};
        matlabbatch{1}.spm.spatial.preproc.tissue(k).ngaus  = ngaus_list(k);
        if k <= 3
            matlabbatch{1}.spm.spatial.preproc.tissue(k).native = [1 0];
            % warped = [normalised, modulated_normalised]
            matlabbatch{1}.spm.spatial.preproc.tissue(k).warped = [0 write_mwc];
        else
            matlabbatch{1}.spm.spatial.preproc.tissue(k).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(k).warped = [0 0];
        end
    end

    matlabbatch{1}.spm.spatial.preproc.warp.mrf      = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup  = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg      = [0 0.001 0.5 0.05 0.2]*0.1;
    matlabbatch{1}.spm.spatial.preproc.warp.affreg   = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm     = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp     = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write    = [0 write_def];

    spm_jobman('run', matlabbatch);
end
