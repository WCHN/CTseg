function Nii = coreg(Nii,ref_ix)
fprintf('Co-register...')
N = numel(Nii{1});
if N==1
    return;
end
if nargin < 2, ref_ix = 1; end

% Set options
matlabbatch{1}.spm.spatial.coreg.estimate.ref               = {Nii{1}(ref_ix).dat.fname};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];

% Co-register
ixs       = 1:N;
source_ix = ixs(ixs~=ref_ix);
R         = cell(1,N);
for n=source_ix
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {Nii{1}(n).dat.fname};         

    res  = spm_jobman('run',matlabbatch);
    R{n} = res{1}.M;
end

% Update NIIs
for n=1:N
    f         = Nii{1}(n).dat.fname;
    Nii{1}(n) = nifti(f);
end

if numel(Nii) > 1
    % Keep labels in alignment
    for n=source_ix
        if isempty(Nii{2}(n).dat), continue; end
        
        f         = Nii{2}(n).dat.fname;
        mat0      = Nii{2}(n).mat;
        spm_get_space(f,R{n}\mat0); 
        Nii{2}(n) = nifti(f); 
    end    
end
fprintf('done!\n')
%=================================================================