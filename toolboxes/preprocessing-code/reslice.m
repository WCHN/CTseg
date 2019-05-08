function Nii = reslice(Nii,ref_ix,deg)
if nargin < 2, ref_ix = 1; end
if nargin < 3, deg    = 4; end

fprintf('Reslicing...')
N = numel(Nii{1});
if N == 1
    return
end

ixs       = 1:N;
source_ix = ixs(ixs ~= ref_ix);

% Use SPM batch job to reslice
fnames = {Nii{1}(ref_ix).dat.fname, Nii{1}(source_ix).dat.fname}';
run_reslice(fnames,deg);

% Update NIIs
for n=source_ix
    f             = Nii{1}(n).dat.fname;
    [pth,nam,ext] = fileparts(f);
    nf            = fullfile(pth,['r' nam ext]);
    
    delete(f);
    Nii{1}(n) = nifti(nf);
end

if numel(Nii) > 1
    % Reslice labels too
    V = spm_vol(Nii{1}(ref_ix).dat.fname);
    for n=1:N
        if isempty(Nii{2}(n).dat), continue; end
        
        Nii{2}(n) = reslice_labels(V,Nii{2}(n));
    end    
end
fprintf('done!\n')
%==========================================================================

%==========================================================================
function V = run_reslice(fnames,deg)
matlabbatch{1}.spm.spatial.realign.write.data            = fnames;
matlabbatch{1}.spm.spatial.realign.write.roptions.which  = [1 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.interp = deg;
matlabbatch{1}.spm.spatial.realign.write.roptions.wrap   = [0 0 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.mask   = 1;
matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);

V = spm_vol(fnames{1});
for n=2:numel(fnames)
    f             = fnames{n};
    [pth,nam,ext] = fileparts(f);
    nf            = fullfile(pth,['r' nam ext]);
    V(n)          = spm_vol(nf);
end
%==========================================================================

%==========================================================================
function Niio = reslice_labels(V_ref,Nii)
% Parameters
deg = 1;
dt  = [spm_type('float32') spm_platform('bigend')];

% Get labels, etc
labels        = Nii.dat(:,:,:);
lkp           = unique(labels)';
K             = numel(lkp); % Number of labels
fname         = Nii.dat.fname;
[pth,nam,ext] = fileparts(fname);

% Iterate over each label and create a resliced label image (nlabels)
dm      = V_ref.dim;
labelso = zeros([dm K],'single');
cnt     = 1;
for k=lkp
    labels_k = single(labels == k);

    % Smooth
    labels_k = convn(labels_k,reshape([0.25 0.5 0.25],[3,1,1]),'same');
    labels_k = convn(labels_k,reshape([0.25 0.5 0.25],[1,3,1]),'same');
    labels_k = convn(labels_k,reshape([0.25 0.5 0.25],[1,1,3]),'same');
                
    fname_k = fullfile(pth,['n' nam ext]);    
    create_nii(fname_k,labels_k,Nii.mat,dt,'Resliced labels');        
        
    % Reslice
    fnames   = {V_ref.fname, fname_k}';
    Vo       = run_reslice(fnames,deg);
    labels_k = single(Vo(2).private.dat(:,:,:));

    labelso(:,:,:,cnt) = cat(4,labels_k);
    
    delete(fname_k);
    
    cnt = cnt + 1;
end
clear labels labels_k
delete(fname);

% Get MLs of resliced labels
[~,ml] = max(labelso,[],4);
ml     = ml - 1;

% Write output
Niio            = nifti(Vo(2).fname);
Niio.dat(:,:,:) = ml;
%==========================================================================

%==========================================================================
function Nii = create_nii(pth,dat,mat,dtype,descrip,offset,scl_slope,scl_inter)
if nargin<6, offset    = 0; end
if nargin<7, scl_slope = 1; end
if nargin<8, scl_inter = 0; end

if exist(pth,'file')==2, delete(pth); end

Nii         = nifti;
dm          = size(dat);
Nii.dat     = file_array(pth,dm,dtype,offset,scl_slope,scl_inter);
Nii.mat     = mat;
Nii.mat0    = mat;
Nii.descrip = descrip;
create(Nii);

if numel(dm)==4
    Nii.dat(:,:,:,:) = dat;
else
    Nii.dat(:,:,:)   = dat;
end
%==========================================================================