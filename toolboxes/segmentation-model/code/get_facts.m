function fct = get_facts(dat,model)
% FORMAT fct = get_facts(dat,model)
% dat   - Subject's data structure (one subject)
% model - Model structure
% fct   - Facts structure with fields:
%         * templ.dm  - Template dimensions 
%         * templ.mat - Template orientation matrix
%         * templ.K   - Number of template classes
%         * subj.dm   - Observed image dimensions
%         * subj.mat  - Observed image orientation matrix
%         * subj.vs   - Observed image voxel size
%         * subj.C    - Number of channels
%         * subj.nam  - Names of the different channels
%         * subj.I    - Number of voxels in observed image
%         * subj.mod  - Name of modality
%         * subj.ff   - Fudge factor
%
% Store a set of useful variables in a structure.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Template
fct.templ.dm  = model.template.nii.dat.dim;   % Template orientation matrix   
fct.templ.mat = model.template.nii.mat;       % Template orientation matrix                  
fct.templ.K   = fct.templ.dm(4);

% Subject
[dm,mat,vs,C,nam,V,fnames,chn_names] = obs_info(dat);
modality                             = dat.modality{1}.name; 

chn_names{end + 1}        = 'Template';
chn_names{end + 1}        = 'Z';

fct.subj.dm  = dm;
fct.subj.mat = mat;
fct.subj.vs  = vs;
fct.subj.C   = C;
fct.subj.nam = chn_names;
fct.subj.I   = prod(dm(1:3));
fct.subj.mod = modality;
fct.subj.ff  = get_ff(vs); % Fudge factor (see spm_preproc8)   
%==========================================================================