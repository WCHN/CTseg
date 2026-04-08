% preproc.m — Preprocess paired MR/CT data for CTseg validation.
%
% Dataset: SynthRAD2025 Challenge, Task 1 (Head & Neck) training data.
%          https://synthrad2025.grand-challenge.org/
%
% Input:   Raw mr.nii.gz and ct.nii.gz per subject folder.
% Output:  pp_mr.nii.gz and pp_ct.nii.gz (preprocessed) in the same folder.
%
% Preprocessing steps (via spm-hospital-preproc):
%   1. Rigid alignment to MNI space
%   2. Co-registration of MR and CT using NMI
%   3. Resampling to SPM12 template bounding box
%
% Usage:
%   preproc              % process all subjects
%   preproc('test')      % process test subjects only (from config)
%   preproc('1HNA001')   % process single subject
%
% Requires:
%   - SPM12 on MATLAB path
%   - spm-hospital-preproc on MATLAB path
%     https://github.com/WTCN-computational-anatomy-group/spm-hospital-preproc

% 

function preproc(subject)

exp_dir = fileparts(fileparts(mfilename('fullpath')));  % experiments/
addpath(exp_dir);
addpath(fullfile(exp_dir, 'utils'));
cfg = config();
addpath(cfg.spm_preproc_dir);

dir_data = cfg.data_dir;
ext      = '.nii';
prefix   = '';

% Preprocessing options
opt             = struct;
opt.do.real_mni = true;           % Rigid alignment to MNI
opt.realign2mni.ix_realign  = 2;
opt.do.coreg = true;
opt.do.bb_spm   = true;           % SPM12 template bounding box

if nargin >= 1 && ~isempty(subject) && strcmpi(subject, 'test')
    % Test mode: use same test subjects as run_all
    cfg.test_mode = true;
    subjects = get_subjects(cfg);
    fprintf('Processing %d test subjects\n', numel(subjects));

    for i = 1:numel(subjects)
        subj_dir_i = fullfile(dir_data, subjects(i).name);
        fprintf('Processing %d/%d: %s\n', i, numel(subjects), subjects(i).name);
        files = spm_select('FPList', subj_dir_i, ['^' prefix '(mr|ct)\' ext '$']);
        files = cellstr(files);
        opt.dir_out = subj_dir_i;
        RunPreproc(files, opt);
    end
elseif nargin >= 1 && ~isempty(subject)
    % Single subject mode
    subj_dir = fullfile(dir_data, subject);
    if ~exist(subj_dir, 'dir')
        error('Subject folder not found: %s', subj_dir);
    end
    fprintf('Processing single subject: %s\n', subject);
    files = spm_select('FPList', subj_dir, ['^' prefix '(mr|ct)\' ext '$']);
    files = cellstr(files);
    opt.dir_out = 'output-test';
    RunPreproc(files, opt);
else
    % All subjects (excluding flagged via config)
    subjects = get_subjects(cfg);
    fprintf('Processing %d subjects\n', numel(subjects));

    parfor i = 1:numel(subjects)
        opt_cpy = opt;
        subj_dir_i = fullfile(dir_data, subjects(i).name);
        fprintf('Processing %d/%d: %s\n', i, numel(subjects), subjects(i).name);
        files = spm_select('FPList', subj_dir_i, ['^' prefix '(mr|ct)\' ext '$']);
        files = cellstr(files);
        opt_cpy.dir_out = subj_dir_i;
        RunPreproc(files, opt_cpy);
    end
end

end
