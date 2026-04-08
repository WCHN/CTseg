% run_all.m — Reproduce all results for CTseg validation paper.
%
% Prerequisites:
%   - SPM12 on MATLAB path (with Shoot and Longitudinal toolboxes)
%   - CTseg on MATLAB path (https://github.com/WCHN/CTseg)
%   - Multi-Brain toolbox (https://github.com/WTCN-computational-anatomy-group/mb)
%   - spm-hospital-preproc (https://github.com/WTCN-computational-anatomy-group/spm-hospital-preproc)
%   - Image Processing Toolbox (for surface distance metrics)
%   - Statistics Toolbox (for signrank, prctile)
%   - Preprocessed data in the path specified in config.m
%
% Usage:
%   >> run_all              % runs everything
%   >> run_all('step4')     % re-run a single step
%   >> run_all('test')      % test mode: full pipeline on limited subjects
%   >> run_all('test4')     % test mode: single step on limited subjects
%
% Steps 1-3 are slow (segmentation, run once per dataset).
% Steps 4-6 compute metrics (moderate, re-run when adding metrics).
% Analysis scripts are fast (iterate freely for figure/table tweaks).

function run_all(mode)
    if nargin < 1, mode = ''; end

    % Setup paths
    exp_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(exp_dir, 'utils'));
    addpath(fullfile(exp_dir, 'steps'));
    addpath(fullfile(exp_dir, 'analysis'));
    addpath(fullfile(exp_dir, 'atlas'));
    addpath(fullfile(exp_dir, 'preprocessing'));
    addpath(fullfile(exp_dir, 'review'));
    addpath(fullfile(exp_dir, '..'));  % CTseg root

    cfg = config();

    % Add spm-hospital-preproc
    if exist(cfg.spm_preproc_dir, 'dir')
        addpath(cfg.spm_preproc_dir);
    end

    % Add PredictPRoNTo and PRoNTo
    if exist(cfg.predict_pronto_dir, 'dir')
        addpath(cfg.predict_pronto_dir);
        % PRoNTo v2 may be bundled as a zip — check for extracted folder
        pronto_dir = fullfile(cfg.predict_pronto_dir, 'PRoNTo_dev-2.1.2');
        if exist(pronto_dir, 'dir')
            addpath(genpath(pronto_dir));
        end
    end

    % Verify dependencies
    check_dependencies(cfg);

    spm('defaults', 'fmri');
    spm_jobman('initcfg');

    % Handle clean modes
    if strcmpi(mode, 'clean')
        clean_outputs(cfg, false);
        return;
    elseif strcmpi(mode, 'testclean')
        cfg.test_mode = true;
        clean_outputs(cfg, true);
        return;
    end

    % Parse mode: 'test', 'test4', 'step4', etc.
    only_step = '';
    if startsWith(mode, 'test', 'IgnoreCase', true)
        cfg.test_mode = true;
        only_step = regexprep(mode, '^test', '', 'ignorecase');
        % 'test' -> '' (run all), 'test4' -> '4' -> 'step4'
        % 'testfigures' -> 'figures' (keep as-is for non-numeric modes)
        if ~isempty(only_step) && ~isnan(str2double(only_step))
            only_step = ['step' only_step];  % numeric: '4' -> 'step4'
        end
        fprintf('*** TEST MODE: processing %d subject(s) ***\n\n', cfg.test_n_subjects);
    else
        only_step = mode;
    end

    % --- Processing (slow, run once) ---
    if should_run('step1', only_step), step1_segment_mr(cfg);         end
    if should_run('step2', only_step), step2_segment_ct_spm(cfg);     end
    if should_run('step3', only_step), step3_segment_ct_ctseg(cfg);   end

    % --- Metrics (moderate) ---
    if should_run('step4', only_step),  step4_compute_metrics(cfg);         end
    if should_run('step4b', only_step), step4b_normalisation_metrics(cfg); end
    if should_run('step5', only_step), step5_average_normalised(cfg); end
    if should_run('step6', only_step), step6_volumetrics(cfg);        end
    if should_run('step7', only_step), step7_prediction(cfg);         end

    % --- Analysis (fast, iterate freely) ---
    if should_run('stats',   only_step), run_stats(cfg);    end
    if should_run('figures', only_step), make_figures(cfg);  end
    if should_run('tables',  only_step), make_tables(cfg);  end

    fprintf('\n=== All done. Results in: %s ===\n', cfg.out_dir);
end

function tf = should_run(name, only)
    tf = isempty(only) || strcmpi(name, only);
end

function clean_outputs(cfg, test_only)
% Delete all generated outputs so steps can be re-run from scratch.
%   test_only: if true, only clean test subjects; if false, clean all.
    if test_only
        subjects = get_subjects(cfg);
        fprintf('Cleaning %d test subject(s)...\n', numel(subjects));
    else
        subjects = get_subjects(cfg);
        fprintf('Cleaning ALL %d subjects...\n', numel(subjects));
    end

    % Patterns to delete per subject (all generated files)
    % MATLAB dir() doesn't support [123] brackets, so list individually
    patterns = { ...
        'c1pp_mr.nii', 'c2pp_mr.nii', 'c3pp_mr.nii', ...
        'c1pp_ct.nii', 'c2pp_ct.nii', 'c3pp_ct.nii', ...
        'c1pp_ct_def.nii', 'c2pp_ct_def.nii', 'c3pp_ct_def.nii', ...
        'mwc1pp_mr.nii', 'mwc2pp_mr.nii', 'mwc3pp_mr.nii', ...
        'mwc1pp_ct.nii', 'mwc2pp_ct.nii', 'mwc3pp_ct.nii', ...
        'mwc1pp_ct_def.nii', 'mwc2pp_ct_def.nii', 'mwc3pp_ct_def.nii', ...
        'mwc1_ctseg_mni.nii', 'mwc2_ctseg_mni.nii', 'mwc3_ctseg_mni.nii', ...
        'c01_*_CTseg.nii', 'c02_*_CTseg.nii', 'c03_*_CTseg.nii', ...
        'mwc01_*_CTseg.nii', 'mwc02_*_CTseg.nii', 'mwc03_*_CTseg.nii', ...
        'wpp_ct.nii', 'wpp_ct_def.nii', 'wmr_pp_ct.nii', 'wmr_pp_ct_def.nii', ...
        'w_ctseg_pp_ct.nii', 'w_ctseg_pp_ct_def.nii', ...
        'wnorm_spm_*.nii', 'wnorm_spmct_*.nii', 'wnorm_mr_*.nii', 'wnorm_ctseg_*.nii', 'wmu_*.nii', ...
        'y_*.nii', 'v_*.nii', ...
        '*_seg8.mat', 'mb_fit_CTseg.mat', ...
        'vol_CTseg.mat', 'runtime_*.mat', ...
        'temp_*.nii', 'tmp_ref_1mm.nii'};

    n_deleted = 0;
    for i = 1:numel(subjects)
        subj_dir = fullfile(cfg.data_dir, subjects(i).name);
        for p = 1:numel(patterns)
            files = dir(fullfile(subj_dir, patterns{p}));
            for f = 1:numel(files)
                delete(fullfile(subj_dir, files(f).name));
                n_deleted = n_deleted + 1;
            end
        end
    end

    % Clean results directory
    results_patterns = {'*.mat', '*.nii', '*.png', '*.pdf', '*.tex'};
    for p = 1:numel(results_patterns)
        files = dir(fullfile(cfg.out_dir, results_patterns{p}));
        for f = 1:numel(files)
            delete(fullfile(cfg.out_dir, files(f).name));
            n_deleted = n_deleted + 1;
        end
    end

    % Clean PRoNTo output directories and mwc link directories
    pronto_dirs = dir(fullfile(cfg.out_dir, 'pronto_*'));
    link_dirs   = dir(fullfile(cfg.out_dir, 'mwc_links_*'));
    for d = [pronto_dirs; link_dirs]'
        if d.isdir
            rmdir(fullfile(cfg.out_dir, d.name), 's');
            n_deleted = n_deleted + 1;
        end
    end

    fprintf('Deleted %d files/directories.\n', n_deleted);
end

function check_dependencies(cfg)
% Verify all required tools are on the MATLAB path.
    errors = {};

    if isempty(which('spm'))
        errors{end+1} = 'SPM12 not found. Download from https://www.fil.ion.ucl.ac.uk/spm/software/download/';
    end
    if isempty(which('spm_shoot3d'))
        errors{end+1} = 'SPM Shoot toolbox not found. Add spm12/toolbox/Shoot to path.';
    end
    if isempty(which('spm_CTseg'))
        errors{end+1} = 'CTseg not found. Download from https://github.com/WCHN/CTseg';
    end
    if isempty(which('spm_mb_fit'))
        errors{end+1} = 'Multi-Brain toolbox not found. Download from https://github.com/WTCN-computational-anatomy-group/mb and place in spm12/toolbox/mb.';
    end
    if isempty(which('RunPreproc'))
        errors{end+1} = sprintf('spm-hospital-preproc not found. Download from https://github.com/WTCN-computational-anatomy-group/spm-hospital-preproc and set cfg.spm_preproc_dir in config.m (currently: %s).', cfg.spm_preproc_dir);
    end
    if isempty(which('PredictPRoNTo'))
        errors{end+1} = sprintf('PredictPRoNTo not found. Clone from https://github.com/brudfors/PredictPRoNTo and set cfg.predict_pronto_dir in config.m (currently: %s).', cfg.predict_pronto_dir);
    end
    if isempty(which('pronto'))
        errors{end+1} = 'PRoNTo v2 not found. Extract PRoNTo_dev-2.1.2.zip from the PredictPRoNTo repo.';
    end
    if ~exist('bwperim', 'file')
        errors{end+1} = 'Image Processing Toolbox not found (needed for bwperim/bwdist in surface distance metrics).';
    end
    if ~exist('signrank', 'file')
        errors{end+1} = 'Statistics Toolbox not found (needed for signrank).';
    end
    if ~exist(cfg.data_dir, 'dir')
        errors{end+1} = sprintf('Data directory not found: %s. Edit config.m.', cfg.data_dir);
    end

    if ~isempty(errors)
        fprintf('\n*** MISSING DEPENDENCIES ***\n');
        for i = 1:numel(errors)
            fprintf('  [%d] %s\n', i, errors{i});
        end
        fprintf('\n');
        error('Please resolve the above dependencies before running.');
    end
end
