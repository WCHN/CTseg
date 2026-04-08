function cfg = config()
% Configuration for CTseg validation experiments.
% Edit this file to match your local setup.

    % Paths
    cfg.data_dir     = 'C:\Users\mbrudfors\Data\HN';
    cfg.exp_dir      = fileparts(mfilename('fullpath'));  % experiments/
    cfg.out_dir      = fullfile(cfg.exp_dir, 'results');  % experiments/results/
    cfg.fig_dir      = fullfile(cfg.exp_dir, '..', 'manuscript', 'figures');  % manuscript/figures/

    % External dependencies (edit if installed elsewhere)
    cfg.spm_preproc_dir    = 'C:\Users\mbrudfors\Code\spm-hospital-preproc';
    cfg.predict_pronto_dir = 'C:\Users\mbrudfors\Code\PredictPRoNTo';

    % CTseg template (SPM-space aligned; change to use a different template)
    cfg.pth_mu = fullfile(fileparts(which('spm_CTseg')), 'mu_CTseg_spm15.nii');

    % CTseg v_settings (spatial regularisation multiplier; empty = use default)
    cfg.v_settings = [];

    % Tissue classes
    cfg.tissue_names = {'GM', 'WM', 'CSF'};
    cfg.n_tissues    = numel(cfg.tissue_names);

    % Analysis settings
    cfg.dice_thresh = 0.5;   % binarisation threshold for tissue probability maps
    cfg.alpha       = 0.05;  % significance level

    % Excluded subjects (incomplete brain coverage, reviewed via review_subjects.m)
    cfg.exclude = { ...
        '1HNA098','1HNA099','1HNA100','1HNA102','1HNA103','1HNA104', ...
        '1HNA105','1HNA106','1HNA107','1HNA108','1HNA109','1HNA110', ...
        '1HNA113','1HNA115','1HNA116','1HNA117','1HNA119','1HNA120', ...
        '1HNA121','1HNA124','1HNA126','1HNA129','1HNA130','1HNA132', ...
        '1HNA133','1HNA135','1HNA136','1HNA138','1HNA139','1HNA141', ...
        '1HNA142','1HNA143', ...
        '1HNC001','1HNC002','1HNC003','1HNC004','1HNC005','1HNC007', ...
        '1HNC008','1HNC012','1HNC014','1HNC017','1HNC019','1HNC020', ...
        '1HNC021','1HNC022','1HNC023','1HNC025','1HNC029','1HNC031', ...
        '1HNC035','1HNC036','1HNC037','1HNC038','1HNC039','1HNC040', ...
        '1HNC043','1HNC045','1HNC046','1HNC050','1HNC061','1HNC066', ...
        '1HNC067','1HNC068','1HNC071','1HNC072','1HNC073','1HNC076', ...
        '1HNC082','1HNC083','1HNC084','1HNC085','1HNC087','1HNC088', ...
        '1HNC094','1HNC098','1HNC099','1HNC101','1HNC102','1HNC103', ...
        '1HNC104','1HNC105','1HNC107','1HNC109','1HNC110','1HNC111', ...
        '1HNC112','1HNC114','1HNC117','1HNC118','1HNC120','1HNC121', ...
        '1HNC124','1HNC125','1HNC127','1HNC128','1HNC130'};

    % Test mode: run on a few subjects for quick debugging
    cfg.test_mode = false;       % set true to limit number of subjects
    cfg.test_n_subjects = 1;     % number of subjects in test mode
    cfg.test_subjects = {};      % specific subjects, e.g. {'1HNA001','1HNA004'}; empty = use first N

    % Create output directory
    if ~exist(cfg.out_dir, 'dir'), mkdir(cfg.out_dir); end
    if ~exist(cfg.fig_dir, 'dir'), mkdir(cfg.fig_dir); end
end
