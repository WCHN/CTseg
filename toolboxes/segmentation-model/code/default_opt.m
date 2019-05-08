function [opt,holly] = default_opt(opt)
% FORMAT [opt,holly] = default_opt(opt)
% opt   - Options structure
% holly - Parallel processing structure
%
% Defines default options.
% Options that were already set in the input structure are not overriden.
%
%--------------------------------------------------------------------------
%
% GENERAL
% -------
% opt.sched
%
% MODEL
% -----
% opt.model.it           - Number of current VEM iterations
% opt.model.tol          - Lower bound tolerance criterion
% opt.model.niter        - Maximum number of VEM iterations
% opt.model.nam_cls      - Names of Template classes
% opt.model.clean_up     - Delete temporary files at end of write_results()
% opt.model.PropPrior.do - Optimise prior on tissue proportions
%
% GMM
% -----
% opt.gmm.niter           - Number of iterations in each GMM sub-loop
% opt.gmm.tol             - Lower bound tolerance criterion
% opt.gmm.pth_GaussPrior  - Path to GMM parameters file on disk
% opt.gmm.pth_PropPrior   - Path to proportion parameters file on disk
% opt.gmm.hist.niter_main - Number of iteration for initial histogram fitting
% opt.gmm.hist.niter_gmm  - Number of sub-iterations
% opt.gmm.hist.init_ix    - Map cluster indices to Template classes
% opt.gmm.hist.verbose    - Verbosity level
% opt.gmm.hist.verbose_gmm- Verbosity level in sub-loop 
% opt.gmm.labels.cm       - Map manual labels to Template classes
% opt.gmm.labels.use      - Use manual labels
% opt.gmm.labels.S        - Rater sensitivity for labelled voxels
% opt.gmm.labels.Su       - Rater sensitivity for unlabelled voxels
% opt.gmm.GaussPrior.constrained - Constrain class variance to be similar
% opt.gmm.GaussPrior.verbose     - Verbosity when updating GMM prior
% opt.gmm.GaussPrior.mods  - If we wanna share, let's say, one prior over
%                            all CTs
%
% TEMPLATE
% -----
% opt.template.do           - Update template
% opt.template.pth_template - Path to template file on disk
% opt.template.K            - Number of template classes
% opt.template.vs           - Template voxel size
% opt.template.reg0         - Base regularisation parameter
% opt.template.reg          - Decreasing regularisation parameters
% opt.template.shrink       - Crop template to bounding box of data
% opt.template.load_a_der   - Read and write derivatives on disk
% opt.template.R            - Rotation matrix towards null space
% opt.template.sym          - Force symmetric template
% opt.template.niter        - Number of Gauss-Newton iterations, will only
%                             iterate if opt.template.load_a_der = false
% opt.template.verbose      - Verbosity level
% opt.template.bg_class     - Label of air background class, used as outside FOV 
%                             prior when warping template (see init_template.bg)
% opt.template.resize       - Crops the template to FOV of the default SPM
%                             template using resize_template()
% opt.template.keep_neck    - Keep or remove neck region in
%                             resize_template()
% opt.template.clean.its    - Iterations at when to clean template
% opt.template.clean.brain  - Brain classes for template cleaning
% opt.template.clean.les    - Lesion class for template cleaning
% opt.template.clean.air    - Air class for template cleaning
% opt.template.tc_miss      - Tell the algorithm to which tissue class the
%                             missing data belong
%
% REGISTRATION
% ------------
% opt.reg.rparam0       - Base regularisation parameter
% opt.reg.rparam        - Decresing regularisation parameters
% opt.reg.int_args      - Number of integration steps
% opt.reg.niter         - Number of Gauss-Newton iteratons
% opt.reg.tol           - Lower bound tolerance criterion
% opt.reg.strt_nl       - Iteration at which to start optimising non-linear warps
% opt.reg.mc_aff        - Mean correct affine parameters
% opt.reg.aff_type      - Type of affine transformation ['translation','rotation','rigid','similitude','affine']
% opt.reg.aff_reg       - Regularisation of affine transformation
% opt.reg.do_aff        - Optimise affine transform?
% opt.reg.do_nl         - Optimise non-linear warp?
% opt.reg.nit_init_aff  - Number of iteration for initial affine
%                         registration of template to subject
% opt.reg.init_aff_tol  - Tolerance for initial affine registration of template to subject
%
% SEGMENTATION
% ------------
% opt.seg.niter         - Maximum number of subject segmentation sub-iterations
% opt.seg.tol           - Lower bound tolerance criterion
% opt.seg.show          - Plot lower bound and images
% opt.seg.samp          - Sample a subset of voxels (makes everything faster)
% opt.seg.bg            - Outside of brain classes, e.g., air, soft tissue and
%                         skull. Used when doing ad-hoc clean-up in
%                         clean_brain().
% opt.seg.mrf.ml        - ML or VB MRF update
% opt.seg.mrf.val_diag  - Value on the diagonal of the MRF confusion matrix
% opt.seg.mrf.alpha     - Parameter of VB MRF
% opt.seg.mrf.niter     - Number of MRF update iterations
% opt.seg.mskonlynan    - Only mask NaN values in observed images
%
% BIAS FIELD
% ----------
% opt.bf.biasfwhm       - Full-width half max of the smallest basis function
% opt.bf.niter          - Maximum number of Gauss-Newton iterations
% opt.bf.tol            - Lower bound tolerance criterion
% opt.bf.mc_bf          - Do population-wise mean correction of bias field
% opt.bf.biasreg        - Regularisation parameter
% opt.bf.do             - Optimise bias field?
% opt.bf.mc_bf_verbose  - Verbose when population-wise mean correcting of bias field
%
% TISSUE PROPORTIONS
% ------------------
% opt.prop.niter        - Maximum number of subject proportion update sub-iterations
% opt.prop.gnniter      - Maximum number of Gauss-Newton iterations
% opt.prop.tol          - Lower bound tolerance criterion
% opt.prop.reg          - Initial regularisation parameter
% opt.prop.do           - Optimise tissue proportions?
%
% LINE SEARCH
% -----------
% opt.nline_search.bf   - Number of line-search iterations for bias field
% opt.nline_search.aff  - Number of line-search iterations for affine
% opt.nline_search.nl   - Number of line-search iterations for non-linear
% opt.nline_search.prop - Number of line-search iterations for proportions
%
% CLEANING
% --------
% opt.clean.brain               - Do ad-hock brain clean-up with clean_brain()
% opt.clean.mrf.strength        - Diagonal value of linear MRF
% opt.clean.mrf.niter           - Number of post MRF iterations
% opt.clean.les.bwlabeln        - Try to extract lesion by connected
%                                 components
% opt.clean.les.val             - Probablity threshold for lesion class
% opt.clean.les.class           - The tissue class containing lesion
% opt.clean.les.cnn_mrf.do      - Post-processing lesions using CNN-MRF
% opt.clean.les.cnn_mrf.pth_net - Path to CNN file on disk
%
% STARTING ITERATION
% ------------------
% opt.start_it.do_mg        - Iteration at when to introduce multiple
%                             Gaussians per tissue
% opt.start_it.do_prop      - Iteration at which to start optimising proportions
% opt.start_it.do_upd_mrf   - Iteration at which to start optimising MRF weights
% opt.start_it.upd_mg       - Iteration at which to start optimising
%                             Gaussian weights
%
% ACTIVATE/DEACTIVATE
% -------------------
% opt.do.mg             - Use more than one Gaussian per tissue
% opt.do.update_mrf     - Update MRF weights
% opt.do.mrf            - Use MRF prior in segmentation model
%
% CT-SPECIFIC
% -----------
% opt.ct.GaussPrior     - GMM prior for CT
%
% LESIONS
% -------
% opt.lesion.hemi       - Split lesion class between hemispheres
%
% FILES TO WRITE
% --------------
% opt.write.tc          - Write various segmentations:
%                         tc(1:K,1): Subject space responsibilities
%                         tc(1:K,2): Subject space responsibilities, post
%                         processed
%                         tc(1:K,3): Template space responsibilities, post
%                         processed
%                         tc(1:K,4): Modulated template space responsibilities, post
%                         processed
%                         tc(1:K,5): DARTEL imports
% opt.write.bf          - Write bias field:
%                         bf(1): Bias field corrected image in subject
%                         space
%                         bf(2): Bias field
%                         bf(3): Bias field corrected image in template
%                         space
% opt.write.df          - Write initial velocities?
% opt.write.ml          - Write maximum-likelihood labels?
% opt.write.les         - Write lesion masks?
% opt.dir_output_train  - Folder to write model (= population) data
% opt.dir_output_seg    - Folder to write segmentation (= subject) data
%
% DICTIONARIES
% ------------
% opt.dict.lkp          - Mapping between GMM clusters and Template classes
% opt.dict.prop_exc     - Classes to force responsibilities to zero, e.g, 
%                         prop_exc = [0 0 1 0], forces the responsibilities
%                         of class 3 to zero. Good when a population has
%                         lesions and one population does not.
%
% VERBOSITY
% ---------
% opt.verbose.level     - Set overall verbosity of algorithm:
%                         0: No verbose at all
%                         1: Only printed verbose
%                         2: Printed verbose + graphics
% opt.verbose.model     - Population-wise verbose, e.g., template, global
%                         tissue prior
% opt.verbose.gmm       - Subject GMM verbose
% opt.verbose.reg       - Subject registration verbose
% opt.verbose.bf        - Subject bias field verbose
% opt.verbose.prop      - Subject tissue proportion verbose
% opt.verbose.mrf       - Subject MRF verbose
% opt.verbose.mx_rows   - Number of rows for the sample figures
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 1, opt = struct; end

def = spm_shoot_defaults;

% opt
if ~isfield(opt,'dir_output')
    opt.dir_output = './output/';
end

% opt.model
if ~isfield(opt,'model') 
    opt.model         = struct;
end
if ~isfield(opt.model,'it') 
    opt.model.it      = 1;
end
if ~isfield(opt.model,'tol') 
    opt.model.tol     = 1e-4;
end
if ~isfield(opt.model,'niter') 
    opt.model.niter   = 20;
end
if ~isfield(opt.model,'nam_cls') 
    opt.model.nam_cls = {};
end
if ~isfield(opt.model,'clean_up') 
    opt.model.clean_up = true;
end

% opt.model.PropPrior
if ~isfield(opt.model,'PropPrior') 
    opt.model.PropPrior    = struct;
end
if ~isfield(opt.model.PropPrior,'do') 
    opt.model.PropPrior.do = true;
end

% opt.gmm
if ~isfield(opt,'gmm') 
    opt.gmm                = struct;
end
if ~isfield(opt.gmm,'niter') 
    opt.gmm.niter          = 20;
end
if ~isfield(opt.gmm,'tol') 
    opt.gmm.tol            = 1e-4;
end
if ~isfield(opt.gmm,'pth_GaussPrior') 
    opt.gmm.pth_GaussPrior = '';
end
if ~isfield(opt.gmm,'pth_PropPrior') 
    opt.gmm.pth_PropPrior  = '';
end

% opt.gmm.hist
if ~isfield(opt.gmm,'hist') 
    opt.gmm.hist             = struct;
end
if ~isfield(opt.gmm.hist,'niter_main') 
    opt.gmm.hist.niter_main  = 5;
end
if ~isfield(opt.gmm.hist,'niter_gmm') 
    opt.gmm.hist.niter_gmm   = 10; 
end
if ~isfield(opt.gmm.hist,'init_ix') 
    % For setting indices of classes (e.g. to match different modalities)
    map                      = containers.Map;
    opt.gmm.hist.init_ix     = map; 
end
if ~isfield(opt.gmm.hist,'verbose') 
    opt.gmm.hist.verbose     = true; % [true,false]
end
if ~isfield(opt.gmm.hist,'verbose_gmm') 
    opt.gmm.hist.verbose_gmm = 0; % [0,1,2,3]
end

% opt.gmm.GaussPrior
if ~isfield(opt.gmm,'GaussPrior') 
    opt.gmm.GaussPrior             = struct;
end
if ~isfield(opt.gmm.GaussPrior,'constrained') 
    opt.gmm.GaussPrior.constrained = false;
end
if ~isfield(opt.gmm.GaussPrior,'verbose') 
    opt.gmm.GaussPrior.verbose     = true; % [true,false]
end
if ~isfield(opt.gmm.GaussPrior,'mods') 
    opt.gmm.GaussPrior.mods        = {};
end
if ~iscell(opt.gmm.GaussPrior.mods)
    opt.gmm.GaussPrior.mods = {opt.gmm.GaussPrior.mods};
end

% opt.gmm.labels
if ~isfield(opt.gmm,'labels') 
    opt.gmm.labels     = struct;
end
if ~isfield(opt.gmm.labels,'cm') 
    % For including labels (if provided)
    map                = containers.Map;
    opt.gmm.labels.cm  = map;
end
if ~isfield(opt.gmm.labels,'use') 
    opt.gmm.labels.use = false;
end
if ~isfield(opt.gmm.labels,'S') 
    opt.gmm.labels.S   = 0.99;
end
if ~isfield(opt.gmm.labels,'Su') 
    opt.gmm.labels.Su  = 0.6;
end

% opt.template
if ~isfield(opt,'template') 
    opt.template              = struct;
end
if ~isfield(opt.template,'do')
    opt.template.do           = true;
end
if ~isfield(opt,'sched') 
    opt.sched = get_sched(opt);
end
if ~isfield(opt.template,'pth_template')
    opt.template.pth_template = '';
end
if ~isfield(opt.template,'K')
    if exist(opt.template.pth_template,'file')
        Nii                   = nifti(opt.template.pth_template);
        opt.template.K        = Nii.dat.dim(4);        
    else
        opt.template.K        = 6;
    end
end
if ~isfield(opt.template,'vs')
    opt.template.vs           = 1.5;
end
if numel(opt.template.vs) == 1
    opt.template.vs           = opt.template.vs*ones(1,3);
end
if ~isfield(opt.template,'reg0')
    opt.template.reg0         = def.sparam;
end
if ~isfield(opt.template,'reg')
    opt.template.reg          = [opt.template.reg0(1:2) opt.sched.a(1)*opt.template.reg0(3)];
end
if ~isfield(opt.template,'shrink')
    opt.template.shrink       = false;
end
if ~isfield(opt.template,'load_a_der')
    opt.template.load_a_der   = true;
end
if ~isfield(opt.template,'R')
    opt.template.R            = null(ones(1,opt.template.K));
end
if ~isfield(opt.template,'sym')
    opt.template.sym          = false;
end
if ~isfield(opt.template,'niter')
    opt.template.niter        = 16;
end
if ~isfield(opt.template,'verbose')
    opt.template.verbose      = 0; % [0,1,2]
end
if ~isfield(opt.template,'bg_class')
    opt.template.bg_class     = 0;
end
if ~isfield(opt.template,'resize')
    opt.template.resize       = true;
end
if ~isfield(opt.template,'keep_neck')
    opt.template.keep_neck    = false;
end
if ~isfield(opt.template,'clean')
    opt.template.clean        = struct;
end
if ~isfield(opt.template.clean,'its')
    opt.template.clean.its    = [];
end
if ~isfield(opt.template.clean,'brain')
    opt.template.clean.brain  = [];
end
if ~isfield(opt.template.clean,'les')
    opt.template.clean.les    = [];
end
if ~isfield(opt.template.clean,'air')
    opt.template.clean.air    = [];
end
if ~isfield(opt.template,'tc_miss')
    opt.template.tc_miss      = [];
end
if ~isfield(opt.template.clean,'val_brain')
    opt.template.clean.val_brain = 0.3;
end
if ~isfield(opt.template.clean,'val_air')
    opt.template.clean.val_air   = 0.3;
end
if ~isfield(opt.template.clean,'dil_er')
    opt.template.clean.dil_er   = false;
end
if ~isfield(opt.template.clean,'it_dil_er')
    opt.template.clean.it_dil_er   = 8;
end

% opt.reg
if ~isfield(opt,'reg') 
    opt.reg              = struct;
end
if ~isfield(opt.reg,'rparam0') 
    opt.reg.rparam0 = def.rparam;
%     opt.reg.rparam0      = [1e-4  1e-1 2 0.25 0.5]*0.01;%[0 0.005 0.2 0.025 0.05];
end
if ~isfield(opt.reg,'rparam') 
    % [absolute displacements, laplacian, bending energy, linear elasticity mu, linear elasticity lambda]
    opt.reg.rparam       = opt.sched.reg(1)*opt.reg.rparam0;
end
if ~isfield(opt.reg,'int_args') 
    opt.reg.int_args     = opt.sched.eul(1);
end
if ~isfield(opt.reg,'niter') 
    opt.reg.niter        = 1;
end
if ~isfield(opt.reg,'tol') 
    opt.reg.tol          = 1e-4;
end
if ~isfield(opt.reg,'strt_nl') 
    if opt.template.do
        opt.reg.strt_nl      = 4;
    else
        opt.reg.strt_nl      = 1;
    end
end
if ~isfield(opt.reg,'mc_aff') 
    opt.reg.mc_aff       = true;
end
if ~isfield(opt.reg,'aff_type') 
    opt.reg.aff_type     = 'similitude'; % ['translation','rotation','rigid','similitude','affine']
end
if ~isfield(opt.reg,'aff_reg') 
    opt.reg.aff_reg      = 0;
end
if ~isfield(opt.reg,'do_aff') 
    opt.reg.do_aff       = true;
end
if ~isfield(opt.reg,'do_nl') 
    opt.reg.do_nl        = true;
end
if ~isfield(opt.reg,'nit_init_aff') 
    opt.reg.nit_init_aff = 12;
end
if ~isfield(opt.reg,'init_aff_tol') 
    opt.reg.init_aff_tol = 1e-4;
end

% holly
if ~isfield(opt,'holly') 
    opt.holly          = struct;
end
if ~isfield(opt.holly,'mode') 
    if opt.template.do
        opt.holly.mode = 'parfor';
    else
        opt.holly.mode = 'for';
    end
end
holly                  = distribute_default(opt.holly);

% opt.seg
if ~isfield(opt,'seg') 
    opt.seg                = struct;
end
if ~isfield(opt.seg,'niter')
    opt.seg.niter          = 20;
end
if ~isfield(opt.seg,'tol')
    opt.seg.tol            = 1e-4;
end
if ~isfield(opt.seg,'show')
    opt.seg.show           = false;
end
if ~isfield(opt.seg,'samp')
    opt.seg.samp           = 0;
end
if ~isfield(opt.seg,'bg')
    opt.seg.bg             = 1;
end
if ~isfield(opt.seg,'mskonlynan')
    opt.seg.mskonlynan     = false;
end

% opt.seg.mrf
if ~isfield(opt.seg,'mrf') 
    opt.seg.mrf          = struct;
end
if ~isfield(opt.seg.mrf,'ml')
    opt.seg.mrf.ml       = true;
end
if ~isfield(opt.seg.mrf,'val_diag')
    opt.seg.mrf.val_diag = 0.5;
end
if ~isfield(opt.seg.mrf,'alpha')
    opt.seg.mrf.alpha    = 1e5;
end
if ~isfield(opt.seg.mrf,'niter')
    opt.seg.mrf.niter    = 1;
end

% opt.bf
if ~isfield(opt,'bf') 
    opt.bf               = struct;
end
if ~isfield(opt.bf,'biasfwhm')
    opt.bf.biasfwhm      = 60;
end
if ~isfield(opt.bf,'niter')
    opt.bf.niter         = 8;
end
if ~isfield(opt.bf,'tol')
    opt.bf.tol           = 1e-4;
end
if ~isfield(opt.bf,'mc_bf')
    opt.bf.mc_bf         = false;
end
if ~isfield(opt.bf,'biasreg')
    opt.bf.biasreg       = 1e4;
end
if ~isfield(opt.bf,'do')
    opt.bf.do            = true;
end
if ~isfield(opt.bf,'mc_bf_verbose')
    opt.bf.mc_bf_verbose = false; % [true,false]
end

% opt.prop
if ~isfield(opt,'prop') 
    opt.prop         = struct;
end
if ~isfield(opt.prop,'niter')
    opt.prop.niter   = 1;
end
if ~isfield(opt.prop,'gnniter')
    opt.prop.gnniter = 1;
end
if ~isfield(opt.prop,'tol')
    opt.prop.tol     = 1e-4;
end
if ~isfield(opt.prop,'reg')     
    opt.prop.reg     = 1;
end
if opt.prop.reg <= 0
    error('opt.prop.reg <= 0');
end
if ~isfield(opt.prop,'do')
    opt.prop.do      = true;
end

% opt.nline_search
if ~isfield(opt,'nline_search') 
    opt.nline_search      = struct;
end
if ~isfield(opt.nline_search,'bf')
    opt.nline_search.bf   = 6;
end
if ~isfield(opt.nline_search,'aff')
    opt.nline_search.aff  = 6;
end
if ~isfield(opt.nline_search,'nl')
    opt.nline_search.nl   = 6;
end
if ~isfield(opt.nline_search,'prop')
    opt.nline_search.prop = 6;
end

% opt.clean
if ~isfield(opt,'clean') 
    opt.clean          = struct;
end
if ~isfield(opt.clean,'brain')
    opt.clean.brain = true;
end

% opt.clean.mrf
if ~isfield(opt.clean,'mrf') 
    opt.clean.mrf          = struct;
end
if ~isfield(opt.clean.mrf,'strength')
    opt.clean.mrf.strength = 1;
end
if ~isfield(opt.clean.mrf,'niter')
    opt.clean.mrf.niter    = 10;
end

% opt.clean.les
if ~isfield(opt.clean,'les') 
    opt.clean.les          = struct;
end
if ~isfield(opt.clean.les,'bwlabeln')
    opt.clean.les.bwlabeln = false;
end
if ~isfield(opt.clean.les,'val')
    opt.clean.les.val      = 0.1;
end
if ~isfield(opt.clean.les,'class')
    opt.clean.les.class    = 0;
end

% opt.clean.les.cnn_mrf
if ~isfield(opt.clean.les,'cnn_mrf') 
    opt.clean.les.cnn_mrf         = struct;
end
if ~isfield(opt.clean.les.cnn_mrf,'do')
    opt.clean.les.cnn_mrf.do      = false;
end
if ~isfield(opt.clean.les.cnn_mrf,'pth_net')
    opt.clean.les.cnn_mrf.pth_net = '';
end

if (opt.clean.les.class <= 0 || opt.clean.les.class > opt.template.K) && (opt.clean.les.bwlabeln || opt.clean.les.cnn_mrf.do)
    error('Value error: opt.clean.les.class')
end

% opt.start_it
if ~isfield(opt,'start_it') 
    opt.start_it            = struct;
end
if ~isfield(opt.start_it,'do_mg')
    opt.start_it.do_mg      = 1;
end
if ~isfield(opt.start_it,'do_prop')
    opt.start_it.do_prop    = opt.reg.strt_nl + 1;
end
if ~isfield(opt.start_it,'do_upd_mrf')
    opt.start_it.do_upd_mrf = opt.reg.strt_nl + 1;
end
if ~isfield(opt.start_it,'upd_mg')
    opt.start_it.upd_mg = opt.reg.strt_nl + 1;
end

% opt.do
if ~isfield(opt,'do') 
    opt.do            = struct;
end
if ~isfield(opt.do,'mg')
    opt.do.mg         = true;
end
if ~isfield(opt.do,'update_mrf')
    opt.do.update_mrf = false;
end
if ~isfield(opt.do,'mrf')
    opt.do.mrf        = false;
end
if opt.do.mrf && opt.seg.samp > 0
    error('opt.do.mrf = true does not work at the moment, need to implement resizing of Z...')
end

% opt.ct
if ~isfield(opt,'ct') 
    opt.ct            = struct;
end
if ~isfield(opt.ct,'GaussPrior')
    opt.ct.GaussPrior = [];
end
if isempty(opt.ct.GaussPrior)
    lb_prW         = struct;                                        
    lb_prW.KL_qVpV = 0;
    lb_prW.ElnDetV = zeros(1,opt.template.K);

    b  = 1e4*ones(1,opt.template.K);
    n  = 0.04*ones(1,opt.template.K);    
    MU = [-995 -50 10 linspace(20,50,opt.template.K - 6) 70 100 500];
    W  = ones([1 1 opt.template.K]);

    opt.ct.GaussPrior = {MU,b,W,n,'CT',lb_prW,1:opt.template.K};   
end

% opt.lesion
if ~isfield(opt,'lesion') 
    opt.lesion      = struct;
end
if ~isfield(opt.lesion,'hemi')
    opt.lesion.hemi = [];
end

% opt.write
if ~isfield(opt,'write') 
    opt.write        = struct;
end
if ~isfield(opt.write,'tc')
    opt.write.tc     = true(opt.template.K,5); % [native-orig,native-final,template-final,modulated-template-final,dartel]
elseif size(opt.write.tc,2) ~=5
    opt.write.tc     = true(opt.template.K,5);
end
if ~isfield(opt.write,'bf')
    opt.write.bf     = true(1,3); % [native-im,bf,template-im]
end
if ~isfield(opt.write,'df')
    opt.write.df     = true;
end
if ~isfield(opt.write,'ml')
    opt.write.ml     = false;
end
if ~isfield(opt.write,'les')
    opt.write.les    = true(1,2);
end

% opt.dict
if ~isfield(opt,'dict') 
    opt.dict          = struct;
end
if ~isfield(opt.dict,'lkp')
    % For using multiple Gaussians per tissue
    map                = containers.Map;
    opt.dict.lkp       = map;
end
if ~isfield(opt.dict,'prop_excl')
    % For constraining a class to have zero resps by setting its proportion to
    % a large negative value
    map                = containers.Map;
    opt.dict.prop_excl = map;
end

% opt.verbose
if ~isfield(opt,'verbose') 
    opt.verbose         = struct;
end
if ~isfield(opt.verbose,'level') 
    opt.verbose.level   = 2;
end
if ~isfield(opt.verbose,'mx_rows') 
    opt.verbose.mx_rows = 8;
end
if ~isfield(opt.verbose,'model') 
    opt.verbose.model = 0; % [0,1,2,3]
end
if opt.verbose.level == 2
    opt.verbose.gmm   = 4; % [0,1,2,3,4,5]
    opt.verbose.reg   = 3; % [0,1,2,3]
    opt.verbose.bf    = 3; % [0,1,2,3]   
    opt.verbose.prop  = 3; % [0,1,2,3]   
    opt.verbose.mrf   = 3; % [0,1,2,3]  
elseif opt.verbose.level == 1
    opt.verbose.gmm   = 1;
    opt.verbose.reg   = 1;
    opt.verbose.bf    = 1;    
    opt.verbose.prop  = 1;
    opt.verbose.mrf   = 1;
else
    opt.verbose.gmm   = 0;
    opt.verbose.reg   = 0;
    opt.verbose.bf    = 0;    
    opt.verbose.prop  = 0;
    opt.verbose.mrf   = 0;    
end
if ~strcmpi(holly.mode,'for')
    opt.verbose.gmm  = 0;
    opt.verbose.reg  = 0;
    opt.verbose.bf   = 0;
    opt.verbose.prop = 0;
    opt.verbose.mrf  = 0;
end

if opt.template.do 
    % Some options when training                                        
    if opt.template.load_a_der
        opt.template.verbose = 1;        
    else        
        opt.template.verbose = 2;        
    end
    opt.verbose.model    = 3; 
    opt.bf.mc_bf_verbose = true;    
end

opt.dir_output_train = fullfile(opt.dir_output,'train');
opt.dir_output_seg   = fullfile(opt.dir_output,'subject-results');
%==========================================================================