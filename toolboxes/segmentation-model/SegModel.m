function varargout = SegModel(varargin)
%__________________________________________________________________________
%
% Code for segmentation model
%
%--------------------------------------------------------------------------
%
% opt = SegModel('train',dat,opt)
% > Train segmentation model from a bunch of images
%
% opt = SegModel_continue(pth_dat,pth_model,pth_opt,pth_holly)
% > Continue training segmentation model from a bunch of images
%
% res = SegModel('segment',dat,opt)
% > Segment images with trained segmentation model
%
% res = SegModel('write',pth_dat,pth_model,pth_opt,pth_holly)
% > Segment images with trained segmentation model
%
% Arguments
% ----------
% dat - [1xS] Cell array with subject data
% opt -       Model options
% res - [1xS] Cell array with segmentation results
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin == 0
    help run_SegModel
    error('Not enough argument. Type ''help spm_gmm_lib'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'train'
        [varargout{1:nargout}] = SegModel_train(varargin{:});        
    case 'continue'
        [varargout{1:nargout}] = SegModel_continue(varargin{:});         
    case 'segment'
        [varargout{1:nargout}] = SegModel_segment(varargin{:});           
    case 'write'
        [varargout{1:nargout}] = SegModel_write(varargin{:});         
    otherwise
        help run_SegModel
        error('Unknown function %s. Type ''help run_SegModel'' for help.', id)
end
%==========================================================================

%==========================================================================
function opt = SegModel_train(dat,opt)
% Train segmentation model from a bunch of images
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

%--------------------------------------------------------------------------
% Init paths, etc.
%--------------------------------------------------------------------------

dat = SegModel_init(dat,opt);

%--------------------------------------------------------------------------
% Parse options
%--------------------------------------------------------------------------

[opt,holly] = default_opt(opt);

%--------------------------------------------------------------------------
% Initialise model
%--------------------------------------------------------------------------

make_animations('clear');

[dat,model,opt] = init_all(dat,opt);
model.lb        = -Inf;

%--------------------------------------------------------------------------
% Train model
%--------------------------------------------------------------------------

save_stuff(opt.dir_model,dat,opt,model,holly);

for it_mod=1:opt.model.niter                            
    
    [dat,model,opt,holly] = SegModel_iter(dat,model,opt,holly,it_mod);
    
end
%==========================================================================

%==========================================================================
function opt = SegModel_continue(pth_dat,pth_model,pth_opt,pth_holly,mem_new)
% Continue training segmentation model from a bunch of images
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

%--------------------------------------------------------------------------
% Load variables (dat, model, opt, holly)
%--------------------------------------------------------------------------

if nargin < 5, mem_new = []; end

var   = load(pth_dat);
dat   = var.dat;
var   = load(pth_model);
model = var.model;
var   = load(pth_opt);
opt   = var.opt;
var   = load(pth_holly);
holly = var.holly;
clear var

if ~isempty(mem_new)
    holly.job.mem{1} = mem_new;
end

%--------------------------------------------------------------------------
% Init paths, etc.
%--------------------------------------------------------------------------

SegModel_init(dat,opt);

%--------------------------------------------------------------------------
% Continue training model
%--------------------------------------------------------------------------

for it_mod=(opt.model.it + 1):opt.model.niter   
    
    [dat,model,opt,holly] = SegModel_iter(dat,model,opt,holly,it_mod);
    
end
%==========================================================================

%==========================================================================
function opt = SegModel_write(pth_dat,pth_model,pth_opt,pth_holly)
% Continue training segmentation model from a bunch of images
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

%--------------------------------------------------------------------------
% Load variables (dat, model, opt, holly)
%--------------------------------------------------------------------------

var   = load(pth_dat);
dat   = var.dat;
var   = load(pth_model);
model = var.model;
var   = load(pth_opt);
opt   = var.opt;
var   = load(pth_holly);
holly = var.holly;
clear var

%--------------------------------------------------------------------------
% Init paths, etc.
%--------------------------------------------------------------------------

SegModel_init(dat,opt);

%--------------------------------------------------------------------------
% Write segmentation results
%--------------------------------------------------------------------------

[~,~] = distribute(holly,'write_results','inplace',dat,model,opt);
%==========================================================================

%==========================================================================
function opt = SegModel_segment(dat,opt)
% Segment images with trained segmentation model
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

%--------------------------------------------------------------------------
% Init paths, etc.
%--------------------------------------------------------------------------

dat = SegModel_init(dat,opt);

%--------------------------------------------------------------------------
% Parse options
%--------------------------------------------------------------------------

[opt,holly] = default_opt(opt);

%--------------------------------------------------------------------------
% Initialise model
%--------------------------------------------------------------------------

[dat,model,opt] = init_all(dat,opt);

%--------------------------------------------------------------------------
% Segment subject(s)
%--------------------------------------------------------------------------
               
[~,dat] = distribute(holly,'segment_subject','inplace',dat,'iter',model,opt);  
[~,~]   = distribute(holly,'write_results','inplace',dat,'iter',model,opt);      
%==========================================================================

%==========================================================================
function [dat,model,opt,holly] = SegModel_iter(dat,model,opt,holly,it_mod)
% Run one iteration of segmentation model
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

opt.model.it = it_mod;

% Make animation of algorithm progress
make_animations('savefigs',opt,it_mod);

% Set-up loading template derivatives from disk       
dat = init_load_a_der(dat,opt);

% Introduce using multiple Gaussians per tissue
[dat,model] = introduce_lkp(dat,model,opt,it_mod);

% Do segmentation
[holly,dat] = distribute(holly,'segment_subject','inplace',dat,model,opt);      

% Some model parameters changes with iteration number
opt = modify_opt(opt,it_mod);

% Mean correct the rigid-body transforms
dat = meancorrect_aff(dat,opt);

% Mean correct bias field (also updates posteriors)
dat = meancorrect_bf(dat,model.GaussPrior,opt);          

% Update template	    
model = update_template(dat,model,opt);           

if ~isempty(opt.template.clean.its) && ismember(it_mod,opt.template.clean.its)
    % Ad-hoc template clean-up routine
    model = clean_template(model,opt,true);
end

% Update Gauss-Wishart hyper-parameters	    
model = update_GaussPrior(dat,model,opt); 

% Update proportions hyper-parameter
model = update_PropPrior(dat,model,opt,it_mod);

% Compute model lower bound
model = compute_model_lb(dat,model,opt);     

% Check convergence
gain = spm_misc('get_gain',model.lb);
if opt.verbose.model >= 1
    fprintf('%3s | %3d | lb = %10.6f | gain = %1.5f\n', 'mod', it_mod, model.lb(end), gain);
end        

% Some verbose
if opt.verbose.model >= 2, plot_model_lb(dat,model,it_mod,opt); end
if opt.verbose.model >= 3, show_tissues(dat,model,opt); end
if opt.verbose.model >= 3, show_params(dat,model,opt); end

% Save some variables
save_stuff(opt.dir_model,dat,opt,model,holly);

if it_mod == opt.model.niter 
    
    % Write segmentation results
    [~,~] = distribute(holly,'write_results','inplace',dat,model,opt);           

    % Make animation of algorithm progress
    make_animations('savefigs',opt,it_mod);
    make_animations('make',opt,it_mod);
    
    if opt.model.clean_up
        % Clean-up temporary files
        rmdir(opt.dir_vel,'s');
        rmdir(opt.dir_a_der,'s');
        rmdir(opt.dir_seg2d,'s');
    end
end
%==========================================================================

%==========================================================================
function dat = SegModel_init(dat,opt)
% Initialise segmentation model
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% ---
% SPM
if ~isfield(opt.dep, 'spm') || isempty(opt.dep.spm)
    opt.dep.spm = fileparts(which('spm'));
end
if (exist(opt.dep.spm, 'dir')==7)
    if ~isdeployed
        addpath(opt.dep.spm);
        addpath(fullfile(opt.dep.spm,'toolbox','Shoot'));
        addpath(fullfile(opt.dep.spm,'toolbox','Longitudinal'));
    end
end
if ~(exist('spm','file') == 2)
    error('SPM is not on the MATLAB path.'); 
end

% -------------------
% auxiliary-functions
if isfield(opt.dep,'aux_toolbox') && (exist(opt.dep.aux_toolbox, 'dir')==7)
    if ~isdeployed, addpath(opt.dep.aux_toolbox);  end
end
if ~(exist('spm_gmm','file') == 2)
    error('The "auxiliary" toolbox is not on the MATLAB path.'); 
end

% ------------------
% distribute-toolbox
if isfield(opt.dep,'dist_toolbox') && (exist(opt.dep.dist_toolbox, 'dir')==7)
    if ~isdeployed, addpath(opt.dep.dist_toolbox); end
end
if ~(exist('distribute','file') == 2)
    error('The "distribute" toolbox is not on the MATLAB path.'); 
end

% -------
% CNN-MRF
if isfield(opt,'clean') && isfield(opt.clean,'cnn_mrf') && isfield(opt.clean.cnn_mrf,'do')
    if opt.clean.cnn_mrf.do
        if ~isfield(opt.dep,'cnn_mrf')
            error('~isfield(opt.dep,''cnn_mrf'')')
        else
            if ~isdeployed, addpath(opt.dep.cnn_mrf); end
        end
    end
end

% ------------
% subfunctions
if ~isdeployed, addpath(fullfile(fileparts(mfilename('fullpath')),'code')); end

% -------
% globals
set_globals;

% --------------
% data structure
if ischar(dat)
    dat = spm_json_manager('init_dat',dat);
    if ~(isfield(opt,'template') && isfield(opt.template,'do') && opt.template.do)
        dat = dat(1);
    end
end
%==========================================================================

%==========================================================================
function save_stuff(dir_save,dat,opt,model,holly)
fname = fullfile(dir_save,'dat.mat');
save(fname,'dat','-v7.3');    
fname = fullfile(dir_save,'opt.mat');
save(fname,'opt','-v7.3');
fname = fullfile(dir_save,'model.mat');
save(fname,'model','-v7.3');
fname = fullfile(dir_save,'holly.mat');
save(fname,'holly','-v7.3');
%==========================================================================