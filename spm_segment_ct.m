function results = spm_segment_ct(Image,DirOut,CleanBrain,Write,Samp,MRF)
% Combined segmentation and spatial normalisation of computed tomography (CT) images
%
% FORMAT results = spm_segment_ct(Image,DirOut,CleanBrain,Write,Samp,MRF)
%
% INPUT
% PthImage   - A CT image, given as:
%                1. A nifti object (nifti(pth))             
%                2. A string with a nifti filename
%                3. The path to a folder with DICOM files, corresponding to 
%                   one CT image
% DirOut     - The directory where to write all of the results ['CTseg-Results']
% CleanBrain - Run an ad-hoc brain clean-up routine [false]
% Write      - What results to write to DirOut:
%                Write.image  = [native image, warped image]          [1 0]
%                Write.native = [native seg, dartel seg]              [1 0]
%                Write.warped = [] [warped seg, warped modulated seg] [0 0]
% Samp       - Sub-sampling distance, to speed things up [3]
% MRF        - Run post-processing MRF [1]
%
% OUTPUT
% results  - A struct with paths to all results
%_______________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Image = '/home/mbrud/dev/packages/matlab/spm/trunk/toolbox/CTseg/data/PLORAS-lesion-cerebellum.nii';
% Image = '/data/mbrud/populations/original/ATLAS-NOLABELS/c0004s0007t01.nii';
% Image = '/data/mbrud/populations/original/CROMIS/sCROMIS2ICH_26003-0002-00001-000001.nii';
% Image = '/data/mbrud/populations/original/CROMIS/sCROMIS2ICH_24036-0005-00003-000001.nii';
% Image = '/data/mbrud/populations/original/DELIRIUM/780_s99021516-0003-00002-000001.nii';
if nargin < 1
    Image = nifti(spm_select(1,'nifti','Select CT image'));   
end
if nargin < 2, DirOut     = 'CTseg-Results'; end
if nargin < 3, CleanBrain = false; end
if nargin < 4
Write        = struct;
Write.image  = [1 0];
Write.native = [1 0];
Write.warped = [0 0];
end
if nargin < 5, Samp = 3; end
if nargin < 6, MRF  = 1; end

% Some parameters
VoxSize    = [];     % Set VoxSize = [], to work in native resolution
DoDenoise  = false;
VerboseDen = 1;
VerboseSeg = 1;
DoPreproc  = true;

%--------------------------------------------------------------------------
% Add required toolboxes to path, and see if model file exists (if not d/l)
%--------------------------------------------------------------------------

spm_check_path('Shoot','Longitudinal','pull');

PthToolboxes = add2path;

get_model;
            
%--------------------------------------------------------------------------
% Read image to NIfTI
%--------------------------------------------------------------------------

Nii = get_nii(Image,DirOut);

% Copy so to not overwrite originals
Nii = make_copies(Nii,DirOut);

%--------------------------------------------------------------------------
% Preprocess
%--------------------------------------------------------------------------

if DoPreproc
    % Reset origin, and set voxels smaller than VoxSize to VoxSize.
    Nii = reset_origin(Nii,VoxSize);
    
    % Realing to MNI space
    Nii = realign2mni(Nii);

    % Crop air
    Nii = crop(Nii);

    if DoDenoise
        % Denoise
        Nii = denoise(Nii,'CT',VerboseDen);
    end
end

%--------------------------------------------------------------------------
% Segment
%--------------------------------------------------------------------------

opt = segment_ct(Nii,DirOut,PthToolboxes,VerboseSeg,CleanBrain,Write,Samp,MRF);

%--------------------------------------------------------------------------
% Clean-up
%--------------------------------------------------------------------------

results = clean_up(Nii,DirOut,opt);

if 1
    spm_check_registration(char({results.i{1},results.c{:}}))
end

return
%==========================================================================

%==========================================================================
function results = clean_up(Nii,DirOut,opt)

% Re-organise files
delete(Nii{1}(1).dat.fname);
[~,nam]  = fileparts(Nii{1}(1).dat.fname);
dir_res  = fullfile(opt.dir_output_seg,nam);
ndir_res = fullfile(DirOut,nam);
copyfile(dir_res,ndir_res);
rmdir(opt.dir_output_seg,'s')

% Create output struct
results     = struct;
results.c   = {};
results.rc  = {};
results.wc  = {};
results.mwc = {};
results.i   = {};
results.bf  = {};
results.wi  = {};
results.def = {};

if any(opt.write.tc(:,2))
    % c
    files = spm_select('FPListRec',ndir_res,'^c.*\.nii$');
    for k=1:size(files,1)
        results.c{k} = deblank(files(k,:));
    end
end

if any(opt.write.tc(:,2))
    % rc
    files = spm_select('FPListRec',ndir_res,'^rc.*\.nii$');
    for k=1:size(files,1)
        results.rc{k} = deblank(files(k,:));
    end
end

if any(opt.write.tc(:,2))
    % wc
    files = spm_select('FPListRec',ndir_res,'^wc.*\.nii$');
    for k=1:size(files,1)
        results.wc{k} = deblank(files(k,:));
    end
end

if any(opt.write.tc(:,2))
    % mwc
    files = spm_select('FPListRec',ndir_res,'^mwc.*\.nii$');
    for k=1:size(files,1)
        results.mwc{k} = deblank(files(k,:));
    end
end

if any(opt.write.tc(:,2))
    % i
    files = spm_select('FPListRec',ndir_res,'^i.*\.nii$');
    for k=1:size(files,1)
        results.i{k} = deblank(files(k,:));
    end
end

if any(opt.write.tc(:,2))
    % bf
    files = spm_select('FPListRec',ndir_res,'^bf.*\.nii$');
    for k=1:size(files,1)
        results.bf{k} = deblank(files(k,:));
    end
end

if any(opt.write.tc(:,2))
    % wi
    files = spm_select('FPListRec',ndir_res,'^wi.*\.nii$');
    for k=1:size(files,1)
        results.wi{k} = deblank(files(k,:));
    end
end

if any(opt.write.tc(:,2))
    % def
    files = spm_select('FPListRec',ndir_res,'^def.*\.nii$');
    for k=1:size(files,1)
        results.def{k} = deblank(files(k,:));
    end
end

return
%==========================================================================

%==========================================================================
function opt = segment_ct(Nii,DirOut,PthToolboxes,VerboseSeg,CleanBrain,Write,Samp,MRF)

if iscell(Nii)
    Nii = Nii{1};
end

isCT = min(Nii.dat(:)) < 0;

[~,nam] = fileparts(Nii.dat.fname);
d       = fullfile(DirOut,nam);
if exist(d,'dir'), rmdir(d,'s'); end

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------

dat    = cell(1);
dat{1} = struct;

[~,nam]     = fileparts(Nii.dat.fname);
dat{1}.name = nam;
if isCT
    dat{1}.modality{1}.nii  = Nii;
    dat{1}.modality{1}.name = 'CT';
    dat{1}.population       = 'CROMIS-LABELS';
else
    dat{1}.modality{1}.name            = 'MRI';
    dat{1}.modality{1}.channel{1}.name = 'T1';
    dat{1}.modality{1}.channel{1}.nii  = Nii;
    dat{1}.population                  = 'ATLAS-LABELS';
end

%--------------------------------------------------------------------------
% Model
%--------------------------------------------------------------------------

DirModel 	  = fullfile(spm('dir'),'toolbox','CTseg','model');
PthTemplate   = fullfile(DirModel,'template.nii');
PthGaussPrior = fullfile(DirModel,'GaussPrior.mat');
PthPropPrior  = fullfile(DirModel,'PropPrior.mat');
    
%--------------------------------------------------------------------------
% Options
%--------------------------------------------------------------------------

opt             = struct;
opt.dir_output  = DirOut;
opt.template.do = false;

opt.template.pth_template = PthTemplate;
opt.gmm.pth_GaussPrior    = PthGaussPrior;
% opt.gmm.pth_PropPrior     = PthPropPrior;

opt.verbose.level       = VerboseSeg;
opt.seg.samp            = Samp;
opt.prop.do             = true;
opt.start_it.upd_mg     = 2;
opt.start_it.do_prop    = 2;
opt.start_it.do_upd_mrf = 2;
opt.do.mrf              = false;
opt.do.update_mrf       = false;
opt.seg.mrf.val_diag    = 0.6;
opt.model.clean_up      = false;
opt.write.df            = false;

% Write results
opt.write.tc        = false(8,5);
if Write.native(1)
    opt.write.tc(:,2) = true;
end
if Write.native(2)
    opt.write.tc(:,5) = true;
end
if Write.warped(1)
    opt.write.tc(:,3) = true;
end
if Write.warped(2)
    opt.write.tc(:,4) = true;
end

opt.write.bf = false(1,3);
if Write.image(1)
    opt.write.bf(1) = true;
end
if Write.image(2)
    opt.write.bf(3) = true;
end

% Resps cleaning options
opt.clean.mrf.strength = MRF;

opt.clean.brain              = CleanBrain;
opt.template.clean.brain     = [4 5];
opt.template.clean.les       = [3 7];
opt.template.clean.val_brain = 0.5;
opt.template.clean.val_air   = 0.5;
opt.template.clean.dil_er    = true;
opt.template.clean.it_dil_er = 4;

% For using multiple Gaussians per tissue
map          = containers.Map;
map('CT')    = [1 1 1 2 2 2 3 4 5 6 7 8 8 8];
map('MRI')   = repelem(1:8,2);
opt.dict.lkp = map;

% These two are mandatory (for now)
opt.dep.aux_toolbox  = PthToolboxes{2};
opt.dep.dist_toolbox = PthToolboxes{3};
    
%--------------------------------------------------------------------------
% Segment
%--------------------------------------------------------------------------

opt = SegModel('segment',dat,opt);
%==========================================================================

%==========================================================================
function PthToolboxes = add2path

% DirCode      = './code';
% PthToolboxes = {'/home/mbrud/dev/mbrud/code/matlab/preprocessing-code', ...
%                 '/home/mbrud/dev/mbrud/code/matlab/auxiliary-functions', ...
%                 '/home/mbrud/dev/mbrud/code/matlab/distributed-computing', ...
%                 '/home/mbrud/dev/mbrud/code/matlab/MTV-preproc', ...
%                 '/home/mbrud/dev/mbrud/code/matlab/segmentation-model'};
PthToolboxes = {fullfile('toolboxes','preprocessing-code'), ...
                fullfile('toolboxes','auxiliary-functions'), ...
                fullfile('toolboxes','distributed-computing'), ...
                fullfile('toolboxes','MTV-preproc'), ...
                fullfile('toolboxes','segmentation-model')};
                        
% if (exist(DirCode,'dir') == 7)  
%    rmdir(DirCode,'s') ;
% end
% mkdir(DirCode);  
for i1=1:numel(PthToolboxes)            

    addpath(genpath(PthToolboxes{i1}));
%     nam = strsplit(PthToolboxes{i1},filesep);
%     nam = nam{end};
%     trg = fullfile(DirCode,nam);
%     mkdir(trg);
% 
%     copyfile(PthToolboxes{i1},trg);
% 
% %     res = dir(PthToolboxes{i1});
% %     for i2=1:numel(res)
% %         if strcmp(res(i2).name(1),'.')
% %             continue; 
% %         end        
% %         [~,~,ext] = fileparts(res(i2).name);
% %         if ~strcmp(ext,'.m')
% %             continue; 
% %         end
% %         
% %         src = fullfile(res(i2).folder,res(i2).name);        
% %         copyfile(src,trg);
% %     end
%     
%     % Add to path
%     addpath(genpath(trg));
end
%==========================================================================

%==========================================================================
function spm_check_path(varargin)
% Check that SPM is on the MATLAB path. Can also check to see if some other
% tools are avaiable.
%
% EXAMPLE USAGE
%
% Call as spm_check_path('Shoot','Longitudinal','pull') to check if these
% toolboxes and/or functions are available. If you just want to check if,
% for example, the Shoot toolbox is available, just do spm_check_path('Shoot').
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
   
% Check that SPM is on the MATLAB path
if ~(exist('spm','file') == 2)
    error('SPM is not on the MATLAB path, get the latest SPM version from: https://www.fil.ion.ucl.ac.uk/spm/software/ and add it to the MATLAB path'); 
end

if ~isempty(varargin)
    
    if any(strcmpi(varargin,'Shoot'))
        % Check that the Shoot toolbox is on the MATLAB path
        
        try
            spm_shoot3d;
        catch e
            if strcmp(e.message,'Undefined function or variable ''spm_shoot3d''.')
                error('Add the Shoot toolbox, from /spm/toolbox/Shoot, to the MATLAB path')
            end
        end
    end

    if any(strcmpi(varargin,'Longitudinal'))
        % Check that the Longitudinal toolbox is on the MATLAB path
        
        try
            spm_groupwise_ls;
        catch e
            if strcmp(e.message,'Undefined function or variable ''spm_groupwise_ls''.')
                error('Add the Longitudinal toolbox, from /spm/toolbox/Longitudinal, to the MATLAB path')
            end
        end
    end
    
    if any(strcmpi(varargin,'pull'))
        % Check that spm_diffeo('pull') is available
        
        try
            spm_diffeo('pull');
        catch e
            if strcmp(e.message,'Option not recognised.')
                error('The function spm_diffeo(''pull'') is not available, update to the latest SPM version from: https://www.fil.ion.ucl.ac.uk/spm/software/')
            end
        end
    end
end
%==========================================================================

%==========================================================================
function get_model
DirModel = 'model';
if ~(exist(DirModel,'dir') == 7)  
    mkdir(DirModel);  
end
n = 'template.nii';
f = fullfile(fileparts(mfilename('fullpath')),'model',n);
if ~isfile(f)
    fprintf('Downloading %s...',n);
    urlwrite('https://ndownloader.figshare.com/files/15103274',f);
    fprintf('done!\n');
end
n = 'GaussPrior.mat';
f = fullfile(fileparts(mfilename('fullpath')),'model',n);
if ~isfile(f)
    fprintf('Downloading %s...',n);
    urlwrite('https://ndownloader.figshare.com/files/15103268',f);
    fprintf('done!\n');    
end
n = 'PropPrior.mat';
f = fullfile(fileparts(mfilename('fullpath')),'model',n);
if ~isfile(f)
    fprintf('Downloading %s...',n);
    urlwrite('https://ndownloader.figshare.com/files/15103271',f);
    fprintf('done!\n');    
end
%==========================================================================