function [dat,opt] = build_dir_structure(dat,opt)
% FORMAT [dat,opt] = build_dir_structure(dat,opt)
% dat   - Subjects data structure
% opt   - Options structure
% 
% Build directory structure on disk to store all temporary variables:
%
% <dir_output_train>
% |-- model
% |   |-- animations
% |-- vel
% |-- a_der
% |-- seg2d
%
% <dir_output_seg>
% |-- <subject>
% |   |-- images
% |   |-- segmentations
% |   |-- bias-field
% |   |-- registration
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if opt.template.do
    dir_output_train = opt.dir_output_train;
    if exist(dir_output_train,'dir'), rmdir(dir_output_train,'s'); end; mkdir(dir_output_train);   

    % Create directory that will store all that is model related
    dir_model     = fullfile(dir_output_train,'model');
    mkdir(dir_model);
    opt.dir_model = dir_model;

    % Create directory that will store all that animation related
    dir_anim           = fullfile(dir_model,'animations');
    mkdir(dir_anim);
    opt.dir_animations = dir_anim;
    
    % Create directory that will store all initial velocities
    dir_vel     = fullfile(dir_output_train,'vel');
    mkdir(dir_vel);   
    opt.dir_vel = dir_vel;   

    % Create directory that will store template derivatives
    dir_a_der     = fullfile(dir_output_train,'a_der');
    mkdir(dir_a_der);
    opt.dir_a_der = dir_a_der;
    
    if opt.verbose.model >= 3
        % Create directory that will store 2D segmentations
        dir_seg2d     = fullfile(dir_output_train,'seg2d');
        mkdir(dir_seg2d);
        opt.dir_seg2d = dir_seg2d;
    end
end

dir_output_seg = opt.dir_output_seg;
if ~exist(dir_output_seg,'dir'), mkdir(dir_output_seg); end;

S0 = numel(dat);
for s=1:S0
    nam = dat{s}.name;
    d   = fullfile(dir_output_seg,nam);
    
    if exist(d,'dir'), rmdir(d,'s'); end; mkdir(d);       

    % Create directory that stores segmentations
    dir_seg0 = fullfile(d,'segmentations');
    mkdir(dir_seg0);         

    if any(opt.write.tc(:,1) == true) || ((opt.clean.les.cnn_mrf.do || opt.clean.les.bwlabeln) && opt.write.les(1))
        % Create directory that stores segmentations before cleaning up
        dir_seg_orig        = dir_seg0;    
        dat{s}.dir.seg_orig = dir_seg_orig;   
    end
    
    if any(opt.write.tc(:,2) == true) || any(opt.write.tc(:,3) == true) || ((opt.clean.les.cnn_mrf.do || opt.clean.les.bwlabeln) && opt.write.les(2))
        % Create directory that stores final segmentations
        dir_seg        = dir_seg0;
        dat{s}.dir.seg = dir_seg;                
    end
    
    if opt.write.ml
        % Create directory that ML estimates from segmentations
        dir_ml     = dir_seg0;
        dat{s}.dir.ml = dir_ml;   
    end

    if opt.write.df
        % Create directory that stores initial velocities
        dir_def        = fullfile(d,'registration');
        mkdir(dir_def);
        dat{s}.dir.def = dir_def;   
    end
    
    if opt.write.bf(1) || opt.write.bf(3)
        % Create directory that stores images
        dir_img        = fullfile(d,'images');
        mkdir(dir_img);
        dat{s}.dir.img = dir_img;   
    end
    
    if opt.write.bf(2)
        % Create directory that stores bias-field
        dir_bf        = fullfile(d,'bias-field');
        mkdir(dir_bf);
        dat{s}.dir.bf = dir_bf; 
    end
    
    if S0 > 1 && ~opt.template.do
        % Create directory that stores temporary initial velocities
        dir_vel        = fullfile(d,'temp');
        mkdir(dir_vel);   
        dat{s}.dir.vel = dir_vel;   
    end
end
%==========================================================================