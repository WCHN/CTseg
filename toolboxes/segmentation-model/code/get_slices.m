function [slices,miss,lblnDetbf] = get_slices(dat,model,opt,do,reg)
% FORMAT [slices,miss,lblnDetbf] = get_slices(dat,model,opt,do,reg)
% dat       - Subject's data structure (one subject)
% model     - Model structure
% opt       - Options structure
% do        - Do structure with field bf 
% reg       - Registration structure obtained from `get_reg`
% slices    - Structure array with `nz` elements and fields
%             * obs      - One slice of observations
%             * bf       - One slice of bias fields
%             * template - One slice of warped + softmaxed template
%             * code     - One slice of code image
%             * bin_var  - Bin variance
%             * labels   - One slice of manual labels
% miss      - Missing data structure
% lblnDetbf - Bias field normalisation constant
%
% Split all useful arrays into slices.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

modality = dat.modality{1}.name; % Imaging modality

%--------------------------------------------------------------------------
% Get image data
%--------------------------------------------------------------------------

[obs,dm,~,~,scl] = get_obs(dat,'mskonlynan',opt.seg.mskonlynan,'samp',opt.seg.samp);

%--------------------------------------------------------------------------
% Missing data struct
%--------------------------------------------------------------------------

miss = get_par('missing_struct',obs);

%--------------------------------------------------------------------------
% Bias field
%--------------------------------------------------------------------------

if do.bf
    % Compute bias-field
    bf = get_bf(dat.bf.chan,dm);      
    
    lblnDetbf                   = bf;
    lblnDetbf(isnan(lblnDetbf)) = 1;        
    lblnDetbf                   = log(prod(lblnDetbf,2));  
    lblnDetbf                   = sum(lblnDetbf);
else  
    % Bias field not estimated
    bf = 1;    
    
    lblnDetbf = 0;
end

%--------------------------------------------------------------------------
% Labels (if provided)
%--------------------------------------------------------------------------

labels = get_labels(dat,opt); 

%--------------------------------------------------------------------------
% Warped template
%--------------------------------------------------------------------------

template = warp_template(model,reg.y,reg.Affine);  

%--------------------------------------------------------------------------
% Build struct holding slice data
%--------------------------------------------------------------------------

cl     = cell(dm(3),1);
slices = struct('obs',cl,'bf',cl,'template',cl,'code',cl,'bin_var',cl,'labels',cl);

for z=1:dm(3)
    slices(z) = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels,scl);
end
%==========================================================================