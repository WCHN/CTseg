function [labels,mn,mx] = get_labels(dat,opt,samp,subsmp,grd)
% FORMAT [labels,mn,mx] = get_labels(dat,opt)
% dat    - Subject's data structure (one subject)
% opt    - Options structure
% grd    - Subsampling info struct
% grd    - Subsampling grid struct
% labels - {1} Image of labels (in uint8)
%          {2} Confusion matrix
% mn     - Minimum label value
% mx     - Maximum label value
%
% Load image of labels from disk + post-process + get confusion matrix
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 3, samp   = 0;  end
if nargin < 4, subsmp = []; end
if nargin < 5, grd    = []; end

mn     = 0;
mx     = 0;
labels = {}; 
if isfield(dat,'label') && opt.gmm.labels.use
           
    if ~opt.gmm.labels.cm.isKey(dat.population)
        return
    else
        ix = opt.gmm.labels.cm(dat.population);
    end        
    
    ix_bg  = numel(ix);
    dm     = dat.label{1}.nii.dat.dim;
    V      = spm_vol(dat.label{1}.nii.dat.fname);
    
    % We may need to sub-sample the labels
    if samp > 0
        dm     = subsmp.dm;
        dm     = [dm 1];
        labels = zeros(dm(1:3),'uint8');
        for z=1:numel(grd.z0)
            labels(:,:,z) = spm_sample_vol(V,grd.x0,grd.y0,grd.z0(z)*grd.o,0);
        end
        labels = labels(:);
    else
       labels = uint8(dat.label{1}.nii.dat(:)); 
    end            
%     figure(666); imshow3D(reshape(labels,dm))
    
    unq = [];
    for l=1:numel(ix) - 1
        if ~isempty(ix{l})
            unq = [unq l];
        end
    end
    
    msk          = ismember(labels,unq);
    labels(~msk) = ix_bg;
%     figure(666); imshow3D(reshape(msk,dm))
    
    mn = min(labels);
    mx = max(labels);
            
    CM = get_label_cm(dat,opt);
    
    labels = {labels,CM};   
end
%==========================================================================