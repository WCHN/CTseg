function Z = get_resp(obs,bf,dat,template,labels,BinWidth,miss,dm,opt,varargin)
% FORMAT Z = get_resp(obs,bf,dat,template,labels,BinWidth,miss,dm,opt,['use_uint',...,'z',...])
% obs      - Observed image
% bf       - Reconstructed bias field
% dat      - Subject's data structure
% template - Warped and softmaxed template
% labels   - Manual labels
% BinWidth - Observed values precision (adds a bit of regularisation)
% miss     - Missing data structure
% dm       - Observed image dimensions
% opt      - Options structure
% use_uint - Store responsibilities as uint8 (default: single)
% z        - Get a specific slice (default: all)
% Z        - (Slice of) cluster responsibilities
%
% Compute (a slice of) responsibilities. Because the number of GMM clusters
% can be quite large, full responsibilities can't be held in memory. We
% thus recompute slices on the fly from its dependencies each time it is
% needed. This slows things down but saves tremendous amounts of memory.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Parse inputs
cluster = dat.gmm.cluster;
prop    = dat.gmm.prop;
part    = dat.gmm.part;

p              = inputParser;
p.FunctionName = 'get_resp';
p.addParameter('use_uint', false, @islogical);
p.addParameter('z',        0,     @isnumeric);
p.parse(varargin{:});
use_uint       = p.Results.use_uint;
z              = p.Results.z;

% Parameters
lkp     = part.lkp;
K       = max(lkp);
const   = spm_gmm_lib('Const', cluster{1}, cluster{2}, miss.L);
ix_tiny = get_par('ix_tiny',dat.population,part.lkp,opt);

% Get current responsibilities
if z > 0 && z <= dm(3)       
    %----------------------------------------------------------------------
    % Get responsibilities for a slice
    %----------------------------------------------------------------------
        
    dm(3) = 1;
     
    % Neighborhood part
    lnPzN  = gmm_mrf('apply',dat.mrf,z);
    lnPzNz = double(lnPzN);  
    
    % Get slice data
    slice = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels,BinWidth);
   
    % Get responsibilities for a slice
    Z_slice = gmm_img('slice_resp_and_lb',slice,cluster{1},cluster{2},prop,part,miss,const,lnPzNz,ix_tiny);

    % Go from cluster to tissue responsibilities
    Z_slice = cluster2template(Z_slice,part);    

    if use_uint
        Z_slice = uint8((2^8)*Z_slice);
    end

    Z = reshape(Z_slice,[dm(1:2) 1 K]);
else
    %----------------------------------------------------------------------
    % Get full responsibilities
    %----------------------------------------------------------------------
        
    % Neighborhood part
    lnPzN = gmm_mrf('apply',dat.mrf);

    if use_uint
        Z = zeros([dm K],'uint8');
    else
        Z = zeros([dm K],'single');
    end

    for z=1:dm(3)        
        % Get slice data
        [slice,ix] = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels,BinWidth);

        if dat.mrf.do && numel(lnPzN) > K    
            lnPzNz = double(lnPzN(ix,:));
        else
            lnPzNz = lnPzN;
        end

        % Get responsibilities for a slice
        Z_slice = gmm_img('slice_resp_and_lb',slice,cluster{1},cluster{2},prop,part,miss,const,lnPzNz,ix_tiny);

        % Go from cluster to tissue responsibilities
        Z_slice = cluster2template(Z_slice,part);    

        if use_uint
            Z_slice = uint8((2^8)*Z_slice);
        end

        Z(:,:,z,:) = reshape(Z_slice,[dm(1:2) 1 K]);
    end
end
%==========================================================================