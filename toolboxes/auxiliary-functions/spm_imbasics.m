function varargout = spm_imbasics(varargin)
%__________________________________________________________________________
% Collection of tools for image calculation (gradient, suff stat, ...).
%
% FORMAT [V,W,C]    = spm_imbasics('hist',X,...)
% FORMAT div        = spm_imbasics('dive',Dx,Dy,Dz,vx)
% FORMAT [Dx,Dy,Dz] = spm_imbasics('grad',X,vx)
% FORMAT spm_imbasics('smooth_img_in_mem',img,fwhm) 
% FORMAT nfname1 = spm_imbasics('create_2d_slice',fname,axis_2d,clean_up)
% FORMAT [BB,vx] = spm_imbasics('compute_bb',img,mat,dm,thr,premul)
%
% FORMAT help spm_imbasics>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging
if nargin == 0
    help spm_imbasics
    error('Not enough argument. Type ''help spm_imcalc'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'dive'
        [varargout{1:nargout}] = dive(varargin{:});        
    case 'grad'
        [varargout{1:nargout}] = grad(varargin{:});
    case 'smooth_img'
        [varargout{1:nargout}] = smooth_img_in_mem(varargin{:});           
    case 'hist'
        [varargout{1:nargout}] = spm_hist(varargin{:});            
    case 'create_2d_slice'
        [varargout{1:nargout}] = create_2d_slice(varargin{:});        
    case 'compute_bb'
        [varargout{1:nargout}] = compute_bb(varargin{:});          
    otherwise
        help spm_imcalc
        error('Unknown function %s. Type ''help spm_imcalc'' for help.', id)
end
%==========================================================================

%==========================================================================
function div = dive(Dx,Dy,Dz,vx)  
% Computes the divergence of an image (with voxel size)
% FORMAT div = dive(Dx,Dy,Dz,vx) 
% [Dx,Dy,Dz] - Gradients in x-,y- and z-direction
% vx         - Voxel size
% div        - Divergence
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin < 4, vx = ones([1 3],'like',Dx); end

if size(Dx,3) == 1
    Du  = [-Dx(:,1), -diff(Dx(:,1:end-1),1,2), Dx(:,end-1)];
    Dv  = [-Dy(1,:); -diff(Dy(1:end-1,:),1,1); Dy(end-1,:)];
    div = Du./vx(1) + Dv./vx(2);
else
    Du  = cat(2, -Dx(:,1,:), -diff(Dx(:,1:end-1,:),1,2), Dx(:,end-1,:)); 
    Dv  = cat(1, -Dy(1,:,:), -diff(Dy(1:end-1,:,:),1,1), Dy(end-1,:,:));
    Dw  = cat(3, -Dz(:,:,1), -diff(Dz(:,:,1:end-1),1,3), Dz(:,:,end-1));
    div = Du./vx(1) + Dv./vx(2) + Dw./vx(3);
end
%==========================================================================

%==========================================================================
function [Dx,Dy,Dz] = grad(X,vx) 
% Calculate 2D or 3D gradient of an image (with voxel size)
% FORMAT [Dx,Dy,Dz] = grad(X,vx)
% X          - Image
% vx         - voxel size
% [Dx,Dy,Dz] - Gradients in x-,y- and z-direction
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging   
if nargin < 2, vx = ones([1 3],'like',X); end

if size(X,3)==1
    Dx = [diff(X,1,2),zeros(size(X,1),1,'like',X)]./vx(1);
    Dy = [diff(X,1,1);zeros(1,size(X,2),'like',X)]./vx(2);
    Dz = 0;
else
    Dx = cat(2,diff(X,1,2),zeros(size(X,1),1,size(X,3),'like',X))./vx(1);
    Dy = cat(1,diff(X,1,1),zeros(1,size(X,2),size(X,3),'like',X))./vx(2);
    Dz = cat(3,diff(X,1,3),zeros(size(X,1),size(X,2),1,'like',X))./vx(3);  
end
%==========================================================================

%==========================================================================
function simg = smooth_img_in_mem(img,fwhm,VoxelSize) 
% Smooth an image with a Gaussian kernel
% FORMAT smooth_img_in_mem(img,fwhm) 
% img          - Image
% fwhm         - Full-width at half maximum
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging   
if nargin<2, fwhm      = 10; end
if nargin<3, VoxelSize = 1; end

if numel(fwhm) == 1
    fwhm = fwhm*ones(1,3);
end
if numel(VoxelSize) == 1
    VoxelSize = VoxelSize*ones(1,3);
end

simg = zeros(size(img));

fwhm = fwhm./VoxelSize;                        % voxel anisotropy
s1   = fwhm/sqrt(8*log(2));              % FWHM -> Gaussian parameter

x  = round(6*s1(1)); x = -x:x; x = spm_smoothkern(fwhm(1),x,1); x  = x/sum(x);
y  = round(6*s1(2)); y = -y:y; y = spm_smoothkern(fwhm(2),y,1); y  = y/sum(y);
z  = round(6*s1(3)); z = -z:z; z = spm_smoothkern(fwhm(3),z,1); z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;

spm_conv_vol(img,simg,x,y,z,-[i,j,k]);
%==========================================================================

%==========================================================================
function [V,W,C,BW,El] = spm_hist(X,varargin)
% _________________________________________________________________________
%
% Compute the (joint) histogram of a (multidimensional) dataset
%
% FORMAT [V,W,C] = spm_misc('hist',X,B..)
% FORMAT [V,W]   = spm_misc('hist',X,C..)
%
% MANDATORY
% ---------
% X - NxP matrix of observed values
% 
% OPTIONAL
% --------
% B - 1x1 or 1xP number of bins [64]
%   or
% C - Bx1 ordered bin centres (or 1xP cell of bin centres)
%
% KEYWORD
% -------
% KeepZero - Keep bins with zero observations [true]
% Missing  - Keep rows with missing data [false]
%            Additional bins are created for missing values.
% Reshape  - Reshape W and V so that their lattice is B1xB2x... [false]
% Smooth   - FWHM of the smoothing kernel (in bins) [0]
% Verbose  - Verbosity level [0]
%
% OUTPUT
% ------
% V  - prod(Bp) x P matrix of multidimensional values (bin centres)
% W  - prod(Bp) x 1 vector of weights (bin counts)
% C  - 1xP cell of Bx1 bin centres
% BW - 1xP bin widths
%
% (B can be smaller that the specified number of bins if KeepZero = false)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    
% -------------------------------------------------------------------------
% Parse inputs
p = inputParser;
p.FunctionName = 'spm_imbasics(''hist'')';
p.addRequired('X',                  @isnumeric);
p.addOptional('B',           64,    @(X) isnumeric(X) || iscell(X));
p.addParameter('KeepZero',   true,  @isscalar);
p.addParameter('Missing',    false, @isscalar);
p.addParameter('Reshape',    false, @isscalar);
p.addParameter('Smooth',     0,     @isnumeric);
p.addParameter('Verbose',    0,     @isscalar);
p.addParameter('Labels',     {},    @iscell);
p.parse(X, varargin{:});
B      = p.Results.B;
Labels = p.Results.Labels;

% -------------------------------------------------------------------------
% Discard missing values
if ~p.Results.Missing
    missing = any(isnan(X),2);
    X       = X(~missing,:);
end

% -------------------------------------------------------------------------
% Compute bin centres / edges
P      = size(X,2); % Number of channels
N      = size(X,1); % Number of observations
minval = min(X, [], 1, 'omitnan'); % Min value / channel
maxval = max(X, [], 1, 'omitnan'); % Max value / channel
if ~iscell(B) && size(B,1) == 1
% Number of bins provided
    E = B;
    if numel(B) < P
        E = padarray(E, [0, P-numel(B)], 'replicate', 'post');
    end
    BW = (maxval - minval)./B;
    E  = num2cell(E);
else
% Bin centres provided
    if ~iscell(B)
        if size(B,2) == 1
            B = repmat(B(:),1,P);
        end
        B = num2cell(B, 1);
    end
        
    E  = cell(1,P);
    BW = cell(1,P);
    for c=1:P
        E{c}  = (B{c}(2:end) + B{c}(1:end-1))/2;
        E{c}  = [minval(c); E{c}; maxval(c)]';
        BW{c} = E{c}(2:end) - E{c}(1:end-1);
    end
end
clear B

% -------------------------------------------------------------------------
% Discretize data
I  = cell(1,P);
V  = cell(1,P);
dim = zeros(1,P);
hasnan = zeros(1,P,'logical');
for c=1:P
    [I{c},V{c}]       = discretize(X(:,c),E{c});    
    I{c}              = single(I{c});
    I{c}(isnan(I{c})) = numel(V{c});
    V{c}              = (V{c}(2:end) + V{c}(1:end-1))/2;
    hasnan(c)         = any(isnan(X(:,c)));
    dim(c)            = numel(V{c}) + hasnan(c);
    if hasnan(c)
        V{c}(end+1)   = NaN;
    end
end
clear E X

% -------------------------------------------------------------------------
% Count
if numel(dim) == 1
    linI = [I{:}];
else
    linI = sub2ind(dim, I{:});
end
clear I

% if ~isempty(Labels)
%     nlabels = size(Labels{2},1);
%     El      = zeros([prod(dim) nlabels],'single');
%     
%     if 0
%         [a,b] = unique(linI);
%         for i=1:numel(a)
%             msk_rhs  = linI==a(i);
%             msk_lhs  = find(msk_rhs);        
%             labels_i = Labels{1}(msk_rhs);
%             for l=1:nlabels
%                 El(msk_lhs,l) = sum(labels_i==(l - 1));
%             end        
%         end
%     else    
%         for i=1:prod(dim)
%             for l=1:nlabels                
%                 El(i,l) = sum(Labels{1}(linI==i)==(l - 1));
%             end  
%             El(i,:)
%         end
%     end
% 
%     El = bsxfun(@rdivide,El,sum(El,2));    
% else
%     El = [];
% end

W = histcounts(linI, 1:prod(dim)+1); clear linI
C = V;
V = combvec(V{:});
V = V.';
W = W.';

if p.Results.Reshape && ~p.Results.KeepZero
    error('spm_imbasics::hist - Cannot Reshape and not KeepZero')
end


% -------------------------------------------------------------------------
% Smooth
if p.Results.Smooth
    W = reshape(W, dim);
    lim = ceil(4/2.355*p.Results.Smooth);
    ker = spm_smoothkern(p.Results.Smooth, -lim:lim, 0);
    ker = ker(ker~=0);
    for c=1:P
        if hasnan(c)
            W1        = W;
            subs      = cell(1,P);
            [subs{:}] = deal(':');
            subs{c}   = 1:size(W,c)-1;
            W = subsref(W1, struct('type', '()', 'subs', {subs}));
        end
        W = convn(W, reshape(ker, [ones(1,c-1) numel(ker) 1]), 'same');
        if hasnan(c)
            [W1,W] = deal(W,W1);
            W = subsasgn(W, struct('type', '()', 'subs', {subs}), W1);
            clear W1
        end
        
    end
    W = W(:);
end

% -------------------------------------------------------------------------
% Reshape
if p.Results.Reshape
    W = reshape(W, dim);
    V = reshape(V, [dim P]);
end

% -------------------------------------------------------------------------
% Remove empty bins
if ~p.Results.KeepZero
    empty = W == 0;
    W     = W(~empty);
    V     = V(~empty,:);
end
%==========================================================================

%==========================================================================
function nfname1 = create_2d_slice(fname,deg,axis_2d,clean_up)
% Extract the central 2D slice from a 3D volume.
% FORMAT nfname1 = create_2d_slice(fname,axis_2d,clean_up)
% fname   - Input filename
% axis_2d - Axis to extract along [axis_2d=3]
% nfname1 - Filename of 2d image
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging   

if nargin<2, deg      = 0;    end
if nargin<3, axis_2d  = 3;    end
if nargin<4, clean_up = true; end

% Create bounding box
V  = spm_vol(fname);
dm = V.dim;
if axis_2d     == 1
    d1 = floor(dm(1)/2) + 1;
    bb = [d1 d1;-inf inf;-inf inf];   
elseif axis_2d == 2
    d1 = floor(dm(2)/2) + 1;
    bb = [-inf inf;d1 d1;-inf inf];
elseif axis_2d == 3 
    d1 = floor(dm(3)/2) + 1;
    bb = [-inf inf;-inf inf;d1 d1];
end                

% Crop according to bounding-box
spm_impreproc('subvol',V,bb','2d_',deg);      
[pth,nam,ext] = fileparts(fname);   
nfname1       = fullfile(pth,['2d_' nam ext]);   
if clean_up, delete(fname); end

if axis_2d == 1 || axis_2d == 2
    % Make sure 1D plane is in z dimension
    Nii  = nifti(nfname1);
    mat  = Nii.mat;
    
    % Permute image data and apply permutation matrix to orientation matrix
    if axis_2d == 1
        img = permute(Nii.dat(:,:,:),[2 3 1]);            
        P   = [0 1 0 0; 0 0 1 0; 1 0 0 0; 0 0 0 1];
    else
        img = permute(Nii.dat(:,:,:),[1 3 2]);        
        P   = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
    end   
    mat     = P*mat*P';
    dm      = [size(img) 1];
    
    % Overwrite image data
    VO             = spm_vol(nfname1);
    VO.dim(1:3)    = dm(1:3);        
    VO.mat         = mat;
    VO             = spm_create_vol(VO);        
    Nii            = nifti(VO.fname);    
    Nii.dat(:,:,:) = img; 
end
%==========================================================================

%==========================================================================
function [BB,vx] = compute_bb(img,mat,dm,thr,premul)
% Compute volume's bounding box, for full field of view or object bounds
% Modified version of Ged's spm_get_bbox.
% FORMAT [BB,vx] = compute_bb(img,mat,dm,thr,premul)
% img - image volume
% mat - orientation matrix
% dm  - image dimensions
% thr - threshold, such that BB contains voxels with intensities > thr
%       or strings 'nz', 'nn', fv', for non-zero, non-NaN, or field of view
%       where 'fv' (the default) uses only the image's header information.
%
% BB  - a [2 x 3] array of the min and max X, Y, and Z coordinates {mm},
%       i.e. BB = [minX minY minZ; maxX maxY maxZ].
% vx  - a [1 x 3] vector of voxel dimensions {mm}.
%__________________________________________________________________________
% Copyright (C) 2011-2013 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: spm_get_bbox.m 5398 2013-04-12 12:37:00Z ged $

% Undocumented expert options:
% V           - can be a 4D @nifti object (but not 5D), image-based BBs
%               will be computed using "all" along the 4th dimension.
% thr = 'old' - reproduce spm_write_sn/bbvox_from_V (and elsewhere)
% premul      - a matrix that premultiplies V.mat, as used in spm_orthviews

%-Compute voxel dimensions (for compatibility with bbvox_from_V)
%--------------------------------------------------------------------------
P  = spm_imatrix(mat);
vx = P(7:9);
% the above agrees with sqrt(sum(V.mat(1:3,1:3).^2)) for simple rotations,
% and seems more appropriate if there are reflections and/or skews.
% Note that spm_imatrix(diag([-1 1 1 1])) is [-1 1 1] as expected.

%-Compute bounding box
%--------------------------------------------------------------------------
if nargin < 2 || isempty(thr) || strcmpi(thr, 'fv')
    % overall field-of-view bounding box from header information
    corners = [
        1    1    1    1
        1    1    dm(3) 1
        1    dm(2) 1    1
        1    dm(2) dm(3) 1
        dm(1) 1    1    1
        dm(1) 1    dm(3) 1
        dm(1) dm(2) 1    1
        dm(1) dm(2) dm(3) 1
        ]';
    XYZ = mat(1:3, :) * corners;
elseif strcmpi(thr, 'old')
    % code from spm_write_sn/bbvox_from_V (and other places)
    % NB: main difference is that vx(1)<0 gives descending BB(:,1),
    % shouldn't be used if V.mat contains rotations or skews.
    o  = mat\[0 0 0 1]';
    o  = o(1:3)';
    BB = [-vx.*(o-1) ; vx.*(dm(1:3)-o)];
    if exist('premul', 'var')
        warning('spm_get_bbox:old_and_premul', 'old method ignores premul')
    end
else
    % image-based bounding box using voxel intensities
    if ischar(thr)
        switch lower(thr)
            case 'nn'  % non-NaN, though include +/- Inf in computation
                img = ~isnan(img);
            case 'nz'  % special case of non-zero (rather than > 0)
                img = ~isnan(img) & img ~= 0;
            otherwise
                error('Unknown threshold type %s', thr)
        end
    else
        % treat thr as numeric threshold
        img = img > thr;
    end
    if ndims(img) == 4
        img = all(img, 4);
    end
    if nnz(img) == 0
        warning('spm_get_bbox:nothing', ...
            'Threshold leaves no voxels, returning full field of view');
        if exist('premul', 'var')
            [BB,vx] = compute_bb(img, mat, dm, 'fv', premul);
        else
            [BB,vx] = compute_bb(img, mat, dm, 'fv');
        end
        return
    else
        img     = find(img); % (clears img to save memory)
        [X Y Z] = ind2sub(dm, img);
        XYZ     = mat(1:3, :) * [X Y Z ones(size(X))]';
    end
end

if ~exist('BB', 'var') % exists already if 'old' case chosen above
    if exist('premul', 'var')
        XYZ = premul(1:3, :) * [XYZ; ones(1, size(XYZ, 2))];
    end
    BB = [
        min(XYZ, [], 2)'
        max(XYZ, [], 2)'
        ];
end
%==========================================================================

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

%==========================================================================
function out = get_type(var)
tmp = whos('var');
out = tmp.class;
%==========================================================================