function varargout = spm_misc(varargin)
%__________________________________________________________________________
% Collection of miscellaneous functions.
%
% FORMAT [M_avg,d] = spm_misc('compute_avg_mat',Mat0,dims)
% FORMAT Nii       = spm_misc('create_nii',pth,dat,mat,dtype,descrip,scl)
% FORMAT y         = spm_misc('linspace_vec',x1,x2,n)
% FORMAT spm_misc('manage_parpool',num_workers)
% FORMAT nw        = spm_misc('nbr_parfor_workers')
% FORMAT vx        = spm_misc('vxsize',M)
% FORMAT msk       = spm_misc('msk_modality',f,modality)
% FORMAT gain      = spm_misc('get_gain',vals)
% FORMAT [B, rind] = spm_misc('affine_basis', type, flat)
% FORMAT             spm_misc('clean_holly_mbrud')
%
% FORMAT help spm_parfor>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging
if nargin == 0
    help spm_parfor
    error('Not enough argument. Type ''help spm_parfor'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'compute_avg_mat'
        [varargout{1:nargout}] = compute_avg_mat(varargin{:});  
    case 'create_nii'
        [varargout{1:nargout}] = create_nii(varargin{:});    
    case 'linspace_vec'
        [varargout{1:nargout}] = linspace_vec(varargin{:});        
    case 'manage_parpool'
        [varargout{1:nargout}] = manage_parpool(varargin{:});
    case 'nbr_parfor_workers'
        [varargout{1:nargout}] = nbr_parfor_workers(varargin{:});            
    case 'vxsize'
        [varargout{1:nargout}] = vxsize(varargin{:});              
    case 'msk_modality'
        [varargout{1:nargout}] = msk_modality(varargin{:});                
    case 'get_gain'
        [varargout{1:nargout}] = get_gain(varargin{:});          
    case 'affine_basis'
        [varargout{1:nargout}] = affine_basis(varargin{:});            
    case 'clean_holly_mbrud'
        [varargout{1:nargout}] = clean_holly_mbrud(varargin{:});              
    otherwise
        help spm_parfor
        error('Unknown function %s. Type ''help spm_parfor'' for help.', id)
end
%==========================================================================

%==========================================================================
function Nii = create_nii(pth,dat,mat,dtype,descrip,offset,scl_slope,scl_inter)
% Create a NIfTI file
% FORMAT Nii = create_nii(pth,dat,mat,dtype,descrip,scl)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<6, offset    = 0; end
if nargin<7, scl_slope = 1; end
if nargin<8, scl_inter = 0; end

if exist(pth,'file')==2, delete(pth); end

Nii         = nifti;
dm          = size(dat);
Nii.dat     = file_array(pth,dm,dtype,offset,scl_slope,scl_inter);
Nii.mat     = mat;
Nii.mat0    = mat;
Nii.descrip = descrip;
create(Nii);

if numel(dm)==4
    Nii.dat(:,:,:,:) = dat;
else
    Nii.dat(:,:,:)   = dat;
end
%==========================================================================

%==========================================================================
function y = linspace_vec(x1,x2,n)
% Generalisation of MATLAB's linspace to vectors
% FORMAT y = linspace_vec(x1,x2,n)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if numel(x1)==1 && numel(x2)==1
    y = linspace(x1,x2,n);
    return
end

x1 = squeeze(x1); x2 = squeeze(x2);

if ndims(x1)~= ndims(x2) || any(size(x1)~= size(x2))
    error('d1 and d2 must have the same number of dimension and the same size'),
end

NDim = ndims(x1);
if NDim==2 && any(size(x1)==1)
    NDim = NDim-1;
    if all(size(x1)==1)
        NDim = 0;
    end
end

pp      = (0:n-2)./(floor(n)-1);

Sum1 = kron(x1, ones(1,n-1));
Sum2 = kron((x2-x1), pp);
y    = cat(NDim+1, Sum1  + Sum2, shiftdim(x2, size(x1, 1)==1 ));
%==========================================================================

%==========================================================================
function manage_parpool(num_workers)
% Start/stop parallel pool
% FORMAT manage_parpool(num_workers)
% num_workers - Number of parfor workers
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
nw = spm_misc('nbr_parfor_workers');
if num_workers>nw
    num_workers = nw;
end

poolobj = gcp('nocreate');
if ~isempty(poolobj) && num_workers==0
    delete(poolobj);
elseif ~isempty(poolobj) && poolobj.NumWorkers~=num_workers
    delete(poolobj);
    parpool('local',num_workers);
elseif isempty(poolobj) && num_workers
    parpool('local',num_workers);
end
%==========================================================================

%==========================================================================
function nw = nbr_parfor_workers
% Get number of CPU cores
% FORMAT nw = nbr_parfor_workers
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
c  = parcluster('local');
nw = c.NumWorkers;
%==========================================================================

%==========================================================================
function vx = vxsize(M)
% Get voxel size
% FORMAT vx = vxsize(M)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
M  = M(1:3,1:3);
vx = sqrt(sum(M.^2));
%==========================================================================

%==========================================================================
function [M_avg,d] = compute_avg_mat(Mat0,dims)
% Compute an average voxel-to-world mapping and suitable dimensions
% FORMAT [M_avg,d] = compute_avg_mat(Mat0,dims)
% Mat0  - array of matrices (4x4xN)
% dims  - image dimensions (Nx3)
% M_avg - voxel-to-world mapping
% d     - dimensions for average image
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Rigid-body matrices computed from exp(p(1)*B(:,:,1)+p(2)+B(:,:,2)...)
%--------------------------------------------------------------------------
B = se3_basis;

% Find combination of 90 degree rotations and flips that brings all
% the matrices closest to axial
%--------------------------------------------------------------------------
Matrices = Mat0;
pmatrix  = [1,2,3; 2,1,3; 3,1,2; 3,2,1; 1,3,2; 2,3,1];
for i=1:size(Matrices,3)
    vx    = sqrt(sum(Matrices(1:3,1:3,i).^2));
    tmp   = Matrices(:,:,i)/diag([vx 1]);
    R     = tmp(1:3,1:3);
    minss = Inf;
    minR  = eye(3);
    for i1=1:6,
        R1 = zeros(3);
        R1(pmatrix(i1,1),1)=1;
        R1(pmatrix(i1,2),2)=1;
        R1(pmatrix(i1,3),3)=1;
        for i2=0:7,
            F  = diag([bitand(i2,1)*2-1, bitand(i2,2)-1, bitand(i2,4)/2-1]);            
            R2 = F*R1;
            ss = sum(sum((R/R2-eye(3)).^2));
            if ss<minss,
                minss = ss;
                minR  = R2;
            end
        end
    end
    rdim = abs(minR*dims(i,:)');
    R2   = inv(minR);
    minR = [R2 R2*((sum(R2,1)'-1)/2.*(rdim+1)); 0 0 0 1];
    Matrices(:,:,i) = Matrices(:,:,i)*minR;
end

% Average of these matrices
%--------------------------------------------------------------------------
M_avg = spm_meanm(Matrices);

% If average involves shears, then find the closest matrix that does not
% require them
%--------------------------------------------------------------------------
p = spm_imatrix(M_avg);
if sum(p(10:12).^2)>1e-8,

    % Zooms computed from exp(p(7)*B2(:,:,1)+p(8)*B2(:,:,2)+p(9)*B2(:,:,3))
    %-----------------------------------------------------------------------
    B2        = zeros(4,4,3);
    B2(1,1,1) = 1;
    B2(2,2,2) = 1;
    B2(3,3,3) = 1;

    p      = zeros(9,1); % Parameters
    for it=1:10000,
        [R,dR] = spm_dexpm(p(1:6),B);  % Rotations + Translations
        [Z,dZ] = spm_dexpm(p(7:9),B2); % Zooms

        M  = R*Z; % Voxel-to-world estimate
        dM = zeros(4,4,6);
        for i=1:6, dM(:,:,i)   = dR(:,:,i)*Z; end
        for i=1:3, dM(:,:,i+6) = R*dZ(:,:,i); end
        dM = reshape(dM,[16,9]);

        d   = M(:)-M_avg(:); % Difference
        gr  = dM'*d;         % Gradient
        Hes = dM'*dM;        % Hessian
        p   = p - Hes\gr;    % Gauss-Newton update
        if sum(gr.^2)<1e-8, break; end
    end
    M_avg = M;
end

% Ensure that the FoV covers all images, with a few voxels to spare
%--------------------------------------------------------------------------
mn    =  Inf*ones(3,1);
mx    = -Inf*ones(3,1);
for i=1:size(Mat0,3),
    dm      = [dims(i,:) 1 1];
    corners = [
        1 dm(1)    1  dm(1)   1  dm(1)    1  dm(1)
        1    1  dm(2) dm(2)   1     1  dm(2) dm(2)
        1    1     1     1 dm(3) dm(3) dm(3) dm(3)
        1    1     1     1    1     1     1     1];
    M  = M_avg\Mat0(:,:,i);
    vx = M(1:3,:)*corners;
    mx = max(mx,max(vx,[],2));
    mn = min(mn,min(vx,[],2));
end
mx    = ceil(mx-(mx-mn)*0.05);
mn    = floor(mn+([mx(1:2)-mn(1:2);0])*0.05);
d     = (mx-mn+1)';
% d     = (mx-mn)';
if ~any(dims(:,3)>1)
   d(3)  = 1;
   mn(3) = 1;
end
M_avg = M_avg * [eye(3) mn-1; 0 0 0 1];
M_avg(4,:)=[0 0 0 1];
%==========================================================================

%==========================================================================
function msk = msk_modality(f,modality,mskonlynan)
% Get a mask for masking voxels in different imaging modalities
% FORMAT msk = msk_modality(f,modality)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin < 3, mskonlynan = false; end

if mskonlynan
    if strcmpi(modality,'mri')
        msk = isfinite(f) & (f>=0);
    elseif strcmp(modality,'CT')
        msk = isfinite(f);
    end
else
    if strcmpi(modality,'mri')
        msk = isfinite(f) & (f>0);
    elseif strcmp(modality,'CT')
%         msk = isfinite(f) & (f~=0) & (f>-200) & (f<3000);
%         msk = msk | (f>=-1010) & (f<=-990);
        msk = isfinite(f) & (f~=0) & (f>-1020) & (f<3000);
    end
end
%==========================================================================

%==========================================================================
function gain = get_gain(vals)
% FORMAT gain = get_gain(vals)
%
% vals - A vector of values
%
% gain - Computed gain
%
% Compute gain --- usually used to determine a stopping criteria when
% optimising
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
vals = vals(:);
gain = abs((vals(end - 1) - vals(end))/(max(vals(isfinite(vals))) - min(vals(isfinite(vals)))));   
%==========================================================================

%==========================================================================
function clean_holly_mbrud
fprintf('Start cleaning Holly trash... ')
system('rm -rf /data/mbrud/.Trash-1904');
system('rm -rf /data-scratch/mbrud/.Trash-1904');
fprintf('done!\n')
%==========================================================================

%==========================================================================
function [B, rind] = affine_basis(type, flat)
% FORMAT [B, rind] = affine_basis(type, ('2d'))
% type - * 'translation'
%        * 'rotation'
%        * 'rigid'      or 6
%        * 'similitude' or 7
%        * 'affine'     or 12 [default]
% 
% B    - 4x4xQ array.
% rind - Indices of basis that shoudl be reularised (all but tr/rot)
%
% Returns a Lie algebra basis system encoding for one of the above
% transformation types.

if nargin < 2
    flat = '3d';
    if nargin < 1
        type = 'affine';
    end
end

flat = ischar(flat) && strcmpi(flat, '2d');
if ~ischar(type)
    type = num2str(type);
end

type = deblank(lower(type));

% --- Define basis vectors

Bt = [ 0 0 0 1   0 0 0 0   0 0 0 0 ;
       0 0 0 0   0 0 0 1   0 0 0 0 ;
       0 0 0 0   0 0 0 0   0 0 0 1 ;
       0 0 0 0   0 0 0 0   0 0 0 0 ];

Br = [ 0 -1 0 0   0 0 -1 0   0 0  0 0 ;
       1  0 0 0   0 0  0 0   0 0 -1 0 ;
       0  0 0 0   1 0  0 0   0 1  0 0 ;
       0  0 0 0   0 0  0 0   0 0  0 0 ];

Bsim = [ 1 0 0 0 ;
         0 1 0 0 ;
         0 0 1 0 ;
         0 0 0 0 ];

Bscl = [ 1 0 0 0   0 0 0 0   0 0 0 0 ;
         0 0 0 0   0 1 0 0   0 0 0 0 ;
         0 0 0 0   0 0 0 0   0 0 1 0 ;
         0 0 0 0   0 0 0 0   0 0 0 0 ];

Bshr = [ 0 1 0 0   0 0 1 0   0 0 0 0 ;
         1 0 0 0   0 0 0 0   0 0 1 0 ;
         0 0 0 0   1 0 0 0   0 1 0 0 ;
         0 0 0 0   0 0 0 0   0 0 0 0 ];

% --- Remove 3D basis if 2D
if flat
    Bt = Bt(:,1:8);
    Br = Br(:,1:4);
    Bsim = [ 1 0 0 0 ;
             0 1 0 0 ;
             0 0 0 0 ;
             0 0 0 0 ];
    Bscl = Bscl(:,1:8);
    Bshr = Bshr(:,1:4);
end

% --- Build complete basis
switch type
    case 'translation'
        B = Bt;
        rind = [];
    case 'rotation'
        B = Br;
        rind = [];
    case {'rigid', '6'}
        B = [Bt Br];
        rind = [];
    case {'similitude', '7'}
        B = [Bt Br Bsim];
        if flat,  rind = [4 5];
        else      rind = [7 8 9];  end
    case {'affine', '12'}
        B = [Bt Br Bscl Bshr];
        if flat,  rind = [4 5 6];
        else      rind = [7 8 9 10 11 12];  end
end
B = reshape(B, 4, 4, []);
%==========================================================================

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

%==========================================================================
function B = se3_basis
% Basis functions for the lie algebra of the special Eucliden group
% (SE(3)).
B        = zeros(4,4,6);
B(1,4,1) = 1;
B(2,4,2) = 1;
B(3,4,3) = 1;
B([1,2],[1,2],4) = [0 1;-1 0];
B([3,1],[3,1],5) = [0 1;-1 0];
B([2,3],[2,3],6) = [0 1;-1 0];
%==========================================================================