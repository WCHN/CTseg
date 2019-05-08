function varargout = spm_impreproc(varargin)
%__________________________________________________________________________
% Collection of tools for image pre-processing.
%
% FORMAT [Affine,bb] = spm_impreproc('atlas_crop',P,Affine,prefix,rem_neck)
% FORMAT R      = spm_impreproc('rigid_align',P)
% FORMAT V      = spm_impreproc('reg_and_reslice',V)
% FORMAT nfname = spm_impreproc('downsample_inplane',fname)
% FORMAT nfname = spm_impreproc('downsample_throughplane',fname)
% FORMAT nii    = spm_impreproc('mult_bb_crop',nii,bb,verbose)
% FORMAT V      = spm_impreproc('resize_ims',V,V_ref,vx,prefix,deg)
% FORMAT spm_impreproc('subvol',V,bb,prefix)
% FORMAT spm_impreproc('nm_reorient',Vin,vx,deg)
% FORMAT spm_impreproc('reset_origin',P)
% FORMAT pth = change_vx_size(pth,vx,deg,prefix)
%
% FORMAT help spm_impreproc>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin == 0
    help spm_impreproc
    error('Not enough argument. Type ''help spm_impreproc'' for help.');
end
id       = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'atlas_crop'
        [varargout{1:nargout}] = atlas_crop(varargin{:});  
    case 'nm_reorient'
        [varargout{1:nargout}] = nm_reorient(varargin{:});                   
    case 'coreg'
        [varargout{1:nargout}] = coreg(varargin{:});                 
    case 'reset_origin'
        [varargout{1:nargout}] = reset_origin(varargin{:});              
    case 'reslice'
        [varargout{1:nargout}] = reslice(varargin{:});                     
    case 'rigid_align'
        [varargout{1:nargout}] = rigid_align(varargin{:});                  
    case 'subvol'
        [varargout{1:nargout}] = subvol(varargin{:});    
    case 'downsample_inplane'
        [varargout{1:nargout}] = downsample_inplane(varargin{:});      
    case 'downsample_throughplane'
        [varargout{1:nargout}] = downsample_throughplane(varargin{:});              
    case 'mult_bb_crop'
        [varargout{1:nargout}] = mult_bb_crop(varargin{:});    
    case 'resize_ims'
        [varargout{1:nargout}] = resize_ims(varargin{:});            
    case 'change_vx_size'
        [varargout{1:nargout}] = change_vx_size(varargin{:});               
    otherwise
        help spm_impreproc
        error('Unknown function %s. Type ''help spm_impreproc'' for help.', id)
end
%==========================================================================

%==========================================================================
function [Affine,bb] = atlas_crop(P,prefix,rem_neck)
% Removes air outside of head
% FORMAT [Affine,bb] = atlas_crop(P,Affine,prefix,rem_neck)
% P        - Path to NIfTI file
% prefix   - File prefix (if empty -> overwrites) ['']
% rem_neck - Remove neck/spine [false]
% bb - Computed bounding box
%
% This function rigidly registers the SPM atlas to an image and then
% removes image data outside of the head.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<2, prefix   = ''; end
if nargin<3, rem_neck = 0;  end

% Locate TPM.nii in SPM
pth_tpm = fullfile(spm('dir'),'tpm','TPM.nii,');

Vin  = spm_vol(P);
Vtpm = spm_vol(pth_tpm);

mat    = Vin.mat;
mattpm = Vtpm.mat;

tpm = spm_load_priors8(Vtpm);    

V = spm_vol(P);

M               = V(1).mat;
c               = (V(1).dim+1)/2;
V(1).mat(1:3,4) = -M(1:3,1:3)*c(:);
[Affine1,ll1]   = spm_maff8(V(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid
Affine1         = Affine1*(V(1).mat/M);

% Run using the origin from the header
V(1).mat      = M;
[Affine2,ll2] = spm_maff8(V(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid

% Pick the result with the best fit
if ll1>ll2, Affine  = Affine1; else Affine  = Affine2; end

Affine = spm_maff8(P,2,32,tpm,Affine,'mni');
Affine = spm_maff8(P,2,1,tpm,Affine,'mni');

% Voxel locations in TPM.nii
Ltpm1 = [120 72.2 37.3 1]'; Ltpm2 = [120 72.2 75.9 1]';
Rtpm1 = [3  72.2 37.3 1]'; Rtpm2 = [3  72.2 75.9 1]';

Stpm1 = [58.6 42.6 119 1]'; Stpm2 = [58.60 99.4 119 1]';
if rem_neck==2
    Itpm1 = [61 42 31   1]'; Itpm2 = [61 107 29 1]';   
elseif rem_neck==1
    Itpm1 = [58.6 39.4 2.5   1]'; Itpm2 = [58.60 99.4 2.5  1]';    
else
    Itpm1 = [58.6 39.4 -200  1]'; Itpm2 = [58.60 99.4 -200 1]';
end
Atpm1 = [58.6 144 28.4 1]'; Atpm2 = [58.60 144 82.3 1]';
Ptpm1 = [58.6 3.5  28.4 1]'; Ptpm2 = [58.60 3.5 82.3 1]'; 

% Voxel locations in input
T  = mat\(Affine\mattpm);
L1 = T*Ltpm1; L2 = T*Ltpm2;
R1 = T*Rtpm1; R2 = T*Rtpm2;
U1 = T*Stpm1; U2 = T*Stpm2;
D1 = T*Itpm1; D2 = T*Itpm2;
A1 = T*Atpm1; A2 = T*Atpm2;
P1 = T*Ptpm1; P2 = T*Ptpm2;

% Bounding-box
bb = zeros(2,3);
for i=1:3
    X       = [L1(i) R1(i) U1(i) D1(i) A1(i) P1(i)...
               L2(i) R2(i) U2(i) D2(i) A2(i) P2(i)];
    bb(1,i) = max(X);
    bb(2,i) = min(X);
end

% Do cropping
spm_impreproc('subvol',Vin,bb,prefix);      
%==========================================================================

%==========================================================================
function VO = nm_reorient(Vin,vx,prefix,deg)
% Re-orient images
% FORMAT nm_reorient(Vin,vx,type,deg)
% Vin  - SPM volume objects
% vx   - New voxel size
% type - Order of interpolation
% prefix - Prefix of file to write
% def    - Degree of interpolation [0]
%
% The function reslices the input images to a resolution of vx mm.
% Output images (with the prefix "pn_r") are written in the transverse
% orientation (using information from the ".mat" files).
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
spm('defaults', 'FMRI');

if length(vx)<3
    vx=[vx vx vx];
end
if nargin<3, prefix = 'ro_'; end
if nargin<4, deg    = 4; end

% If no arguments, then prompt for images
%PP = spm_get([1 Inf],'*.img','Select files to reorient');

% Get information about the image volumes
VV = spm_vol(Vin);

for V=VV', % Loop over images

    % The corners of the current volume
    d = V.dim(1:3);
    c = [	1    1    1    1
        1    1    d(3) 1
        1    d(2) 1    1
        1    d(2) d(3) 1
        d(1) 1    1    1
        d(1) 1    d(3) 1
        d(1) d(2) 1    1
        d(1) d(2) d(3) 1]';

    % The corners of the volume in mm space
    tc = V.mat(1:3,1:4)*c;
    if spm_flip_analyze_images, tc(1,:) = -tc(1,:); end;

    % Max and min co-ordinates for determining a bounding-box
    mx = round(max(tc,[],2)');
    mn = round(min(tc,[],2)');

    % Translate so that minimum moves to [1,1,1]
    % This is the key bit for changing voxel sizes,
    % output orientations etc.
    mat = spm_matrix(mn)*diag([vx 1])*spm_matrix(-[1 1 1]);

    % Dimensions in mm
    dim = ceil((mat\[mx 1]')');

    % Output image based on information from the original
    VO               = V;

    % Create a filename for the output image (prefixed by 'r')
    [lpath,name,ext] = fileparts(V.fname);
    VO.fname         = fullfile(lpath,[prefix name ext]);

    % Dimensions of output image
    VO.dim(1:3)      = dim(1:3);

    % Voxel-to-world transform of output image
    if spm_flip_analyze_images, mat = diag([-1 1 1 1])*mat; end;
    VO.mat           = mat;

    % Initialise plot of how far reslicing has gone
    %spm_progress_bar('Init',dim(3),'reslicing...','planes completed');

    % Create .hdr and open output .img
    VO = spm_create_vol(VO);

    for i=1:dim(3), % Loop over slices of output image

        % Mapping from slice i of the output image,
        % to voxels of the input image
        M   = inv(spm_matrix([0 0 -i])*inv(VO.mat)*V.mat);

        % Extract this slice according to the mapping
        img = spm_slice_vol(V,M,dim(1:2),deg);

        % Write this slice to output image
        spm_write_plane(VO,img,i);

        % Update the progress bar
        %spm_progress_bar('Set',i);

    end; % End loop over output slices

    % Get rid of the progress bar
    %spm_progress_bar('Clear');

end; % End loop over images
%==========================================================================

%==========================================================================
function reset_origin(P,orig)
% Reset origin of image
% FORMAT reset_origin(P)
% P = Path to NIfTI image
%
% OBS: Image will have the matrix in its header adjusted.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<2, orig = []; end

V   = spm_vol(P);
M   = V.mat;
dim = V.dim;
vx  = sqrt(sum(M(1:3,1:3).^2));

if det(M(1:3,1:3))<0
    vx(1) = -vx(1); 
end

if isempty(orig)
    orig = (dim(1:3)+1)/2;
end

off  = -vx.*orig;
M1   = [vx(1) 0      0         off(1)
           0      vx(2) 0      off(2)
           0      0      vx(3) off(3)
           0      0      0      1];

spm_get_space(P,M1);   
%==========================================================================

%==========================================================================
function R = rigid_align(P)
% Reposition an image by affine aligning to MNI space and Procrustes adjustment
% FORMAT rigid_align(P)
% P - name of NIfTI image
% R - Affine matrix
%
% OBS: Image will have the matrix in its header adjusted.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Load tissue probability data
tpm = fullfile(spm('dir'),'tpm','TPM.nii,');
tpm = [repmat(tpm,[6 1]) num2str((1:6)')];
tpm = spm_load_priors8(tpm);

% Do the affine registration
V = spm_vol(P);

M               = V(1).mat;
c               = (V(1).dim+1)/2;
V(1).mat(1:3,4) = -M(1:3,1:3)*c(:);
[Affine1,ll1]   = spm_maff8(V(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid
Affine1         = Affine1*(V(1).mat/M);

% Run using the origin from the header
V(1).mat      = M;
[Affine2,ll2] = spm_maff8(V(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid

% Pick the result with the best fit
if ll1>ll2, Affine  = Affine1; else Affine  = Affine2; end

Affine = spm_maff8(P,2,32,tpm,Affine,'mni'); % Heavily regularised
Affine = spm_maff8(P,2,1 ,tpm,Affine,'mni'); % Lightly regularised

% Load header
Nii    = nifti(P);

% Generate mm coordinates of where deformations map from
x      = affind(rgrid(size(tpm.dat{1})),tpm.M);

% Generate mm coordinates of where deformation maps to
y1     = affind(x,inv(Affine));

% Weight the transform via GM+WM
weight = single(exp(tpm.dat{1})+exp(tpm.dat{2}));

% Weighted Procrustes analysis
[~,R]  = spm_get_closest_affine(x,y1,weight);

% Invert
% R      = inv(R);

% Write the new matrix to the header
Nii.mat = R\Nii.mat;
create(Nii);
%==========================================================================

%==========================================================================
function [V,res,source_ix] = coreg(V,ref_ix)
% Co-register images
% FORMAT V = coreg(V)
% V - SPM volume object that can contain N different modalities (e.g. T1- 
% and T2-weighted MRIs.
%
% WARNING: This function overwrites orientation matrices!
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
N = numel(V);
if N==1
    return;
end

if nargin < 2
    % Get image with smallest voxel size and pick this image as reference
    prod_vx = zeros(1,N);
    for n=1:N
        vx         = spm_misc('vxsize',V(n).mat);
        prod_vx(n) = prod(vx);
    end
    [~,ref_ix] = min(prod_vx);
end

% Set options
matlabbatch{1}.spm.spatial.coreg.estimate.ref               = {V(ref_ix).fname};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];

% Co-register
ixs       = 1:N;
source_ix = ixs(ixs~=ref_ix);
res       = cell(1,numel(source_ix));
cnt       = 1;
for n=source_ix
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {V(n).fname};         

    res{cnt} = spm_jobman('run',matlabbatch);
    cnt      = cnt + 1;
end
%==========================================================================

%==========================================================================
function [V,ref_ix] = reslice(V,deg,ref_ix)
% Re-slice images
% FORMAT V = reslice(V)
% V - SPM volume object that can contain N different modalities (e.g. T1- 
% and T2-weighted MRIs.
% ref_ix - index of reference image in V
%
% Takes medical images of the same subject and re-slices the images to the
% same dimensions. If no reference index is given, the image with the largest field of view is chosen as
% reference for the re-slicing. First order interpolation is used not to
% introduce any negative values.
%
% WARNING: This function overwrites the input data!
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Interpolation degree
if nargin<2, deg = 4; end

C = numel(V); % Number of channels
if C == 1
    % Only one channel -> No need to reslice
    return;
end

% Store all dimensions and orientation matrices
dm  = zeros(C,3);
mat = zeros(4,4,C);
for n=1:C
    dm(n,:)    = V(n).dim;
    mat(:,:,n) = V(n).mat;
end

if nargin<3
    % Get image with largest volume (for reslicing using this image as
    % reference)        
    vol = zeros(C,3);
    for n=1:C
        vx       = spm_misc('vxsize',V(n).mat);
        vol(n,:) = vx.*V(n).dim;
    end
    vol        = prod(vol,2);
    [~,ref_ix] = max(vol);
end
ixs       = 1:C;
source_ix = ixs(ixs ~= ref_ix);

% Create mask images
msk        = cell(1,C);
[x0,y0,z0] = ndgrid(1:dm(ref_ix,1),1:dm(ref_ix,2),1:dm(ref_ix,3));
for n=source_ix
    msk{n} = V(n).private.dat(:,:,:);     
    
    T  = mat(:,:,n)\mat(:,:,ref_ix);    
    x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
    y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
    z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);
    
    msk{n}                    = spm_bsplins(msk{n},x1,y1,z1,[0 0 0  0 0 0]);    
    msk{n}(~isfinite(msk{n})) = 0;
    msk{n}                    = msk{n} ~= 0;
end
    
% Use SPM batch job to reslice
matlabbatch{1}.spm.spatial.realign.write.data = {V(ref_ix).fname, V(source_ix).fname}';
matlabbatch{1}.spm.spatial.realign.write.roptions.which = [1 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.interp = deg;
matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'res_';

spm_jobman('run',matlabbatch);

for n=source_ix
    % Delete old data
    delete(V(n).fname);
    
    % Update spm_vol object
    [pth,nam,ext] = fileparts(V(n).fname);    
    V(n)          = spm_vol(fullfile(pth,['res_' nam ext]));
    
    % Mask
    img                     = single(V(n).private.dat(:,:,:));
    V(n).private.dat(:,:,:) = msk{n}.*img;     
end
%==========================================================================

%==========================================================================
function VO = subvol(V,bb,prefix,deg,constrain_mx)
% Extract a subvolume
% FORMAT VO = subvol(V,bb,prefix)
% V      - SPM volume object
% bb     - bounding box
% prefix - file prefix (if empty -> overwrites)
% VO     - resized image
%
% Example:
%     V = spm_vol(spm_select(1,'image'));
%     subvol(V,[32 64 ; 1 64 ; 1 48]');
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin < 3, prefix       = 'sv_'; end
if nargin < 4, deg          = 0;     end
if nargin < 5, constrain_mx = true;     end

bb      = round(bb);
bb      = sort(bb);
bb(1,:) = max(bb(1,:),[1 1 1]);
if constrain_mx
    bb(2,:) = min(bb(2,:),V.dim(1:3));
end

VO            = V;
[pth,nam,ext] = fileparts(V.fname);
VO.fname      = fullfile(pth,[prefix nam ext]);
VO.dim(1:3)   = diff(bb)+1;
VO.mat        = V.mat*spm_matrix((bb(1,:)-1));

VO = spm_create_vol(VO);
for z=1:VO.dim(3)
    M   = V.mat\VO.mat*spm_matrix([0 0 z]);
    img = spm_slice_vol(V,M,VO.dim(1:2),deg);
    VO  = spm_write_plane(VO,img,z);
end
%==========================================================================

%==========================================================================
function nfname = downsample_inplane(fname)
% Down-sample a NIfTI image in the high-resolution plane
% FORMAT nfname = downsample_inplane(fname)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging  
% Get image data

Nii = nifti(fname);
M0  = Nii.mat;             
X   = Nii.dat(:,:,:);  
dm0 = size(X);                                                   
vx0 = spm_misc('vxsize',M0);

% Get down-sampling factor
d = vx0(1:2);
if d(1)>=1, d(1) = 1; end
if d(2)>=1, d(2) = 1; end
d(3) = 1;

if round(d(1),3)<1 || round(d(2),3)<1        
    % NN downsampling    
    D = diag([d, 1]);          
           
    dm1 = floor(D(1:3,1:3)*dm0')';
    M1  = M0/D;       
    
    T = M0\M1;
    y = make_deformation(T,dm1);
                   
    X = spm_bsplins(X,y(:,:,:,1),y(:,:,:,2),y(:,:,:,3),[0 0 0 0 0 0]);
    clear y

    X(~isfinite(X)) = 0;
else
    M1 = M0;
end

fname         = Nii.dat.fname;
[pth,nam,ext] = fileparts(fname);
nfname        = fullfile(pth,['ds_' nam ext]);
        
spm_misc('create_nii',nfname,X,M1,Nii.dat.dtype,Nii.descrip,Nii.dat.offset,Nii.dat.scl_slope,Nii.dat.scl_inter);
%==========================================================================

%==========================================================================
function nfname = downsample_throughplane(fname)
% Down-sample a NIfTI image in the through-plane to 1 mm voxel size
% FORMAT nfname = downsample_throughplane(fname)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging  

% Get image data
Nii = nifti(fname);
M0  = Nii.mat;             
X   = Nii.dat(:,:,:);  
dm0 = size(X);                                                   
vx0 = spm_misc('vxsize',M0);

% Get down-sampling factor
d = vx0(3);

if round(d,3)<1
    % NN downsampling
    d = [1 1 d];
    D = diag([d, 1]);          
           
    dm1 = floor(D(1:3,1:3)*dm0')';
    M1  = M0/D;       
    
    T = M0\M1;
    y = make_deformation(T,dm1);
                   
    X = spm_bsplins(X,y(:,:,:,1),y(:,:,:,2),y(:,:,:,3),[0 0 0 0 0 0]);
    clear y

    X(~isfinite(X)) = 0;
else
    M1 = M0;
end

fname         = Nii.dat.fname;
[pth,nam,ext] = fileparts(fname);
nfname        = fullfile(pth,['dsz_' nam ext]);
        
spm_misc('create_nii',nfname,X,M1,Nii.dat.dtype,Nii.descrip,Nii.dat.offset,Nii.dat.scl_slope,Nii.dat.scl_inter);
%==========================================================================

%==========================================================================
function nii = mult_bb_crop(nii,BB,verbose)
% Crop image(s) according to a bunch of bounding-boxes.
% FORMAT nii = mult_bb_crop(nii,bb,verbose)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<3, verbose = true; end

fname = nii.dat.fname;

mn  = min(BB,[],3);
mx  = max(BB,[],3);
nbb = [mn(:,1) mx(:,2)];

V0 = spm_vol(fname);
od = V0(1).dim;
for k=1:numel(V0)
    spm_impreproc('subvol',V0(k),nbb','tmp');        
end

delete(fname);
[pth,nam,ext] = fileparts(V0(1).fname);
fname_tmp     = fullfile(pth,['tmp' nam ext]);
movefile(fname_tmp,fname);

V  = spm_vol(fname);
nd = V(1).dim;   

nii = nifti(fname);

if verbose
    fprintf('spm_impreproc(''mult_bb_crop'') | odm = [%d %d %d] | ndm = [%d %d %d]\n',od(1),od(2),od(3),nd(1),nd(2),nd(3));
end
%==========================================================================

%==========================================================================
function V = resize_ims(V,V_ref,vx,prefix,deg)
% Resize a bunch of images
% FORMAT V = resize_ims(V,V_ref,vx,prefix,deg)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 2, V_ref  = 'spm'; end
if nargin < 3, vx     = 'same'; end
if nargin < 4, prefix = 'res_'; end
if nargin < 5, deg    = 0;      end

% Get bounding-box
if strcmp(V_ref,'spm')
    
else
    BB = world_bb(V_ref);
end

if strcmp(vx,'same')
    vx = [];
end

% Resize image(s)
V = resize_img(V,BB,vx,prefix,deg);
%==========================================================================

%==========================================================================
function pth = change_vx_size(pth,vx,deg,prefix)
% Resize a bunch of images
% FORMAT pth = change_vx_size(pth,vx,deg,prefix)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 3, deg    = 0;     end
if nargin < 4, prefix = 'vx_'; end

if numel(vx) == 1, vx = vx*ones(1,3); end

bs = [deg deg deg  0 0 0];

% Image data
Nii  = nifti(pth);
img0 = Nii.dat(:,:,:);
mat0 = Nii.mat;
dm0  = size(img0);
vx0  = sqrt(sum(Nii.mat(1:3,1:3).^2));

%
ds   = vx0./vx;
D    = diag([ds 1]);
mat  = mat0/D;
dm   = floor(D(1:3,1:3)*dm0')';

[x0,y0,z0] = ndgrid(1:dm(1),1:dm(2),1:dm(3));

T = mat0\mat;    

x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);

coeff               = spm_bsplinc(img0,bs);
img                 = spm_bsplins(coeff,x1,y1,z1,bs);    
img(~isfinite(img)) = 0;

% Save down-sampled image
[pth,nam,ext] = fileparts(Nii.dat.fname);
nfname        = fullfile(pth,[prefix nam ext]);
Nii           = spm_misc('create_nii',nfname,img,mat,[spm_type('float32') spm_platform('bigend')],'downsampled-bs');
pth           = nfname;
%==========================================================================

%==========================================================================
function img = clean_fov(img,mat_ref,dm_ref,mat_source,dm_source)
% Zeros voxels outside of the field of view
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Get the mapping from Mref to Mmod
T = mat_source\mat_ref;

% Use ndgrid to give an array of voxel indices
[x0,y0,z0] = ndgrid(single(1:dm_ref(1)),...
                    single(1:dm_ref(2)),...
                    single(1:dm_ref(3)));

% Transform these indices to the indices that they point to in the reference image
D = cat(4,T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4), ...
          T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4), ...
          T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4));

% Mask according to whether these are < 1 or > than the dimensions of the reference image.        
msk = cell(1,3);
for i=1:3
    msk{i} = D(:,:,:,i) >= 1 & D(:,:,:,i) <= dm_source(i);
end

% Generate cleaned up image
for i=1:3
    img = msk{i}.*img;
end
%==========================================================================

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

%==========================================================================
function x = rgrid(d)
x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3)
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
%==========================================================================

%==========================================================================
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
%==========================================================================

%==========================================================================
function y = make_deformation(M,dm)
[x0,y0,z0] = ndgrid(1:dm(1),...
                    1:dm(2),...
                    1:dm(3));
y          = cat(4,M(1,1)*x0 + M(1,2)*y0 + M(1,3)*z0 + M(1,4), ...
                   M(2,1)*x0 + M(2,2)*y0 + M(2,3)*z0 + M(2,4), ...
                   M(3,1)*x0 + M(3,2)*y0 + M(3,3)*z0 + M(3,4));
%==========================================================================  

%==========================================================================
function Vo = resize_img(Vi,BB,vx,prefix,deg)

% reslice images one-by-one
Vo = spm_vol;
c  = 1;
for V=Vi'
    % (copy to allow defaulting of NaNs differently for each volume)
    if isempty(vx)
        vx = sqrt(sum(V.mat(1:3,1:3).^2));        
    end
    
    voxdim = vx;
    bb     = BB;
    % default voxdim to current volume's voxdim, (from mat parameters)
    if any(isnan(voxdim))
        vprm = spm_imatrix(V.mat);
        vvoxdim = vprm(7:9);
        voxdim(isnan(voxdim)) = vvoxdim(isnan(voxdim));
    end
    voxdim = voxdim(:)';

    mn = bb(1,:);
    mx = bb(2,:);
    % default BB to current volume's
    if any(isnan(bb(:)))
        vbb = world_bb(V);
        vmn = vbb(1,:);
        vmx = vbb(2,:);
        mn(isnan(mn)) = vmn(isnan(mn));
        mx(isnan(mx)) = vmx(isnan(mx));
    end

    if sum(bb(:,3)) == 0
        offset = 20;
        mn(2)  =  mn(2) + offset;
        mx(2)  =  mx(2) + offset;
    end
    
    % voxel [1 1 1] of output should map to BB mn
    % (the combination of matrices below first maps [1 1 1] to [0 0 0])
    mat = spm_matrix([mn 0 0 0 voxdim])*spm_matrix([-1 -1 -1]);
    % voxel-coords of BB mx gives number of voxels required
    % (round up if more than a tenth of a voxel over)
    imgdim = ceil(mat \ [mx 1]' - 0.1)';

    % output image
    Vo(c)       = V;
    [pth,nam,ext] = fileparts(V.fname);
    Vo(c).fname      = fullfile(pth,[prefix nam ext]);
    Vo(c).dim(1:3)   = imgdim(1:3);
    Vo(c).mat        = mat;
    Vo(c) = spm_create_vol(Vo(c));
    for i = 1:imgdim(3)
        
        D = diag([-1 1 1 1]);
        if det(V.mat(1:3,1:3)) < 0
            D = diag([1 1 1 1]);
        end

        M = inv(spm_matrix([0 0 -i])*inv(Vo(c).mat)*(D*V.mat));
        img = spm_slice_vol(V, M, imgdim(1:2), deg);
        
        spm_write_plane(Vo(c), img, i);
    end
    
    c = c + 1;
end
%==========================================================================

%==========================================================================
function bb = world_bb(V)
%  world-bb -- get bounding box in world (mm) coordinates

d = V.dim(1:3);
% corners in voxel-space
c = [ 1    1    1    1
    1    1    d(3) 1
    1    d(2) 1    1
    1    d(2) d(3) 1
    d(1) 1    1    1
    d(1) 1    d(3) 1
    d(1) d(2) 1    1
    d(1) d(2) d(3) 1 ]';
% corners in world-space
tc = V.mat(1:3,1:4)*c;

% bounding box (world) min and max
mn = min(tc,[],2)';
mx = max(tc,[],2)';
bb = [mn; mx];
%==========================================================================