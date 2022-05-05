function varargout = spm_CTseg_util(varargin)
% CTseg utility functions.
%
% FORMAT y = spm_CTseg_util('affine',d,Mat)
% FORMAT id = spm_CTseg_util('identity',d)
% FORMAT spm_CTseg_util('write_nii',pth,dat,M,descrip,typ)
% FORMAT pth_out = spm_CTseg_util('modify_img_vx',pth_in,vx,odir,deg,bc)
%__________________________________________________________________________

[varargout{1:nargout}] = spm_subfun(localfunctions,varargin{:});
%==========================================================================

%==========================================================================
function pth_out = modify_img_vx(pth_in,vx,odir,deg,bc)
% Modify image voxel size.
% 
% PARAMETERS
% ----------
% pth_in : char
%     Image input path
% vx : float|[float,float,float]
%     Output image voxel size
% odir : char, default=[same directory as input image, prefixed 'r']
%     Directory where to write output image
% deg : int, default=1
%     Interpolation order
% bc : int, default=0
%     Interpolation boundary condition
% 
% RETURNS
% ----------
% pth_out : char
%     Modifed image output path
% 
%_______________________________________________________________________
if nargin < 3, odir = ''; end
if nargin < 4, deg = 1; end
if nargin < 5, bc  = 0; end

if numel(vx)  == 1, vx = vx*ones([1 3]); end
if numel(deg) == 1, deg = deg*ones([1 3]);  end
if numel(bc)  == 1, bc = bc*ones([1 3]);   end

vx(vx == 0) = 1;

% Input image properties
Nii  = nifti(pth_in);
img  = Nii.dat(:,:,:);
mat0 = Nii.mat;
vx0  = sqrt(sum(mat0(1:3,1:3).^2));
dm0  = size(img);
mn   = min(img(:));
mx   = max(img(:));
    
% Output image properties
samp = vx0./vx;
D    = diag([samp 1]);
mat  = mat0/D;
dm   = floor(D(1:3,1:3)*dm0')';

% Make interpolation grid
y = double(affine(dm,mat0\mat));

% Resample
img                 = spm_bsplinc(img, [deg bc]);
img                 = spm_bsplins(img,y(:,:,:,1),y(:,:,:,2),y(:,:,:,3),[deg bc]);    
img(~isfinite(img)) = 0;
img                 = min(mx, max(mn, img));

% Make output
[odir1,nam,ext] = fileparts(pth_in);
if isempty(odir)
    odir = odir1;
end
pth_out       = fullfile(odir,['r' nam ext]);
oNii          = nifti;
oNii.dat      = file_array(pth_out,dm,Nii.dat.dtype,Nii.dat.offset,Nii.dat.scl_slope,Nii.dat.scl_inter);
oNii.mat      = mat;
oNii.mat0     = mat;
oNii.descrip  = 'Resampled';
create(oNii);
oNii.dat(:) = img(:);
%==========================================================================

%==========================================================================
function mask_outside_fov(pth_old, pth_new, val)
% If pth_new has a larger field-of-view (fov) then pth_old, the differing
% fov voxels will be set to the value defined by val.
%
% OBS: Modifies voxel values in new image!
%
%__________________________________________________________________________
if nargin < 3, val = 0; end

% 'old' image
n_old = nifti(pth_old);
mat_old = n_old.mat;
dm_old = n_old.dat.dim(1:3);
% 'new' image
n_new = nifti(pth_new);
mat_new = n_new.mat;
dm_new = n_new.dat.dim(1:3);
% affine mapping from old to new image
mat = mat_old\mat_new;
% affine grid in new image space
grid = affine(dm_new, mat);
% make mask image 
msk = true(dm_new);
for i=1:3
    msk = msk & grid(:, :, :, i) >= 1 & grid(:, :, :, i) < dm_old(i);
end
% mask image data
img_new = n_new.dat();
img_new(~msk) = val;
% modify voxel values in new image
n_new.dat(:, :, :) = img_new;
%==========================================================================

%==========================================================================
function y = affine(d,Mat)
id    = identity(d);
y  = reshape(reshape(id,[prod(d) 3])*Mat(1:3,1:3)' + Mat(1:3,4)',[d 3]);
if d(3) == 1, y(:,:,:,3) = 1; end
%==========================================================================

%==========================================================================
function id = identity(d)
id = zeros([d(:)',3],'single');
[id(:,:,:,1),id(:,:,:,2),id(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));
%==========================================================================

%==========================================================================
function write_nii(pth,dat,M,descrip,typ)
if nargin<5, typ = 'float32'; end

spm_unlink(pth);

switch typ
case 'float32'
    fa = file_array(pth,size(dat),typ,0);
case 'uint8'
    mx = max(dat(isfinite(dat(:))));
    fa = file_array(pth,size(dat),typ,0,mx/255,0);
case 'int16'
    mx = max(dat(isfinite(dat(:))));
    mn = min(dat(isfinite(dat(:))));
    s  = max(mx/32767,-mn/32768);
    fa = file_array(pth,size(dat),typ,0,s,0);
otherwise
    error('Can''t do datatype "%s"', typ);
end

Nii         = nifti;
Nii.dat     = fa;
Nii.mat     = M;
Nii.mat0    = M;
Nii.descrip = descrip;
create(Nii);
Nii.dat(:,:,:,:,:,:) = dat;
%==========================================================================