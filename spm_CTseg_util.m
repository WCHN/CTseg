function varargout = spm_CTseg_util(varargin)
% CTseg utility functions.
%
% FORMAT y = spm_CTseg_util('affine',d,Mat)
% FORMAT id = spm_CTseg_util('identity',d)
% FORMAT spm_CTseg_util('write_nii',pth,dat,M,descrip,typ)
%__________________________________________________________________________

[varargout{1:nargout}] = spm_subfun(localfunctions,varargin{:});
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