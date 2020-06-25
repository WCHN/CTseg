function out = CTseg(in, odir, tc, def, correct_header)
% A CT segmentation+spatial normalisation routine for SPM12. 
% FORMAT out = CTseg(in, odir, tc, def)
%
%   This algorithm produces native|warped|modulated space segmentations of
%   gray matter (GM), white matter (WM) and cerebrospinal fluid (CSF). The
%   outputs are prefixed as the SPM12 unified segmentation routine.
%
% ARGS:
%  in (char|nifti): Input CT scan, either path (char array) or SPM
%                   nifti object.
%  odir (char): Directory where to write outputs, defaults to same as
%               input.
%  tc (logical(3, 3)): Tissue classes to write:
%                      [native_gm,  warped_gm,  modulated_gm;
%                       native_wm,  warped_wm,  modulated_wm;
%                       native_csf, warped_csf, modulated_csf],
%                      defaults to all true.               
%  def (logical): Write deformations? Defaults to true.
%  correct_header (logical): Correct messed up CT header, defaults to
%                            false. 
%                            OBS: This will create copy of the input image data 
%                                 and reslice it (prefixed r*)!
%
% RETURNS:
%   out - A struct with the paths to the native and template space
%         segmentations as:
%         out(1:3).c   = 'c1*.nii',   'c2*.nii',   'c3*.nii'
%         out(1:3).wc  = 'wc1*.nii',  'wc2*.nii',  'wc3*.nii'
%         out(1:3).mwc = 'mwc1*.nii', 'mwc2*.nii', 'mwc3*.nii'
%
% REFERENCE:
%   The algorithm that was used to train this model is described in the paper:
%       Brudfors M, Balbastre Y, Flandin G, Nachev P, Ashburner J. (2020). 
%       Flexible Bayesian Modelling for Nonlinear Image Registration.
%       International Conference on Medical Image Computing and Computer
%       Assisted Intervention.
%   and in the dissertation:
%       Brudfors, M. (2020). 
%       Generative Models for Preprocessing of Hospital Brain Scans.
%       Doctoral dissertation, UCL (University College London).
%   Please consider citing if you find this code useful.A more detailed
%   paper validating the method will hopefully be published soon.
%
% AUTHOR:
%   Mikael Brudfors, brudfors@gmail.com, 2020
%_______________________________________________________________________

if nargin < 2, odir = ''; end
if nargin < 3, tc = true(3, 3); end
if nargin < 4, def = true; end
if nargin < 5, correct_header = false; end
if size(tc,1) == 1
    tc = repmat(tc, 3, 1);
end

if ~(exist(fullfile(fileparts(mfilename('fullpath')),'spm_mb_model.mat'), 'file') == 2)
    % Unzip model file, if has not been done
    unzip(fullfile(fileparts(mfilename('fullpath')),'model.zip'));
end

% Check MATLAB path
if isempty(fileparts(which('spm'))),         error('SPM12 not on the MATLAB path!'); end % download from https://www.fil.ion.ucl.ac.uk/spm/software/download/
if isempty(fileparts(which('spm_mb_fit'))),  error('diffeo-segment not on the MATLAB path!'); end % git clone from https://github.com/WTCN-computational-anatomy-group/diffeo-segment
if isempty(fileparts(which('spm_gmm_lib'))), error('auxiliary-functions not on the MATLAB path!'); end % git clone from https://github.com/WTCN-computational-anatomy-group/auxiliary-functions

% Get nifti
if ~isa(in,'nifti'), Nii = nifti(in);
else,                Nii = in;
end; clear in

% Correct orientation matrix
if correct_header
    Nii = correct_orientation(Nii);
end

% Get image data
F = [];
for i=1:numel(Nii)
    F = cat(4,F,single(Nii(i).dat()));
end

% Quick, rough rigid alignment to MNI
R = Rigid2MNI(Nii.dat.fname);

% Get spm_vol (to integrate R)
V     = spm_vol(Nii.dat.fname);
M0    = V.mat;
V.mat = R\V.mat;

% Get dat struct (for spm_mb)
dat = struct('F',F,'V',V,'is_ct',true,'do_dc',false);

% Settings
sett = struct;
sett.model.init_mu_dm = 8;
sett.nit.init         = 16;
sett.var.v_settings   = [0 0 0.2 0.05 0.2]*4;
% Write options
if isempty(odir)
    odir = fileparts(Nii.dat.fname);
    s    = what(odir); % Get absolute path
    odir = s.path;
end
sett.write.dir_res              = odir;
sett.write.tc                   = false(9,3);
sett.write.tc([1 2 3], [1 2 3]) = tc;
sett.write.df                   = def;
sett.clean_z.mrf                = 1;
sett.clean_z.gwc_tix            = struct('gm',[3],'wm',[2 5],'csf',[4]);
sett.write.vx                   = 1.5;
sett.write.bb                   = NaN(2,3);
sett.write.bb                   = [-90 -126 -72; 90 90 108];
% % Uncomment for testing
% sett.show.figs = {'model','segmentations'};
% sett.nit.init = 1;
% sett.nit.init_mu = 1;
% sett.nit.zm = 1;
% sett.model.init_mu_dm = 32;

% Path to spm_mb model file
PthModel = fullfile('spm_mb_model.mat');

% If SPM has been compiled with OpenMP support then the number of threads
% are here set to speed up the algorithm.
setenv('SPM_NUM_THREADS',sprintf('%d',-1));

% Run Register
[dat,mu,sett] = spm_mb_fit(dat,'PthModel',PthModel,'sett',sett);

% Write results in normalised space
res = spm_mb_output(dat,mu,sett);

% Reset orientation matrix
for i=1:numel(Nii), spm_get_space(res.c{i},M0); end

% Make output
cl = cell(1, 3);
out = struct('c', cl, 'wc', cl, 'mwc', cl);
for k=1:3        
    if ~isempty(res.c)      
        out(k).c = res.c{k};
    end
    if ~isempty(res.wc)      
        out(k).wc = res.wc{k};
    end
    if ~isempty(res.mwc)      
        out(k).mwc = res.mwc{k};
    end
end
%==========================================================================

%==========================================================================
function Nii = correct_orientation(Nii)
f = nm_reorient(Nii.dat.fname);
reset_origin(f);
Nii = nifti(f);
%==========================================================================

%==========================================================================
function Mout = reset_origin(pth)
V   = spm_vol(pth);
M   = V.mat;
dim = V.dim;
vx  = sqrt(sum(M(1:3,1:3).^2));
if det(M(1:3,1:3))<0
    vx(1) = -vx(1); 
end
orig = (dim(1:3)+1)/2;
off  = -vx.*orig;
M1   = [vx(1) 0      0         off(1)
           0      vx(2) 0      off(2)
           0      0      vx(3) off(3)
           0      0      0      1];
V    = spm_vol(pth);
M0   = V.mat;
Mout = M0/M1;
spm_get_space(pth,M1);   
%==========================================================================

function npth = nm_reorient(pth,vx,prefix,deg)
if nargin < 2, vx     = [];   end
if nargin < 3, prefix = 'r'; end
if nargin < 4, deg    = 1;    end

if ~isempty(vx) && length(vx) < 3
    vx=[vx vx vx];
end

% Get information about the image volumes
VV = spm_vol(pth);

for V=VV' % Loop over images

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
    if spm_flip_analyze_images, tc(1,:) = -tc(1,:); end

    % Max and min co-ordinates for determining a bounding-box
    mx = round(max(tc,[],2)');
    mn = round(min(tc,[],2)');

    vx0 = sqrt(sum(V.mat(1:3,1:3).^2));
    if isempty(vx)
        vx = vx0;
    end    
    
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
    if spm_flip_analyze_images, mat = diag([-1 1 1 1])*mat; end
    VO.mat           = mat;

    % Create .hdr and open output .img
    VO = spm_create_vol(VO);

    for i=1:dim(3) % Loop over slices of output image

        % Mapping from slice i of the output image,
        % to voxels of the input image
        M   = inv(spm_matrix([0 0 -i])*inv(VO.mat)*V.mat);

        % Extract this slice according to the mapping
        img = spm_slice_vol(V,M,dim(1:2),deg);

        % Write this slice to output image
        spm_write_plane(VO,img,i);
    end % End loop over output slices

end % End loop over images
npth = VO.fname;
%==========================================================================

%==========================================================================
function R = Rigid2MNI(P)
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

% Affine = spm_maff8(P,2,32,tpm,Affine,'mni'); % Heavily regularised
% Affine = spm_maff8(P,2,1 ,tpm,Affine,'mni'); % Lightly regularised

% % Load header
% Nii    = nifti(P);

% Generate mm coordinates of where deformations map from
x      = affind(rgrid(size(tpm.dat{1})),tpm.M);

% Generate mm coordinates of where deformation maps to
y1     = affind(x,inv(Affine));

% Weight the transform via GM+WM
weight = single(exp(tpm.dat{1})+exp(tpm.dat{2}));

% Weighted Procrustes analysis
[~,R]  = spm_get_closest_affine(x,y1,weight);
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