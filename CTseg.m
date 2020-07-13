function seg = CTseg(in, odir, tc, def, correct_header)
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
%  def (logical): Write deformations? Defaults to false.
%  correct_header (logical): Correct messed up CT header, defaults to
%                            false. 
%                            OBS: This will create copy of the input image data 
%                                 and reslice it (prefixed r*)!
%
% RETURNS:
%   seg - A struct with the paths to the native and template space
%         segmentations as:
%         seg(1:3).c   = 'c1*.nii',   'c2*.nii',   'c3*.nii'
%         seg(1:3).wc  = 'wc1*.nii',  'wc2*.nii',  'wc3*.nii'
%         seg(1:3).mwc = 'mwc1*.nii', 'mwc2*.nii', 'mwc3*.nii'
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
if nargin < 4, def = false; end
if nargin < 5, correct_header = false; end
if size(tc,1) == 1
    tc = repmat(tc, 3, 1);
end

% Get model files
%--------------------------------------------------------------------------
ctseg_dir = fileparts(mfilename('fullpath'));
if ~(exist(fullfile(ctseg_dir,'mu_CTseg.nii'), 'file') == 2)
    % Path to model zip file
    pth_model_zip = fullfile(ctseg_dir, 'model.zip');    
    % Model file not present
    if ~(exist(fullfile(ctseg_dir,'model.zip'), 'file') == 2)
        % Download model file
        url_model = 'https://www.dropbox.com/s/clatdvlbih97mr1/model.zip?dl=1';
        fprintf('Downloading model files (first use only)... ')
        websave(pth_model_zip, url_model);                
        fprintf('done.\n')
    end    
    % Unzip model file, if has not been done
    fprintf('Extracting model files  (first use only)... ')
    unzip(pth_model_zip, ctseg_dir);
    fprintf('done.\n')
end

% Check MATLAB path
%--------------------------------------------------------------------------
if isempty(fileparts(which('spm'))),         error('SPM12 not on the MATLAB path!'); end % download from https://www.fil.ion.ucl.ac.uk/spm/software/download/
if isempty(fileparts(which('spm_mb_fit'))),  error('Multi-Brain not on the MATLAB path!'); end % download/clone from https://github.com/WTCN-computational-anatomy-group/diffeo-segment

% Get nifti
%--------------------------------------------------------------------------
if ~isa(in,'nifti'), Nii = nifti(in);
else,                Nii = in;
end; clear in

% Correct orientation matrix
%--------------------------------------------------------------------------
if correct_header
    Nii = correct_orientation(Nii);
end

% Get model file paths
%--------------------------------------------------------------------------
pth_mu = fullfile(ctseg_dir,'mu_CTseg.nii');
if ~(exist(pth_mu, 'file') == 2)
    error('Atlas file (mu_CTseg.nii) could not be found! Has model.zip not been unzipped?')
end
pth_int_prior = fullfile(ctseg_dir,'prior_CTseg.mat');
if ~(exist(pth_int_prior, 'file') == 2)
    error('Intensity prior file (pth_int_prior.mat) could not be found! Has model.zip not been unzipped?')
end

% Output directory
%--------------------------------------------------------------------------
if isempty(odir)
    odir = fileparts(Nii.dat.fname);
    s    = what(odir); % Get absolute path
    odir = s.path;
end

% Settings
%--------------------------------------------------------------------------
% spm_mb_fit
run = struct;
run.mu.exist = {pth_mu};
run.aff = 'SE(3)';
run.v_settings = 0.5*[0.0001 0 0.4 0.1 0.4];
run.onam = 'mb';
run.odir = {odir};
run.cat = {{}};
run.gmm.chan.images = {Nii(1).dat.fname};
run.gmm.chan.inu.inu_reg = 10000;
run.gmm.chan.inu.inu_co = 40;
run.gmm.chan.modality = 2;
run.gmm.labels.false = [];
run.gmm.pr.file = {pth_int_prior};
run.gmm.pr.hyperpriors = [];
run.gmm.tol_gmm = 0.0005;
run.gmm.nit_gmm_miss = 32;
run.gmm.nit_gmm = 8;
run.gmm.nit_appear = 4;
run.accel = 0.8;
run.min_dim = 16;
run.tol = 0.001;
run.sampdens = 2;
run.save = false;
run.nworker = 0;
% spm_mb_output
out = struct;
out.i = false;
out.mi = false;
out.wi = false;
out.wmi = false;
out.inu = false;
out.c = find(tc(:,1) > 0);
out.wc = find(tc(:,2) > 0);
out.mwc = find(tc(:,3) > 0);
out.v = false;
out.y = def;
out.mrf = 1;

% Run segmentation+normalisation
%--------------------------------------------------------------------------
% Init MB
[dat,sett] = spm_mb_init(run);
if ~isempty(dat)
    % Fit MB
    [dat,sett] = spm_mb_fit(dat,sett);
    dat        = spm_mb_io('save_psi',dat,sett);
    % Save results
    p_res = fullfile(sett.odir,['mb_fit_' sett.onam '.mat']);
    save(p_res,'dat','sett');
    out.result = p_res;
    % Write output
    res = spm_mb_output(out);
    delete(p_res);
end

% Format output
%--------------------------------------------------------------------------
cl = cell(1, 3);
seg = struct('c', cl, 'wc', cl, 'mwc', cl);
for k=1:3        
    if ~isempty(out.c) && out.c(k)
        seg(k).c = res.c{k};
    end
    if ~isempty(out.wc) && out.wc(k)     
        seg(k).wc = res.wc{k};
    end
    if ~isempty(out.mwc) && out.mwc(k)    
        seg(k).mwc = res.mwc{k};
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