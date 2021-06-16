function [res,vol] = spm_CTseg(in, odir, tc, def, correct_header, skullstrip, vox)
% A CT segmentation+spatial normalisation routine for SPM12. 
% FORMAT [res,vol] = spm_CTseg(in, odir, tc, def, correct_header, skullstrip, vox)
%
% This algorithm produces native|warped|modulated space segmentations of:
%     1. Gray matter (GM) hemisphere 1
%     2. Gray matter hemisphere 2
%     3. White matter (WM) hemisphere 1
%     4. White matter hemisphere 2
%     5. Cerebrospinal fluid (CSF)
%     6. Bone (BONE)
%     7. Soft tissue (ST)
%     8. Background (BG)
% the outputs are prefixed as the SPM12 unified segmentation (c*, wc*, mwc*).
%
% ARGS:
% --------------
% in (char|nifti): Input CT scan, either path (char array) or SPM
%                  nifti object.
%
% odir (char): Directory where to write outputs, defaults to same as
%              input CT scan.
%
% tc (logical(8, 3)): Matrix where native, warped and warped modulated are
%                     indexed by columns and tissue classes are indexed by rows 
%                     (in the above order).             
%
% def (logical): Write deformations? Defaults to true.
%
% correct_header (logical): Correct messed up CT header, defaults to true. 
%
% skullstrip (logical): Write skull-stripped CT scan to disk, prefixed 
%                       'ss_'. Defaults to false.
%
% vox (double): Template space voxel size, defaults to voxel size of
%               template.
%
% RETURNS:
% --------------
% res - A struct with paths to algorithm results.
%
% vol - A struct with total brain and intercranial volume (TBV and TIV), in
%       millilitres.
%
% REFERENCES:
% --------------
% The algorithm that was used to train this model is described in the paper:
%
%     Brudfors M, Balbastre Y, Flandin G, Nachev P, Ashburner J. (2020). 
%     Flexible Bayesian Modelling for Nonlinear Image Registration.
%     International Conference on Medical Image Computing and Computer
%     Assisted Intervention.
%
% and in the PhD dissertation:
%
%     Brudfors, M. (2020). 
%     Generative Models for Preprocessing of Hospital Brain Scans.
%     Doctoral dissertation, UCL (University College London).
%
% Please consider citing if you find this code useful. A more detailed
% paper validating the method will hopefully be published soon.
%
% CONTACT:
% --------------
% Mikael Brudfors, brudfors@gmail.com, 2020
%_______________________________________________________________________

if ~nargin
  spm_jobman('interactive','','spm.tools.CTseg');
  return;
end
  
if nargin < 2, odir = ''; end
if nargin < 3, tc   = true; end
<<<<<<< HEAD
=======
end
>>>>>>> hemisphere
if size(tc,2) == 1
    tc = repmat(tc, 1, 3);
end
if nargin < 4, def            = true; end
if nargin < 5, correct_header = true; end
if nargin < 6, skullstrip     = false; end
if nargin < 7, vox            = NaN; end

% check MATLAB path
%--------------------------------------------------------------------------
if isempty(fileparts(which('spm')))
    error('SPM12 not on the MATLAB path! Download from https://www.fil.ion.ucl.ac.uk/spm/software/download/'); 
end
if isempty(fileparts(which('spm_shoot3d')))
    error('Shoot toolbox not on the MATLAB path! Add from spm12/toolbox/Shoot'); 
end
if isempty(fileparts(which('spm_dexpm')))
    error('Longitudinal toolbox not on the MATLAB path! Add from spm12/toolbox/Longitudinal'); 
end
% add MB toolbox
addpath(fullfile(spm('dir'),'toolbox','mb')); 
if isempty(fileparts(which('spm_mb_fit')))
    error('Multi-Brain toolbox not on the MATLAB path! Download/clone from https://github.com/WTCN-computational-anatomy-group/mb and place in the SPM12 toolbox folder.');
end
if ~(exist('spm_gmmlib','file') == 3)
    error('Multi-Brain GMM library is not compiled, please follow the Install instructions on the Multi-Brain GitHub README.')
end

% Get model files
%--------------------------------------------------------------------------
ctseg_dir = fileparts(mfilename('fullpath'));
if ~(exist(fullfile(ctseg_dir,'mu_CTseg.nii'), 'file') == 2)
    % Path to model zip file
    pth_model_zip = fullfile(ctseg_dir, 'model.zip');    
    % Model file not present
    if ~(exist(pth_model_zip, 'file') == 2)
        % Download model file
        url_model = 'https://www.dropbox.com/s/qjdqavysgqqhyzc/model.zip?dl=1';
        fprintf('Downloading model files (first use only)... ')
        websave(pth_model_zip, url_model);                
        fprintf('done.\n')
    end    
    % Unzip model file, if has not been done
    fprintf('Extracting model files  (first use only)... ')
    unzip(pth_model_zip, ctseg_dir);
    fprintf('done.\n')    
    % Delete model.zip
    spm_unlink(pth_model_zip);
end

% Get nifti
%--------------------------------------------------------------------------
Nii = nifti(in);

% Output directory
%--------------------------------------------------------------------------
if isempty(odir)
    odir = fileparts(Nii.dat.fname);
    s    = what(odir); % Get absolute path
    odir = s.path;
elseif ~(exist(odir, 'dir') == 7)
    mkdir(odir)    
end

% Correct orientation matrix
%--------------------------------------------------------------------------
Mc = eye(4);
oNii = Nii;
if correct_header
    [Nii,Mc] = correct_orientation(Nii, odir);
end

% Get model file paths
%--------------------------------------------------------------------------
pth_mu = fullfile(ctseg_dir,'mu_CTseg.nii');
if ~(exist(pth_mu, 'file') == 2)
    error('Atlas file (mu_CTseg.nii) could not be found! Has model.zip not been extracted?')
end
pth_int = fullfile(ctseg_dir,'prior_CTseg.mat');
if ~(exist(pth_int, 'file') == 2)
    error('Intensity prior file (pth_int_prior.mat) could not be found! Has model.zip not been extracted?')
end
pth_Mmni = fullfile(ctseg_dir,'Mmni.mat');
if ~(exist(pth_Mmni, 'file') == 2)
    error('MNI affine (Mmni.mat) could not be found! Has model.zip not been extracted?')
end
% Get number of tissue classes from template
Nii_mu = nifti(pth_mu);
K      = Nii_mu.dat.dim(4) + 1;
if K == 6
    % Modify default model to separate GM and WM into hemispheres
    % get template data
    mu = single(Nii_mu.dat());
    mu = mu(:,:,:,[1, 1, 2, 2, 3, 4, 5]);
    % separate hemispheres
    ix0                   = floor(0.5*size(mu,1));
    ix1                   = ceil(0.5*size(mu,1));
    min_mu                = min(mu(:));
    mu(1:ix0,:,:,[1 3])   = min_mu;
    mu(ix1:end,:,:,[2 4]) = min_mu;
    % write modified template
    write_nii(pth_mu,mu,Nii_mu.mat,'Hemisphere template','float32');
    % load intensity prior
    pr    = load(pth_int);
    mg_ix = pr.mg_ix;
    pr    = pr.pr;
    % modify number of Gaussians
    mg_ix = [1 2 3 4 mg_ix(3:end) + 2];
    ix    = [1 1 2 2 3:size(pr{1},2)];
    pr{1} = pr{1}(:,ix);
    pr{2} = pr{2}(:,ix);
    pr{3} = pr{3}(:,:,ix);
    pr{4} = pr{4}(:,ix);
    pr{5} = pr{5}(:,ix);
    % save modified prior
    save(pth_int,'pr','mg_ix');
    % new number of classes
    K = size(mu,4) + 1;
end
if size(tc,1) == 1
    tc = repmat(tc, K, 1);
end
% Indices of GM, WM and CSF classes in template
ix_gm  = [1 2];
ix_wm  = [3 4];
ix_csf = [5];
ix_gw  = [ix_gm, ix_wm];
ix_gwc = [ix_gm, ix_wm, ix_csf];

% Run MB
%--------------------------------------------------------------------------
% algorithm settings
run            = struct;
run.mu.exist   = {pth_mu};
run.onam       = 'CTseg';
run.odir       = {odir};    
run.v_settings = [0.00001 0 0.4 0.1 0.4]*1;
run.min_dim    = 8;
run.tol        = 0.5*0.001;
% image
run.gmm.pr.file          = {pth_int};
run.gmm.pr.hyperpriors   = [];
run.gmm.chan.images      = {Nii(1).dat.fname};
run.gmm.chan.modality    = 2;
run.gmm.chan.inu.inu_reg = 1e6;
% output settings
out         = struct;
out.result  = {fullfile(run.odir{1},['mb_fit_' run.onam '.mat'])};
out.c       = 1:K;
out.wc      = find(tc(:,2))';
out.mwc     = find(tc(:,3))';
out.vox     = vox;
out.proc_zn = {@(x) clean_gwc(x, struct('gm',ix_gm,'wm',ix_wm,'csf',ix_csf))};

% fit model and write output
jobs{1}.spm.tools.mb.run = run;
jobs{2}.spm.tools.mb.out = out;
res                      = spm_jobman('run', jobs);

% get results
res     = load(res{1}.fit{1});
dat     = res.dat;
Mmu     = res.sett.mu.Mmu;
res.c   = cell(1,K);
res.wc  = cell(1,K);
res.mwc = cell(1,K);
for k=1:K
    res.c{k} = fullfile(dat.odir, ['c0' num2str(k) '_' dat.onam '.nii']);
    if tc(k,2)
        res.wc{k} = fullfile(dat.odir, ['wc0' num2str(k) '_' dat.onam '.nii']);
    end
    if tc(k,3)
        res.mwc{k} = fullfile(dat.odir, ['mwc0' num2str(k) '_' dat.onam '.nii']);
    end
end

% Reslice template space segmentations to MNI space
reslice2mni(res,pth_Mmni,Mmu);

if nargout > 1 || skullstrip
    % Get responsibilities
    Z = [];
    for k=1:K
        Nii_c = nifti(res.c{k});
        Z     = cat(4, Z, single(Nii_c.dat()));
    end    
end

vol = struct('tbv',NaN,'tiv',NaN);
if nargout > 1
    % Compute TBV and TIV    
    vx      = sqrt(sum(Nii(1).mat(1:3,1:3).^2));
    vol.tbv = prod(vx(1:3))*sum(sum(sum(sum(Z(:,:,:,ix_gw)))));
    vol.tiv = prod(vx(1:3))*sum(sum(sum(sum(Z(:,:,:,ix_gwc)))));
end

if correct_header
    % Reslice corrected native space segmentations to original native space.
    M1  = spm_get_space(Nii(1).dat.fname);  % get corrected orientation matrix
    spm_unlink(Nii(1).dat.fname);           % delete corrected image
    Nii = oNii;                             % reset to original input image
    M0  = spm_get_space(Nii(1).dat.fname);  % get corrected orientation matrix
    % new field-of-view
    M = Mc\M1\M0;
    y = affine(Nii.dat.dim,M);     
    % reslice segmentations
    for k=1:K
        if isempty(res.c{k}), continue; end                        
        Nii_c = nifti(res.c{k});        
        rc    = spm_diffeo('bsplins',single(Nii_c.dat()),y,[1 1 1  0 0 0]);
        write_nii(res.c{k},rc,M0,sprintf('Tissue (%d)',k), 'uint8')        
    end
end

res.s = '';
if skullstrip
    % Produce skull-stripped CT scan (prefixed 'ss_')
    %----------------------------------------------------------------------
    if correct_header
        % Get resliced responsibilities
        Z = [];
        for k=1:K
            Nii_c = nifti(res.c{k});
            Z     = cat(4, Z, single(Nii_c.dat()));
        end
        Z = bsxfun(@rdivide, Z, sum(Z,4) + eps('single'));  % renormalise
    end
    % Copy image
    [~,nam,ext] = fileparts(Nii(1).dat.fname);
    nfname      = fullfile(run.odir{1},['ss_' nam ext]);
    copyfile(Nii(1).dat.fname,nfname);
    % Make mask and apply
    Nii_s     = nifti(nfname);
    img       = single(Nii_s.dat());
    msk       = sum(Z(:,:,:,ix_gwc),4) >= 0.5;
    img(~msk) = 0;
    % Modify copied image's data
    Nii_s.dat(:,:,:) = img;
    res.s            = nfname;
    clear msk
end
clear Z

% Delete unrequested native space segmentations
res_c = res.c;
res.c = cell(1,sum(tc(:,1)));
k1 = 1;
for k=1:K
    if ~tc(k,1)
        spm_unlink(res_c{k});
    else
        res.c{k1} = res_c{k};
        k1        = k1 + 1;
    end
end

% Save deformation?
res.y = '';
if def
    res.y = dat(1).psi.dat.fname;
else
    spm_unlink(dat(1).psi.dat.fname);    
end
spm_unlink(dat(1).v.dat.fname); % Delete velocity field
spm_unlink(fullfile(run.odir{1},['mb_fit_' run.onam '.mat'])); % Delete mb_fit_mb.mat
%==========================================================================

%==========================================================================
function [Nii,Mr] = correct_orientation(Nii,odir)
f   = nm_reorient(Nii.dat.fname,odir);
Mr  = reset_origin(f);
Nii = nifti(f);
%==========================================================================

%==========================================================================
function Mr = reset_origin(pth)
V   = spm_vol(pth);
M0  = V.mat;
dim = V.dim;
vx  = sqrt(sum(M0(1:3,1:3).^2));
if det(M0(1:3,1:3))<0
    vx(1) = -vx(1); 
end
orig = (dim(1:3)+1)/2;
off  = -vx.*orig;
M1   = [vx(1)     0     0 off(1)
           0  vx(2)     0 off(2)
           0      0 vx(3) off(3)
           0      0     0      1];
Mr = M1/M0;
spm_get_space(pth,Mr*M0); 
%==========================================================================

%==========================================================================
function npth = nm_reorient(pth,odir,vx,prefix,deg)
if nargin < 3, vx     = [];   end
if nargin < 4, prefix = 'temp_'; end
if nargin < 5, deg    = 1;    end

if ~isempty(vx) && length(vx) < 3
    vx=[vx vx vx];
end

% Get information about the image volumes
V = spm_vol(pth);

% The corners of the current volume
d = V.dim(1:3);
c = [1    1    1    1
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
[~,name,ext] = fileparts(V.fname);
VO.fname     = fullfile(odir,[prefix name ext]);

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
npth = VO.fname;
%==========================================================================

%==========================================================================
function zn = clean_gwc(zn,ixt,level)
if nargin < 2 || isempty(ixt)
    % Default SPM12 template ordering
    ixt = struct('gm',1,'wm',2,'csf',3);
end
if nargin < 3, level = 1; end

b = sum(zn(:,:,:,ixt.wm),4);

% Build a 3x3x3 seperable smoothing kernel
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

% Erosions and conditional dilations
th1 = 0.15;
if level==2, th1 = 0.2; end
niter  = 32;
niter2 = 32;
for j=1:niter
    if j>2
        th       = th1;
    else
        th       = 0.6;
    end  % Dilate after two its of erosion
    for i=1:size(b,3)
        gp       = double(sum(zn(:,:,i,ixt.gm),4));
        wp       = double(sum(zn(:,:,i,ixt.wm),4));
        bp       = double(b(:,:,i));
        bp       = (bp>th).*(wp+gp);
        b(:,:,i) = bp;
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
end

% Also clean up the CSF.
if niter2 > 0
    c = b;
    for j=1:niter2
        for i=1:size(b,3)
            gp       = double(sum(zn(:,:,i,ixt.gm),4));
            wp       = double(sum(zn(:,:,i,ixt.wm),4));
            cp       = double(sum(zn(:,:,i,ixt.csf),4));
            bp       = double(c(:,:,i));
            bp       = (bp>th).*(wp+gp+cp);
            c(:,:,i) = bp;
        end
        spm_conv_vol(c,c,kx,ky,kz,-[1 1 1]);
    end
end

th = 0.05;
for i=1:size(b,3)
    slices = cell(1,size(zn,4));
    for k1=1:size(zn,4)
        slices{k1} = double(zn(:,:,i,k1));
    end
    bp           = double(b(:,:,i));
    bp           = ((bp>th).*(sum(cat(3,slices{ixt.gm}),3)+sum(cat(3,slices{ixt.wm}),3)))>th;
    for i1=1:numel(ixt.gm)
        slices{ixt.gm(i1)} = slices{ixt.gm(i1)}.*bp;
    end
    for i1=1:numel(ixt.wm)
        slices{ixt.wm(i1)} = slices{ixt.wm(i1)}.*bp;
    end

    if niter2>0
        cp           = double(c(:,:,i));
        cp           = ((cp>th).*(sum(cat(3,slices{ixt.gm}),3)+sum(cat(3,slices{ixt.wm}),3)+sum(cat(3,slices{ixt.csf}),3)))>th;

        for i1=1:numel(ixt.csf)
            slices{ixt.csf(i1)} = slices{ixt.csf(i1)}.*cp;
        end
    end
    tot       = zeros(size(bp))+eps;
    for k1=1:size(zn,4)
        tot   = tot + slices{k1};
    end
    for k1=1:size(zn,4)
        zn(:,:,i,k1) = slices{k1}./tot;
    end
end
%==========================================================================

%==========================================================================
function reslice2mni(res, pth_Mmni, Mmu)
% Load affine matrix that aligns MB template with SPM template
load(pth_Mmni, 'Mmni');
% Get SPM template information
Niis = nifti(fullfile(spm('Dir'),'tpm','TPM.nii'));
Ms   = Niis.mat;
ds   = Niis.dat.dim(1:3);
vxs  = sqrt(sum(Ms(1:3,1:3).^2));
% Extract affine transformation from spm_klaff result
Md = inv(Ms\Mmni);
A  = Mmu*Md/Ms;
% Do reslice
if ~isempty(res.wc)
    for k=1:numel(res.wc)
        if isempty(res.wc{k}), continue; end
        reslice_dat(res.wc{k},A,Mmu,Ms,ds,vxs,'uint8');
    end
end
if ~isempty(res.mwc)
    for k=1:numel(res.mwc)
        if isempty(res.mwc{k}), continue; end
        reslice_dat(res.mwc{k},A,Mmu,Ms,ds,vxs,'int16');
    end
end   
%==========================================================================

%==========================================================================
function pth = reslice_dat(pth,A,Mmu0,Ms,ds,vxs,typ)
% Get template-space orientation matrix 
% (possibly with cropped FOV and adjusted voxel size)
Mmu = spm_get_space(pth);
% Get cropping matrix
Mc = Mmu/Mmu0;
% New field of view
vx_out = sqrt(sum(Mmu(1:3,1:3).^2));
D      = diag([vxs./vx_out 1]);
Mout   = Ms/D;
dout   = floor(D(1:3,1:3)*ds')';
% Define sampling grid
M = (Mc*Mmu0)\A*Mout;
y = affine(dout,M);
% Reslice
Nii   = nifti(pth);
dat   = spm_diffeo('bsplins',single(Nii.dat()),y,[1 1 1  0 0 0]);
write_nii(pth,dat,Mout,Nii.descrip,typ);
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

%==========================================================================
function psi0 = affine(d,Mat)
id    = identity(d);
psi0  = reshape(reshape(id,[prod(d) 3])*Mat(1:3,1:3)' + Mat(1:3,4)',[d 3]);
if d(3) == 1, psi0(:,:,:,3) = 1; end
%==========================================================================

%==========================================================================
function id = identity(d)
id = zeros([d(:)',3],'single');
[id(:,:,:,1),id(:,:,:,2),id(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));
%==========================================================================
