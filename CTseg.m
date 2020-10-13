function [res,vol] = CTseg(in, odir, tc, def, correct_header, mni, skullstrip)
% A CT segmentation+spatial normalisation routine for SPM12. 
% FORMAT [res,vol] = CTseg(in, odir, tc, def, correct_header, mni, skullstrip)
%
% This algorithm produces native|warped|modulated space segmentations of:
%     1. Gray matter (GM)
%     2. White matter (WM)
%     3. Meninges, sinuses, calcifications (MEN)
%     4. Bone (BONE)
%     5. Bone (BONE)
%     6. Soft tissue (ST)
%     7. Background (BG),
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
% tc (logical(7, 3)): Matrix where native, warped and warped modulated are
%                     indexed by columns and tissue classes are indexed by rows 
%                     (in the above order).             
%
% def (logical): Write deformations? Defaults to true.
%
% correct_header (logical): Correct messed up CT header, defaults to
%                            true. 
%                            OBS: This will create a copy of the input image 
%                                 data and reslice it (prefixed r*)!
%
% mni (logical): Should normalised space be in MNI? Defaults to true.
%
% skullstrip (logical): Write skull-stripped CT scan to disk, prefixed 
%                        's_'. Defaults to true.
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
% Please consider citing if you find this code useful.A more detailed
% paper validating the method will hopefully be published soon.
%
% CONTACT:
% --------------
% Mikael Brudfors, brudfors@gmail.com, 2020
%_______________________________________________________________________

if nargin < 2, odir = ''; end
if nargin < 3, tc = [[true(3,1); false(4,1)], ...
                     [true(2,1); false(5,1)], ...
                     [true(2,1); false(5,1)]]; 
end
if nargin < 4, def = true; end
if nargin < 5, correct_header = true; end
if nargin < 6, mni = true; end
K = 7; % Number of segmentation classes
if size(tc,1) == 1
    tc = repmat(tc, K, 1);
end
if nargin < 7, skullstrip = true; end

% Check MATLAB path
%--------------------------------------------------------------------------
if isempty(fileparts(which('spm'))), error('SPM12 not on the MATLAB path! Download from https://www.fil.ion.ucl.ac.uk/spm/software/download/'); end
if isempty(fileparts(which('spm_mb_fit'))),error('Multi-Brain not on the MATLAB path! Download/clone from https://github.com/WTCN-computational-anatomy-group/diffeo-segment'); end
if isempty(fileparts(which('spm_shoot3d'))), error('Shoot toolbox not on the MATLAB path! Add from spm/toolbox/Shoot'); end
if isempty(fileparts(which('spm_dexpm'))), error('Longitudinal toolbox not on the MATLAB path! Add from spm/toolbox/Longitudinal'); end

% Get model files
%--------------------------------------------------------------------------
ctseg_dir = fileparts(mfilename('fullpath'));
if ~(exist(fullfile(ctseg_dir,'mu_CTseg.nii'), 'file') == 2)
    % Path to model zip file
    pth_model_zip = fullfile(ctseg_dir, 'model.zip');    
    % Model file not present
    if ~(exist(fullfile(ctseg_dir,'model.zip'), 'file') == 2)
        % Download model file
        url_model = 'https://www.dropbox.com/s/bi0r2t6lcmcl61q/model.zip?dl=1';
        fprintf('Downloading model files (first use only)... ')
        websave(pth_model_zip, url_model);                
        fprintf('done.\n')
    end    
    % Unzip model file, if has not been done
    fprintf('Extracting model files  (first use only)... ')
    unzip(pth_model_zip, ctseg_dir);
    fprintf('done.\n')
end

% Get nifti
%--------------------------------------------------------------------------
if ~isa(in,'nifti'), Nii = nifti(in);
else,                Nii = in;
end; clear in

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
Mr = eye(4);
if correct_header
    [Nii,Mr] = correct_orientation(Nii, odir);
end

% Get model file paths
%--------------------------------------------------------------------------
pth_mu = fullfile(ctseg_dir,'mu_CTseg.nii');
if ~(exist(pth_mu, 'file') == 2)
    error('Atlas file (mu_CTseg.nii) could not be found! Has model.zip not been extracted?')
end
pth_int_prior = fullfile(ctseg_dir,'prior_CTseg.mat');
if ~(exist(pth_int_prior, 'file') == 2)
    error('Intensity prior file (pth_int_prior.mat) could not be found! Has model.zip not been extracted?')
end
if mni
    pth_Mmni = fullfile(ctseg_dir,'Mmni.mat');
    if ~(exist(pth_Mmni, 'file') == 2)
        error('MNI affine (Mmni.mat) could not be found! Has model.zip not been extracted?')
    end
end

% Settings
%--------------------------------------------------------------------------
% spm_mb_fit
run = struct;
run.mu.exist = {pth_mu};
run.aff = 'SE(3)';
run.v_settings = [0.0001 0 0.4 0.1 0.4]*2;
run.onam = 'mb';
run.odir = {odir};
run.cat = {{}};
run.gmm.chan.images = {Nii(1).dat.fname};
run.gmm.chan.inu.inu_reg = 1e4;
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
out.c = true(1,K);
out.wc = false(1,K);
out.mwc = false(1,K);
out.sm = false(1,K);
out.v = false;
out.y = true;
out.mrf = 0;

% Run segmentation+normalisation
%--------------------------------------------------------------------------
% Init MB
[dat,sett] = spm_mb_init(run);

% Fit MB
[dat,sett] = spm_mb_fit(dat,sett);

% Save results
p_res = fullfile(sett.odir,['mb_fit_' sett.onam '.mat']);
save(p_res,'dat','sett');
out.result = p_res;

% Write output
res1    = spm_mb_output(out);    
res     = struct;
res.c1  = res1.c;
clear res1

% Ad-hoc clean-up..
%--------------------------------------------------------------------------
Z = [];
for k=1:K
    Nii_c = nifti(res.c1{k});
    Z     = cat(4, Z, single(Nii_c.dat()));
end
Z = cat(4, Z, 1 - sum(Z,4));
Z = clean_gwc(Z,struct('gm',1,'wm',2,'csf',3),2);
for k=1:K
    Nii_c            = nifti(res.c1{k});
    Nii_c.dat(:,:,:) = Z(:,:,:,k);
end

% Compute TBV and TIV
%--------------------------------------------------------------------------
vol     = struct('tbv',NaN,'tiv',NaN);
vx      = sqrt(sum(Nii(1).mat(1:3,1:3).^2));
vol.tbv = prod(vx(1:3))*sum(Z(:,:,:,[1,2]),'all');
vol.tiv = prod(vx(1:3))*sum(Z(:,:,:,[1,2,3]),'all');

% Makes sure orientation matrix is correct (could have been modfied by call
% to correct_orientation).
M0 = spm_get_space(Nii(1).dat.fname);
spm_get_space(Nii(1).dat.fname, Mr\M0);
for k=1:K
    if isempty(res.c1{k}), continue; end
    spm_get_space(res.c1{k}, Mr\M0);
end

res.s = '';
if skullstrip
    % Produce skull-stripped CT scan (prefixed 's_')
    %----------------------------------------------------------------------
    % Copy image
    [pth,nam,ext] = fileparts(Nii(1).dat.fname);
    nfname        = fullfile(pth,['s_' nam ext]);
    copyfile(Nii(1).dat.fname,nfname);
    % Make mask and apply
    Nii_s = nifti(nfname);
    img   = single(Nii_s.dat());
    msk   = sum(Z(:,:,:,[1,2,3]),4) >= 0.5;
    img(~msk) = 0;
    % Modify copied image's data
    Nii_s.dat(:,:,:) = img;
    res.s = nfname;
    clear msk
end

% Compute template space segmentations
%--------------------------------------------------------------------------
if any(tc(:,2)) || any(tc(:,3))
    dmu     = sett.mu.d;
    Mmu     = sett.mu.Mmu;
    Mn      = dat(1).Mat;
    sd      = spm_mb_shape('samp_dens',Mmu,Mn);
    dir_res = sett.odir;
    onam    = dat(1).onam;
    % Get deformation
    Mat = Mmu\spm_dexpm(double(dat(1).q),sett.B)*Mn;
    psi = spm_mb_io('get_data',dat(1).psi);
    psi = MatDefMul(psi,inv(Mmu));
    psi = spm_mb_shape('compose',psi,spm_mb_shape('affine',dat(1).dm,Mat));
    % Normalise
    if any(tc(:,2)), res.wc  = cell(1,sum(tc(:,2))); end
    if any(tc(:,3)), res.mwc = cell(1,sum(tc(:,3))); end
    kwc  = 0;
    kmwc = 0;
    for k=1:K
        if tc(k,2) || tc(k,3)
            [img,cnt] = spm_mb_shape('push1',Z(:,:,:,k),psi,dmu,sd);
            if tc(k,2)
                % Write normalised segmentation
                kwc  = kwc + 1;
                fpth = fullfile(dir_res, sprintf('wc%.2d_%s.nii',k,onam));
                res.wc{kwc} = fpth;
                write_nii(fpth, img./(cnt + eps('single')),...
                         Mmu, sprintf('Norm. tissue (%d)',k), 'uint8');
            end
            if tc(k,3)
                % Write normalised modulated segmentation
                kmwc = kmwc + 1;
                fpth = fullfile(dir_res,sprintf('mwc%.2d_%s.nii',k,onam));
                res.mwc{kmwc} = fpth;
                img  = img*abs(det(Mn(1:3,1:3))/det(Mmu(1:3,1:3)));
                write_nii(fpth,img, Mmu, sprintf('Norm. mod. tissue (%d)',k), 'int16');
            end
            clear img cnt
        end
    end
end
clear Z
res.c  = cell(1,sum(tc(:,1)));
k1 = 1;
for k=1:K
    if ~tc(k,1)
        delete(res.c1{k});
    else
        res.c{k1} = res.c1{k};
        k1        = k1 + 1;
    end
end
res = rmfield(res,'c1');

% Save deformation?
res.y = '';
if def
    res.y = dat(1).psi.dat.fname;
else
    delete(dat(1).psi.dat.fname);    
end
delete(dat(1).v.dat.fname); % Delete velocity field
delete(p_res);              % Delete mb_fit_mb.mat

if mni
    % Reslice normalised segmentations to SPM's MNI space
    %----------------------------------------------------------------------
    load(pth_Mmni, 'Mmni');  % Generated using DARTEL's spm_klaff
    Nii_spm = nifti(fullfile(spm('Dir'),'tpm','TPM.nii'));
    Mspm = Nii_spm.mat;
    dmspm = Nii_spm.dat.dim(1:3);
    if ~isempty(res.wc)
        for k=1:numel(res.wc)
            if isempty(res.wc{k}), continue; end
            pth = res.wc{k};
            spm_get_space(pth, Mmni);
            reslice_img(pth, Mspm, dmspm, 'uint8');
        end
    end
    if ~isempty(res.mwc)
        for k=1:numel(res.mwc)
            if isempty(res.mwc{k}), continue; end
            pth = res.mwc{k};
            spm_get_space(pth, Mmni);
            reslice_img(pth, Mspm, dmspm, 'int16');
        end
    end    
end

return
%==========================================================================

%==========================================================================
function [Nii,Mr] = correct_orientation(Nii,odir)
f = nm_reorient(Nii.dat.fname,odir);
Mr = reset_origin(f);
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
M1   = [vx(1) 0      0         off(1)
           0      vx(2) 0      off(2)
           0      0      vx(3) off(3)
           0      0      0      1];
Mr = M1/M0;
spm_get_space(pth,Mr*M0); 
%==========================================================================

%==========================================================================
function npth = nm_reorient(pth,odir,vx,prefix,deg)
if nargin < 3, vx     = [];   end
if nargin < 4, prefix = 'r'; end
if nargin < 5, deg    = 1;    end

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

end % End loop over images
npth = VO.fname;
%==========================================================================

%==========================================================================
function reslice_img(pth, Mout, dout, typ, vx, deg, bc)
if nargin < 4, typ = 'float32'; end
if nargin < 5, vx = 1; end
if nargin < 6, deg = 1; end
if nargin < 7, bc  = 0; end

Nii = nifti(pth);
Min = Nii.mat;

vx = vx(1)*ones(1,3);
vx_out = sqrt(sum(Mout(1:3,1:3).^2));
D = diag([vx_out./vx 1]);
Mout = Mout/D;
dout  = floor(D(1:3,1:3)*dout')';
y = affine(dout,Min\Mout);

db    = [repmat(deg, [1 3]) repmat(bc, [1 3])];
dat   = single(Nii.dat());
mn    = min(dat(:));
mx    = max(dat(:));
dat   = spm_diffeo('bsplins',spm_diffeo('bsplinc',dat,db),y,db);
dat(~isfinite(dat)) = 0;
dat   = min(mx, max(mn, dat));

write_nii(pth,dat,Mout,Nii.descrip,typ)
%==========================================================================

%==========================================================================
function write_nii(pth,dat,M,descrip,typ)
if nargin<5, typ = 'float32'; end

if exist(pth,'file'), delete(pth); end

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

% CleanGWC()
function zn = clean_gwc(zn,ixt,level)
if nargin < 2 || isempty(ixt)
    % Default SPM12 template ordering
    ixt = struct('gm',1,'wm',2,'csf',3);
end
if nargin < 3, level = 2; end

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
function phi = MatDefMul(phi,M)
d   = size(phi);
phi = reshape(bsxfun(@plus,reshape(phi,[prod(d(1:3)),3])*M(1:3,1:3)',M(1:3,4)'),d);
%==========================================================================
