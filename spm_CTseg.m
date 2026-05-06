function [res,vol] = spm_CTseg(in, odir, tc, def, correct_header, skullstrip, vox, v_settings, tol, mu, hemisphere, bb)
% A CT segmentation+spatial normalisation routine for SPM12.
% FORMAT [res,vol] = spm_CTseg(in, odir, tc, def, correct_header, skullstrip, vox, v_settings, tol, mu, hemisphere, bb)
%
% This algorithm produces native|warped|modulated space segmentations of:
%     1. Gray matter (GM)        [or GM hemisphere 1 & 2 if hemisphere=true]
%     2. White matter (WM)       [or WM hemisphere 1 & 2 if hemisphere=true]
%     3. Cerebrospinal fluid (CSF)
%     4. Bone (BONE)
%     5. Soft tissue (ST)
%     6. Background (BG)
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
% tc (logical(K, 3)): Matrix where native, warped and warped modulated are
%                     indexed by columns and tissue classes are indexed by rows 
%                     (in the above order).             
%
% def (logical): Write deformations? Defaults to true.
%
% correct_header (logical): Correct messed up CT header, defaults to false. 
%
% skullstrip (logical): Write skull-stripped CT scan to disk, prefixed 
%                       'ss_'. Defaults to false.
%
% vox (double): Template space voxel size, defaults to voxel size of
%               template.
%
% v_settings (int|int(1,5)): Spatial regularisation settings. See Multi-
%                            Brain toolbox. If singleton, acts as a
%                            multiplication factor on the default.
%
% tol (double): Stopping tolerance. Defaults to 0.5*0.001. Larger = faster
%               and less accurate.
%
% mu (char): Path to tissue template, or a shorthand atlas name.
%            Shorthand names:
%              'spm15' - SPM-aligned, 1.5 mm (default)
%              'spm10' - SPM-aligned, 1.0 mm
%              'ctseg' - original groupwise optimal atlas (legacy;
%                        requires spm_CTseg_warp for MNI output)
%            If a shorthand is given, the atlas is auto-downloaded to
%            models/ on first use. If empty, uses 'spm15'. A full file
%            path to a custom atlas is also accepted.
%
% hemisphere (logical): Separate GM and WM into left/right hemisphere
%                       segmentations (8 classes instead of 6). Defaults
%                       to false.
%
% bb (2x3 double | char): Bounding box (in mm) for template-space outputs
%                  (wc*, mwc*). Use to crop outputs to a specific FOV.
%                  Shorthand strings:
%                    'spm'  - SPM TPM bounding box, i.e.
%                             spm_get_bbox(fullfile(spm('Dir'),'tpm','TPM.nii'),'old')
%                             (default; matches SPM TPM FOV for all atlases)
%                    'full' - full input atlas FOV (no cropping)
%                  A numeric 2x3 matrix is also accepted.
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
if size(tc,2) == 1
    tc = repmat(tc, 1, 3);
end
if nargin < 4, def            = true; end
if nargin < 5, correct_header = false; end
if nargin < 6, skullstrip     = false; end
if nargin < 7, vox            = NaN; end
if nargin < 8 || isempty(v_settings)
    v_settings = [0.0001 0 0.4 0.1 0.4] * 3;
elseif numel(v_settings) == 1
    v_settings = [0.0001 0 0.4 0.1 0.4] .* v_settings;
end
if nargin < 9  || isempty(tol), tol = 0.001; end
if nargin < 10 || isempty(mu),  mu  = ''; end
if nargin < 11 || isempty(hemisphere), hemisphere = false; end
if nargin < 12 || isempty(bb), bb = 'spm'; end
if ischar(bb) || isstring(bb)
    switch lower(char(bb))
        case 'spm'
            bb = spm_get_bbox(fullfile(spm('Dir'),'tpm','TPM.nii'),'old');
        case 'full'
            bb = NaN(2,3);
        otherwise
            error('Unknown bb shorthand ''%s''. Use ''spm'', ''full'', or a 2x3 numeric matrix.', char(bb));
    end
end

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
% Add bundled Multi-Brain toolbox, removing any competing copies.
% In a deployed (compiled) standalone, paths are baked in at compile time
% and MEX auto-compilation is unavailable, so skip the path/compile logic.
dir_mb = fullfile(fileparts(mfilename('fullpath')), 'mb');
if ~isdeployed
    if ~exist(fullfile(dir_mb, 'spm_mb_fit.m'), 'file')
        error('Multi-Brain toolbox not found in CTseg/mb/. Run: git submodule update --init');
    end
    mb_on_path = which('spm_mb_fit', '-all');
    for i = 1:numel(mb_on_path)
        mb_dir_i = fileparts(mb_on_path{i});
        if ~strcmp(mb_dir_i, dir_mb)
            rmpath(mb_dir_i);
        end
    end
    addpath(dir_mb);
    if ~exist(fullfile(dir_mb, ['spm_gmmlib.' mexext]), 'file')
        % Try to compile automatically
        fprintf('Compiling Multi-Brain GMM library... ')
        cwd = pwd;
        cd(dir_mb);
        try
            mex -O -largeArrayDims spm_gmmlib.c gmmlib.c
            fprintf('done.\n')
        catch ME
            cd(cwd);
            error(['Failed to compile spm_gmmlib: %s\n' ...
                   'Go to CTseg/mb/ and compile manually (see README).'], ME.message);
        end
        cd(cwd);
    end
end

% Get model files
%--------------------------------------------------------------------------
dir_ctseg  = fileparts(mfilename('fullpath'));
dir_models = fullfile(dir_ctseg, 'models');

% Warn if old file layout detected (files in CTseg root)
old_files = {'mu_CTseg.nii', 'prior_CTseg.mat', 'Mmni.mat'};
old_found = {};
for i = 1:numel(old_files)
    if exist(fullfile(dir_ctseg, old_files{i}), 'file') == 2
        old_found{end+1} = old_files{i}; %#ok<AGROW>
    end
end
if ~isempty(old_found)
    warning('CTseg:oldLayout', ...
        ['Found model files in the CTseg root directory (old layout).\n' ...
         'Model files have moved to the ''models'' subdirectory.\n' ...
         'You can safely delete from %s: %s'], ...
         dir_ctseg, strjoin(old_found, ', '));
end

% Verify intensity prior
pth_int = fullfile(dir_models, 'prior_CTseg.mat');
if ~(exist(pth_int, 'file') == 2)
    error('Intensity prior file (models/prior_CTseg.mat) could not be found!')
end

% Resolve atlas path from mu parameter
%--------------------------------------------------------------------------
registry = get_atlas_registry();
if isempty(mu)
    % Default atlas (SPM-aligned 1.5 mm)
    pth_mu = download_atlas('spm15', dir_models);
    use_default_mu = false;
elseif isfield(registry, mu)
    % Shorthand name
    use_default_mu = strcmp(mu, 'ctseg');
    pth_mu = download_atlas(mu, dir_models);
else
    % Treat as file path
    pth_mu = mu;
    use_default_mu = false;
end
if ~(exist(pth_mu, 'file') == 2)
    error('Atlas file (%s) could not be found!', pth_mu)
end
if use_default_mu
    pth_Mmni = fullfile(dir_models, 'Mmni.mat');
    if ~(exist(pth_Mmni, 'file') == 2)
        error('MNI affine (models/Mmni.mat) could not be found!')
    end
end

% Get nifti
%--------------------------------------------------------------------------
Nii = nifti(in);

% Output directory
%--------------------------------------------------------------------------
if isempty(odir)
    odir = fileparts(Nii.dat.fname);
    odir = spm_file(odir,'cpath'); % Get absolute path
elseif ~(exist(odir, 'dir') == 7)
    mkdir(odir);
end

fprintf('\n--- CTseg started ---\n');
fprintf('Atlas:  %s\n', pth_mu);
fprintf('Output: %s\n', odir);

% Correct orientation matrix
%--------------------------------------------------------------------------
Mc = eye(4);
oNii = Nii;
if correct_header
    [Nii,Mc] = correct_orientation(Nii, odir);
end
% Get number of tissue classes from template
Nii_mu = nifti(pth_mu);
K      = Nii_mu.dat.dim(4) + 1;

% Hemisphere segmentation: split GM and WM into left/right
%--------------------------------------------------------------------------
if hemisphere && K == 6
    % Create temporary modified atlas and prior in output directory
    pth_int_orig = pth_int;
    pth_mu  = fullfile(odir, 'temp_mu_hemisphere.nii');
    pth_int = fullfile(odir, 'temp_prior_hemisphere.mat');
    % Duplicate GM and WM channels
    mu_dat = single(Nii_mu.dat());
    mu_dat = mu_dat(:,:,:,[1, 1, 2, 2, 3, 4, 5]);
    % Separate hemispheres by zeroing out the opposite half
    ix0                       = floor(0.5*size(mu_dat,1));
    ix1                       = ceil(0.5*size(mu_dat,1));
    min_mu                    = min(mu_dat(:));
    mu_dat(1:ix0,:,:,[1 3])   = min_mu;
    mu_dat(ix1:end,:,:,[2 4]) = min_mu;
    % Write modified template
    spm_CTseg_util('write_nii',pth_mu,mu_dat,Nii_mu.mat,'Hemisphere template','float32');
    % Load and modify intensity prior
    pr    = load(pth_int_orig);
    mg_ix = pr.mg_ix;
    pr    = pr.pr;
    mg_ix = [1 2 3 4 mg_ix(3:end) + 2];
    ix    = [1 1 2 2 3:size(pr{1},2)];
    pr{1} = pr{1}(:,ix);
    pr{2} = pr{2}(:,ix);
    pr{3} = pr{3}(:,:,ix);
    pr{4} = pr{4}(:,ix);
    pr{5} = pr{5}(:,ix);
    save(pth_int,'pr','mg_ix');
    % Update number of classes
    K = size(mu_dat,4) + 1;
    clear mu_dat
end

% Tissue class indices
if hemisphere
    ix_gm  = [1 2];
    ix_wm  = [3 4];
    ix_csf = [5];
else
    ix_gm  = [1];
    ix_wm  = [2];
    ix_csf = [3];
end
ix_gwc = [ix_gm, ix_wm, ix_csf];

if size(tc,1) == 1
    tc = repmat(tc, K, 1);
end
% For keeping modulated, if requested
tc0 = tc;
if nargout > 1
    tc(ix_gwc,3) = true;
end

% Ensure SPM can discover the bundled MB toolbox via a directory junction/symlink.
% This is needed because SPM's batch system only discovers toolboxes in spm/toolbox/.
% In a deployed standalone, mb is already bundled as a top-level toolbox at
% compile time (see scripts/build_standalone.m) and batch config is baked in,
% so skip the runtime symlink/addpath dance.
%--------------------------------------------------------------------------
if ~isdeployed
    spm_mb_link = fullfile(spm('Dir'), 'toolbox', 'mb');
    needs_link  = false;
    if ~exist(spm_mb_link, 'dir')
        needs_link = true;
    elseif ~exist(fullfile(spm_mb_link, 'spm_mb_output.m'), 'file')
        % Link exists but is broken or points to wrong directory
        if ispc
            system(sprintf('rmdir "%s"', spm_mb_link));
        else
            system(sprintf('rm -f "%s"', spm_mb_link));
        end
        needs_link = true;
    end
    if needs_link
        if ispc
            [st,msg] = system(sprintf('mklink /J "%s" "%s"', spm_mb_link, dir_mb));
        else
            [st,msg] = system(sprintf('ln -s "%s" "%s"', dir_mb, spm_mb_link));
        end
        if st ~= 0
            error(['Could not create toolbox link: %s\n' ...
                   'Manually create a link/junction from %s to %s'], msg, spm_mb_link, dir_mb);
        end
    end
    addpath(spm_mb_link);
    spm_jobman('initcfg');
    % Re-enforce bundled MB after batch init (SPM may re-add competing copies)
    mb_on_path = which('spm_mb_fit', '-all');
    for i = 1:numel(mb_on_path)
        mb_dir_i = fileparts(mb_on_path{i});
        if ~strcmp(mb_dir_i, dir_mb)
            rmpath(mb_dir_i);
        end
    end
    addpath(dir_mb);
end

% Run MB
%--------------------------------------------------------------------------
% algorithm settings
run              = struct;
run.mu.exist     = {pth_mu};
run.onam         = 'CTseg';
run.odir         = {odir};
run.v_settings   = v_settings;
run.tol          = tol;
run.aff          = 'Aff(3)';
run.del_settings = 1;
% image
run.gmm.pr.file          = {pth_int};
run.gmm.pr.hyperpriors   = [];
run.gmm.chan.images      = {Nii(1).dat.fname};
run.gmm.chan.modality    = 2;
run.gmm.chan.inu.inu_reg = 1e7;
% output settings
out           = struct;
out.result    = {fullfile(odir,['mb_fit_' run.onam '.mat'])};
out.c         = 1:K;
out.wc        = find(tc(:,2))';
out.mwc       = find(tc(:,3))';
out.vox       = vox;
out.mrf       = 1;
out.bb        = bb;
out.proc_zn   = {@(x) clean_gwc(x, struct('gm',ix_gm,'wm',ix_wm,'csf',ix_csf))};

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

vol = struct('tbv',NaN,'tiv',NaN);
if nargout > 1
    % Compute TBV and TIV, note that these are computed using the
    % modulated GM, WM and CSF in template space, with the field-of-view
    % of the SPM12 atlas as the CTseg atlas has a larger FOV.
    % ------------
    if use_default_mu
        % zero voxels outside of SPM12 atlas field-of-view
        % (only needed when template FOV differs from SPM atlas)
        pth_spm = nifti(fullfile(spm('Dir'),'tpm','TPM.nii'));
        for k=ix_gwc
            spm_CTseg_util('mask_outside_fov', pth_spm, res.mwc{k});
        end
    end
    % Compute TBV and TIV from modulated template space segmentations
    vol = struct('tbv',0,'tiv',0);
    for k=ix_gwc
        Nii_mwc = nifti(res.mwc{k});
        sm_dat = sum(Nii_mwc.dat(:));
        if ismember(k, [ix_gm ix_wm])
            vol.tbv = vol.tbv + sm_dat;
        end
        vol.tiv = vol.tiv + sm_dat;
    end
    vx = sqrt(sum(Nii_mwc(1).mat(1:3,1:3).^2));
    vol.tbv = prod(vx(1:3))*vol.tbv / 1000;  % mm^3 -> ml
    vol.tiv = prod(vx(1:3))*vol.tiv / 1000;  % mm^3 -> ml
    for k=ix_gwc
        if ~tc0(k, 3)
            spm_unlink(res.mwc{k});
            res.mwc{k} = [];
        end
    end
    tc = tc0;
end

% % Reslice template space segmentations to MNI space
% % (only needed when using default template, which is not in SPM space)
% if use_default_mu
%     reslice2mni(res,pth_Mmni,Mmu);
% end

if correct_header
    % Reslice corrected native space segmentations to original native space.
    M1  = spm_get_space(Nii(1).dat.fname);  % get corrected orientation matrix
    spm_unlink(Nii(1).dat.fname);           % delete corrected image
    Nii = oNii;                             % reset to original input image
    M0  = spm_get_space(Nii(1).dat.fname);  % get corrected orientation matrix
    % new field-of-view
    M = Mc\M1\M0;
    y = spm_CTseg_util('affine', Nii.dat.dim, M);     
    % reslice segmentations
    for k=1:K
        if isempty(res.c{k}), continue; end                        
        Nii_c = nifti(res.c{k});        
        rc    = spm_diffeo('bsplins',single(Nii_c.dat()),y,[1 1 1  0 0 0]);
        rc    = max(rc, 0);  % clamp negative interpolation artifacts
        spm_CTseg_util('write_nii',res.c{k},rc,M0,sprintf('Tissue (%d)',k), 'uint8')        
    end
end

% Mask out-of-FOV voxels in native space segmentations
% (set to zero where native image extends beyond atlas coverage)
for k=1:K-1
    if ~isempty(res.c{k})
        spm_CTseg_util('mask_outside_fov', pth_mu, res.c{k});
    end
end

res.s = '';
if skullstrip
    % Produce skull-stripped CT scan (prefixed 'ss_')
    %----------------------------------------------------------------------
    % Get native-space responsibilities    
    Z = [];
    for k=1:K
        Nii_c = nifti(res.c{k});
        Z     = cat(4, Z, single(Nii_c.dat()));
    end
    Z = bsxfun(@rdivide, Z, sum(Z,4) + eps('single'));  % renormalise (resps could have been resliced)
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
    clear msk Z
end

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
    if correct_header
        % adjust affine of deformation
        M0 = spm_get_space(dat(1).psi.dat.fname);
        spm_get_space(dat(1).psi.dat.fname, Mc\M0);
    end
    res.y = dat(1).psi.dat.fname;
else
    spm_unlink(dat(1).psi.dat.fname); % Delete deformation
end
spm_unlink(dat(1).v.dat.fname); % Delete velocity field
spm_unlink(fullfile(run.odir{1},['mb_fit_' run.onam '.mat'])); % Delete MB fit results

% Clean up temporary hemisphere files
if hemisphere
    spm_unlink(pth_mu);
    spm_unlink(pth_int);
end

% Print summary
%--------------------------------------------------------------------------
fprintf('\n--- CTseg finished ---\n');
fprintf('Output directory: %s\n', run.odir{1});
for k=1:numel(res.c)
    if ~isempty(res.c{k}),   fprintf('  Native:           %s\n', res.c{k}); end
end
for k=1:numel(res.wc)
    if ~isempty(res.wc{k}),  fprintf('  Warped:           %s\n', res.wc{k}); end
end
for k=1:numel(res.mwc)
    if ~isempty(res.mwc{k}), fprintf('  Warped modulated: %s\n', res.mwc{k}); end
end
if ~isempty(res.s), fprintf('  Skull-stripped:    %s\n', res.s); end
if ~isempty(res.y), fprintf('  Deformation:      %s\n', res.y); end
if nargout > 1
    fprintf('  TBV: %.1f ml, TIV: %.1f ml\n', vol.tbv, vol.tiv);
end
fprintf('----------------------\n');
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
if nargin < 5, deg    = 0;    end

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
    % Use NaN for out-of-bounds voxels to avoid zero-padding artifacts
    img = spm_slice_vol(V,M,dim(1:2),deg);

    % Write this slice to output image
    spm_write_plane(VO,img,i);
end % End loop over output slices

npth = VO.fname;
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
Md = Mmni\Ms;
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
y = spm_CTseg_util('affine',dout,M);
% Reslice
Nii   = nifti(pth);
dat   = spm_diffeo('bsplins',single(Nii.dat()),y,[1 1 1  0 0 0]);
spm_CTseg_util('write_nii',pth,dat,Mout,Nii.descrip,typ);
%==========================================================================

%==========================================================================
function zn = clean_gwc(zn,ixt,level)
if nargin < 2 || isempty(ixt)
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
    end
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
function registry = get_atlas_registry()
% Returns a struct mapping atlas shorthands to filenames and download URLs.
registry.spm15   = struct('file','mu_CTseg_spm15.nii', 'url','https://github.com/WCHN/CTseg/releases/download/v1.0/mu_CTseg_spm15.nii.gz');
registry.spm10   = struct('file','mu_CTseg_spm10.nii', 'url','https://github.com/WCHN/CTseg/releases/download/v1.0/mu_CTseg_spm10.nii.gz');
registry.ctseg   = struct('file','mu_CTseg.nii',       'url','https://github.com/WCHN/CTseg/releases/download/v1.0/mu_CTseg.nii.gz');
%==========================================================================

%==========================================================================
function pth = download_atlas(name, dir_models)
% Download atlas by shorthand name if not already present.
registry = get_atlas_registry();
if ~isfield(registry, name)
    error('Unknown atlas ''%s''. Valid names: %s', ...
        name, strjoin(fieldnames(registry), ', '));
end
entry = registry.(name);
pth   = fullfile(dir_models, entry.file);
if exist(pth, 'file') == 2
    return
end
pth_gz = fullfile(dir_models, [entry.file '.gz']);
fprintf('Downloading atlas ''%s'' (first use only)... ', name)
websave(pth_gz, entry.url);
fprintf('done.\n')
fprintf('Extracting... ')
gunzip(pth_gz, dir_models);
delete(pth_gz);
fprintf('done.\n')
%==========================================================================
