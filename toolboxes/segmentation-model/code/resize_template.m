function model = resize_template(model,do_resiz,vx,keep_neck,padding)
% Reslices input images to have the same size and orientation matrix.
% Also crops a little bit of the FOV.
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 2, do_resiz  = true; end
if nargin < 3, vx        = []; end % Empty -> does not change voxel size
if nargin < 4, keep_neck = false; end
if nargin < 5, padding   = 0;  end

if isempty(vx), vx = [NaN NaN NaN]; end

pth_template  = model.template.nii.dat.fname;
Vin           = spm_vol(pth_template);
[pth,nam,ext] = fileparts(pth_template);

% Get bounding-box of SPM template
PthTemplateSPM = fullfile(spm('dir'),'tpm','TPM.nii,');
V              = spm_vol(PthTemplateSPM);
bb             = world_bb(V(1));

% Pad BB
bb(1,:) = bb(1,:) - padding;
bb(2,:) = bb(2,:) + padding;

if keep_neck
    % Do not crop neck
    bb(1,3) = NaN;
end

if ~do_resiz
    bb = [NaN NaN NaN; NaN NaN NaN];
end

% Reslice 
prefix = 'res_'; % Prefix of resliced NIfTIs
deg    = 0;      % Interpolation degree
reslice_to_bb(Vin,bb,vx,prefix,deg);

% Update template file
movefile(fullfile(pth,[prefix nam ext]),pth_template);
model.template.nii = nifti(pth_template);
%==========================================================================

%==========================================================================
function Vo = reslice_to_bb(Vi,bb,vx,prefix,deg)

% reslice images one-by-one
Vo = spm_vol;
c  = 1;
for V=Vi'
    
    voxdim = vx;
    % default voxdim to current volume's voxdim, (from mat parameters)
    if any(isnan(voxdim))
        vprm = spm_imatrix(V.mat);
        vvoxdim = vprm(7:9);
        voxdim(isnan(voxdim)) = vvoxdim(isnan(voxdim));
    end
    voxdim = voxdim(:)';
    voxdim = round(voxdim,3);
    
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

% function model = resize_template(model,opt)
% % FORMAT model = resize_template(model,opt)
% %
% % Crop the template, removing unecessary background voxels
% %
% %__________________________________________________________________________
% % Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
% 
% % Read input
% pth_template = model.template.nii.dat.fname;
% keep_neck    = opt.template.keep_neck;
% 
% pth_spm_tpm = fullfile(spm('dir'),'tpm','TPM.nii,');
% 
% V_ref = spm_vol(pth_spm_tpm);
% bb    = world_bb(V_ref(1));
% 
% V_in = spm_vol(pth_template);
% dm   = V_in(1).dim;
% if keep_neck
%     bb(1,3) = NaN;
% end
% if dm(3) == 1
%     bb(:,3) = 0;
% end
% 
% resize_img(pth_template,opt.template.vs,bb);
%  
% % Update template in model
% model.template.nii = nifti(pth_template);
% 
% %==========================================================================
% function resize_img(imnames, Voxdim, BB, ismask)
% %  resize_img -- resample images to have specified voxel dims and BBox
% % resize_img(imnames, voxdim, bb, ismask)
% %
% % Output images will be prefixed with 'r', and will have voxel dimensions
% % equal to voxdim. Use NaNs to determine voxdims from transformation matrix
% % of input image(s).
% % If bb == nan(2,3), bounding box will include entire original image
% % Origin will move appropriately. Use world_bb to compute bounding box from
% % a different image.
% %
% % Pass ismask=true to re-round binary mask values (avoid
% % growing/shrinking masks due to linear interp)
% %
% % See also voxdim, world_bb
% 
% % Based on John Ashburner's reorient.m
% % http://www.sph.umich.edu/~nichols/JohnsGems.html#Gem7
% % http://www.sph.umich.edu/~nichols/JohnsGems5.html#Gem2
% % Adapted by Ged Ridgway -- email bugs to drc.spm@gmail.com
% 
% % This version doesn't check spm_flip_analyze_images -- the handedness of
% % the output image and matrix should match those of the input.
% 
% % Check spm version:
% if exist('spm_select','file') % should be true for spm5
%     spm5 = 1;
% elseif exist('spm_get','file') % should be true for spm2
%     spm5 = 0;
% else
%     error('Can''t find spm_get or spm_select; please add SPM to path')
% end
% 
% spm_defaults;
% 
% % prompt for missing arguments
% if ( ~exist('imnames','var') || isempty(char(imnames)) )
%     if spm5
%         imnames = spm_select(inf, 'image', 'Choose images to resize');
%     else
%         imnames = spm_get(inf, 'img', 'Choose images to resize');
%     end
% end
% % check if inter fig already open, don't close later if so...
% Fint = spm_figure('FindWin', 'Interactive'); Fnew = [];
% if ( ~exist('Voxdim', 'var') || isempty(Voxdim) )
%     Fnew = spm_figure('GetWin', 'Interactive');
%     Voxdim = spm_input('Vox Dims (NaN for "as input")? ',...
%         '+1', 'e', '[nan nan nan]', 3);
% end
% if ( ~exist('BB', 'var') || isempty(BB) )
%     Fnew = spm_figure('GetWin', 'Interactive');
%     BB = spm_input('Bound Box (NaN => original)? ',...
%         '+1', 'e', '[nan nan nan; nan nan nan]', [2 3]);
% end
% if ~exist('ismask', 'var')
%     ismask = false;
% end
% if isempty(ismask)
%     ismask = false;
% end
% 
% % reslice images one-by-one
% vols = spm_vol(imnames);
% for V=vols'
%     % (copy to allow defaulting of NaNs differently for each volume)
%     voxdim = Voxdim;
%     bb = BB;
%     % default voxdim to current volume's voxdim, (from mat parameters)
%     if any(isnan(voxdim))
%         vprm = spm_imatrix(V.mat);
%         vvoxdim = vprm(7:9);
%         voxdim(isnan(voxdim)) = vvoxdim(isnan(voxdim));
%     end
%     voxdim = voxdim(:)';
% 
%     mn = bb(1,:);
%     mx = bb(2,:);
%     % default BB to current volume's
%     if any(isnan(bb(:)))
%         vbb = world_bb(V);
%         vmn = vbb(1,:);
%         vmx = vbb(2,:);
%         mn(isnan(mn)) = vmn(isnan(mn));
%         mx(isnan(mx)) = vmx(isnan(mx));
%     end
% 
%     % voxel [1 1 1] of output should map to BB mn
%     % (the combination of matrices below first maps [1 1 1] to [0 0 0])
%     mat = spm_matrix([mn 0 0 0 voxdim])*spm_matrix([-1 -1 -1]);
%     % voxel-coords of BB mx gives number of voxels required
%     % (round up if more than a tenth of a voxel over)
%     imgdim = ceil(mat \ [mx 1]' - 0.1)';
% 
%     % output image
%     VO            = V;
%     [pth,nam,ext] = fileparts(V.fname);
%     VO.fname      = fullfile(pth,['r' nam ext]);
%     VO.dim(1:3)   = imgdim(1:3);
%     VO.mat        = mat;
%     VO = spm_create_vol(VO);
%     spm_progress_bar('Init',imgdim(3),'reslicing...','planes completed');
%     for i = 1:imgdim(3)
%         M = inv(spm_matrix([0 0 -i])*inv(VO.mat)*V.mat);
%         img = spm_slice_vol(V, M, imgdim(1:2), 1); % (linear interp)
%         if ismask
%             img = round(img);
%         end
%         spm_write_plane(VO, img, i);
%         spm_progress_bar('Set', i)
%     end
%     spm_progress_bar('Clear');
% end
% % call spm_close_vol if spm2
% if ~spm5
%     spm_close_vol(VO);
% end
% if (isempty(Fint) && ~isempty(Fnew))
%     % interactive figure was opened by this script, so close it again.
%     close(Fnew);
% end
% disp('Done.')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function bb = world_bb(V)
% %  world-bb -- get bounding box in world (mm) coordinates
% 
% d = V.dim(1:3);
% % corners in voxel-space
% c = [ 1    1    1    1
%     1    1    d(3) 1
%     1    d(2) 1    1
%     1    d(2) d(3) 1
%     d(1) 1    1    1
%     d(1) 1    d(3) 1
%     d(1) d(2) 1    1
%     d(1) d(2) d(3) 1 ]';
% % corners in world-space
% tc = V.mat(1:3,1:4)*c;
% 
% % bounding box (world) min and max
% mn = min(tc,[],2)';
% mx = max(tc,[],2)';
% bb = [mn; mx];