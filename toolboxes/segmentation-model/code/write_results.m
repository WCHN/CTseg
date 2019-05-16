function res = write_results(dat,model,opt)

res = struct;

%--------------------------------------------------------------------------
% Image data, etc
%--------------------------------------------------------------------------

% samp = 0 -> to here work with the original size image (not a sub-sampled version)
[obs,dm_s,mat_s,vs_s,scl] = get_obs(dat,'samp',0); 

[~,~,~,C,nam,~,obs_fnames] = obs_info(dat);
dm_a                       = model.template.nii.dat.dim;
mat_a                      = model.template.nii.mat;
vs_a                       = sqrt(sum(mat_a(1:3,1:3).^2));
K                          = dm_a(4);
modality                   = dat.modality{1}.name; 
ff                         = get_ff(vs_s);   

%--------------------------------------------------------------------------
% Get deformation
%--------------------------------------------------------------------------

% Get velocities
if isnumeric(dat.reg.v)
    v = dat.reg.v;                       
else
    v = single(dat.reg.v.dat(:,:,:,:));   
end

% Registration parameters
subsmp = get_subsampling_grid(dm_s,vs_s,opt.seg.samp);
prm_v  = [subsmp.sk.*vs_s ff*opt.reg.rparam];   

% Green's function
if opt.reg.int_args > 1, Greens = spm_shoot_greens('kernel',subsmp.dm(1:3),prm_v);
else,                    Greens = [];
end

if 0
    show_bf_and_ivel(obs,size(v),v);
end

% Deformation from velocities
y = make_deformation(v,prm_v,opt.reg.int_args,Greens);  

if opt.seg.samp > 0
    % Resize deformations, from sub-sampled size to original subject's image 
    % size     
    y = resize_def(y,dm_s,subsmp);
end

if 0
    show_bf_and_ivel(obs,size(y),y);
end

%--------------------------------------------------------------------------
% Warp template to subject
%--------------------------------------------------------------------------

E        = spm_dexpm(dat.reg.r,opt.reg.B); 
Affine   = model.template.nii.mat\E*mat_s;               
Template = warp_template(model,y,Affine);
y        = spm_warps('transform',Affine,y);
if dm_s(3) == 1
    y(:,:,:,3) = 1;
end

%--------------------------------------------------------------------------
% Write deformation
%--------------------------------------------------------------------------

if opt.write.df
    Nii         = nifti;
    Nii.dat     = file_array(fullfile(dat.dir.def,['def', nam{1}, '.nii']),...
                             size(y),...
                             [spm_type('float32') spm_platform('bigend')],...
                             0,1,0);
    Nii.mat     = mat_s;
    Nii.mat0    = mat_s;
    Nii.descrip = 'Initial velocity';
    create(Nii);

    Nii.dat(:,:,:,:) = y;
end

%--------------------------------------------------------------------------
% Get bias field
%--------------------------------------------------------------------------

do_bf = opt.bf.do && strcmpi(modality,'MRI');
if do_bf
    
    if opt.seg.samp > 0
        % Adjust bias field if images were subsampled
        %------------------------------------------------------------------
        
        [x1,y1] = ndgrid(1:dm_s(1),1:dm_s(2),1);
        z1      = 1:dm_s(3);

        for c=1:C
            d3                = [size(dat.bf.chan(c).T) 1];
            dat.bf.chan(c).B3 = spm_dctmtx(dm_s(3),d3(3),z1);
            dat.bf.chan(c).B2 = spm_dctmtx(dm_s(2),d3(2),y1(1,:)');
            dat.bf.chan(c).B1 = spm_dctmtx(dm_s(1),d3(1),x1(:,1));
        end
        clear x1 y1 z1
    end
    
    bf = get_bf(dat.bf.chan,dm_s);
else     
    bf = 1;
end

if 0
    show_bf_and_ivel(obs,dm_s,bf);
end

%--------------------------------------------------------------------------
% Get responsibilities, using the template to fill in not observed values
%--------------------------------------------------------------------------

miss   = get_par('missing_struct',obs);
labels = get_labels(dat,opt);
Z      = get_final_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt);
clear labels

if 0
    show_seg(obs,Template,dat.gmm.prop,Z,dm_s,modality);
end
clear Template

%--------------------------------------------------------------------------
% Write non-preprocessed responsibilities to disk
%--------------------------------------------------------------------------

for k=find(opt.write.tc(:,1)' == true)
    fname = fullfile(dat.dir.seg_orig,['oc' num2str(k) nam{1} '.nii']);
    spm_misc('create_nii',fname,Z(:,:,:,k),mat_s,[spm_type('float32') spm_platform('bigend')],'Segmentation');
end

%--------------------------------------------------------------------------
% Write bias field
%--------------------------------------------------------------------------
    
if opt.write.bf(2) && do_bf
    for c=1:C
        Nii         = nifti;
        Nii.dat     = file_array(fullfile(dat.dir.bf,['bf', nam{c}, '.nii']),...
                                 dm_s,...
                                 [spm_type('float32') spm_platform('bigend')],...
                                 0,1,0);
        Nii.mat     = mat_s;
        Nii.mat0    = mat_s;
        Nii.descrip = 'Estimated Bias Field';
        create(Nii);

        Nii.dat(:,:,:) = reshape(bf(:,c),dm_s);
    end
end

%--------------------------------------------------------------------------
% Write images (bias-field corrected)
%--------------------------------------------------------------------------

% Multiply observations with bias fields
obs = get_obs(dat,'samp',0,'mskonlynan',true); % We want to write the non-masked images
obs = bf.*obs;
clear bf

if opt.write.bf(1)
    for c=1:C  
        Nii         = nifti;
        Nii.dat     = file_array(fullfile(dat.dir.img,['i', nam{c}, '.nii']),...
                                 dm_s,...
                                 [spm_type('float32') spm_platform('bigend')],...
                                 0,1,0);
        Nii.mat     = mat_s;
        Nii.mat0    = mat_s;
        Nii.descrip = 'Bias Field Corrected Image';
        create(Nii);

        Nii.dat(:,:,:) = reshape(obs(:,c),dm_s);
    end
end

%--------------------------------------------------------------------------
% Write bias-field corrected images, warped to MNI space
%--------------------------------------------------------------------------

if opt.write.bf(3)
    for c=1:C    
        Nii         = nifti;
        Nii.dat     = file_array(fullfile(dat.dir.img,['wi', nam{c}, '.nii']),...
                                 dm_a(1:3),...
                                 [spm_type('float32') spm_platform('bigend')],...
                                 0,1,0);
        Nii.mat     = mat_a;
        Nii.mat0    = mat_a;
        Nii.descrip = 'Bias Field Corrected Image in MNI space';
        create(Nii);

        [x,c]           = spm_diffeo('push',reshape(obs(:,c),dm_s),y,dm_a(1:3));   
        x               = x./c;
        x(~isfinite(x)) = 0;
        if strcmpi(modality,'MRI')
            x(x < 0)    = 0;
        end
        
        Nii.dat(:,:,:) = single(x);        
    end    
end
clear obs x iy c

%--------------------------------------------------------------------------
% Some cleaning up of responsibilities
%--------------------------------------------------------------------------

if opt.clean.brain && dm_a(3) > 1
    % Clean responsibilities using template, deformation and morphological operations
    Z = clean_brain(Z,model,opt,y);                
end

if opt.clean.mrf.strength > 0
    % Ad-hoc MRF clean-up of responsibilities
    Z = mrf_post(Z,vs_s,opt);
end

if 0
    % Some verbose
    figure(672)
    nr = floor(sqrt(K));
    nc = ceil(K/nr);      
    for k=1:K   
        subplot(nr,nc,k)
        imagesc3d(Z(:,:,:,k)); axis off; drawnow
    end
end
    
if isempty(opt.lesion.hemi) && (opt.clean.les.bwlabeln || opt.clean.les.cnn_mrf.do)
    % Try to extract binary lesion representation from lesion class of responsibilities
    
    % Get lesion class
    Z_les = Z(:,:,:,opt.clean.les.class); 
    Z_les(~msk_les) = 0;
    clear msk_les
    
    if opt.clean.les.cnn_mrf.do                        
        % CNN-MRF clean-up
        x    = cell(1,2);
        x{1} = single(Z_les > opt.clean.les.val);
        x{2} = 1 - x{1};

        im_clean = spm_cnn_mrf('predict',opt.clean.les.cnn_mrf.pth_net,x,'speak',0,'w_on',0.9,'nit',40);
        im_clean = im_clean{1};
        clear x
    end
            
    if opt.clean.les.cnn_mrf.do
        bin_les = im_clean > opt.clean.les.val;
    else
        bin_les = Z_les    > opt.clean.les.val;
    end
    clear im_clean Z_les
    
    % Get largest connected component
    L          = bwlabeln(bin_les);
    nL         = unique(L);
    vls        = zeros([1 numel(nL)]);
    for i=1:numel(nL)
        vls(i) = sum(sum(sum(L == nL(i))));
    end
    [~,ix]     = sort(vls);

    if numel(ix) > 1
        bin_les = L == (ix(end - 1) - 1);
        bin_les = imfill(bin_les,'holes');
    end    
    bin_les = single(bin_les); 
        
    if opt.write.les(1)
        % Write to disk
        fname = fullfile(dat.dir.seg,['cles' nam{1} '.nii']);
        spm_misc('create_nii',fname,bin_les,mat_s,[spm_type('float32') spm_platform('bigend')],'Lesion (native)');

        res.cles = nifti(fname);
    end
    
%     if 0 && isfield(dat,'label')
%         gt  = dat.label{1}.nii.dat(:,:,:) > 0;
%         prd = logical(les_msk);        
%         fprintf('ds=%0.3f\n',dice(prd,gt));
%     end
end

%--------------------------------------------------------------------------
% Write DARTEL imports
%--------------------------------------------------------------------------

if any(opt.write.tc(:,5) == true) && dm_s(3) > 1 
    write_dartel(Z,y,mat_s,dm_s,dm_a,mat_a,vs_a,vs_s,dat.dir.seg,nam{1},opt);
end

%--------------------------------------------------------------------------
% Write final segmentation to disk
%--------------------------------------------------------------------------

% Native space (c)
for k=find(opt.write.tc(:,2)' == true)
    fname = fullfile(dat.dir.seg,['c' num2str(k) nam{1} '.nii']);
    spm_misc('create_nii',fname,Z(:,:,:,k),mat_s,[spm_type('float32') spm_platform('bigend')],'Segmentation (native)');
    
    res.c{k} = nifti(fname);
end

if any(opt.write.tc(:,3) == true)
    % Warped to template space (wc)   
    C = zeros(dm_a,'single');
    for k=1:size(Z,4)
        [c,w]      = spm_diffeo('push',Z(:,:,:,k),y,dm_a(1:3));    
        C(:,:,:,k) = spm_field(w,c,[vs_a  1e-6 1e-4 0  3 2]);
        clear w c
    end

    C = max(C,eps);
    s = sum(C,4);
    for k=find(opt.write.tc(:,3)' == true)   
        fname = fullfile(dat.dir.seg,['wc' num2str(k) nam{1} '.nii']);
        spm_misc('create_nii',fname,C(:,:,:,k)./s,mat_a,[spm_type('float32') spm_platform('bigend')],'Segmentation (template)');

        res.wc{k} = nifti(fname);
    end
    clear s C
end

if any(opt.write.tc(:,4) == true)
    % Modulated, warped to template space
    
    for k=find(opt.write.tc(:,4)' == true)   
        c = spm_diffeo('push',Z(:,:,:,k),y,dm_a(1:3));
        c = c*abs(det(mat_s(1:3,1:3))/det(mat_a(1:3,1:3)));
        
        fname = fullfile(dat.dir.seg,['mwc' num2str(k) nam{1} '.nii']);
        spm_misc('create_nii',fname,c,mat_a,[spm_type('float32') spm_platform('bigend')],'Segmentation (modulate template)');

        res.mwc{k} = nifti(fname);
    end
    clear c
end

if ((opt.clean.les.cnn_mrf.do || opt.clean.les.bwlabeln) && opt.write.les(2))
    [c,w]   = spm_diffeo('push',bin_les,y,dm_a(1:3));    
    bin_les = spm_field(w,c,[vs_a  1e-6 1e-4 0  3 2]);
    clear w c
    
    fname = fullfile(dat.dir.seg,['wcles' nam{1} '.nii']);
    spm_misc('create_nii',fname,bin_les > opt.clean.les.val,mat_a,[spm_type('float32') spm_platform('bigend')],'Lesion (template)');
    clear bin_les
    
    res.wcles = nifti(fname);
end
clear y

%--------------------------------------------------------------------------    
% Write ML labels to disk
%--------------------------------------------------------------------------
    
if opt.write.ml    
    [bg,fg] = get_par('bg_fg',opt.seg.bg,K);    
    Z       = cat(4,sum(Z(:,:,:,bg),4),Z(:,:,:,fg));    

    [~,ml] = max(Z,[],4);
    ml     = ml - 1;

    fname = fullfile(dat.dir.ml,'ml.nii');
    spm_misc('create_nii',fname,ml,mat_s,[spm_type('float32') spm_platform('bigend')],'ML-labels');    
    clear ml
end
clear Z
       
%--------------------------------------------------------------------------    
% Merge hemisphere classes into just one class
%--------------------------------------------------------------------------

if ~isempty(opt.lesion.hemi)
    
end

%--------------------------------------------------------------------------    
% Show final segmentations
%--------------------------------------------------------------------------

if opt.seg.show
    files = spm_select('FPList',dat.dir.seg,'^c.*\.nii$');
    spm_check_registration(char(obs_fnames, files));
end

%--------------------------------------------------------------------------    
% Clean-up temporary files
%--------------------------------------------------------------------------

if opt.model.clean_up && isfield(dat,'dir') && isfield(dat.dir,'vel')
    rmdir(dat.dir.vel,'s');
end
%==========================================================================

%==========================================================================
function Z = get_final_resp(obs,bf,dat,template,labels,BinWidth,miss,dm,opt)
% FORMAT Z = get_final_resp(obs,bf,dat,template,labels,BinWidth,miss,dm,opt)
%
% Computes subject responsibilities, filling in unobserved values with the
% template.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Parameters
cluster = dat.gmm.cluster;
prop    = dat.gmm.prop;
part    = dat.gmm.part;
K       = max(part.lkp);
const   = spm_gmm_lib('Const', cluster{1}, cluster{2}, miss.L);
ix_tiny = get_par('ix_tiny',dat.population,part.lkp,opt);

% Neighborhood part
lnPzN = gmm_mrf('apply',dat.mrf);

Z = zeros([dm K],'single');
for z=1:dm(3)        
    % Get slice data
    [slice,ix] = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels,BinWidth);

    if dat.mrf.do && numel(lnPzN) > K    
        lnPzNz = double(lnPzN(ix,:));
    else
        lnPzNz = lnPzN;
    end

    % Get responsibilities for a slice
    Z_slice = get_resp_slice(slice,cluster,prop,part,miss,const,lnPzNz,ix_tiny);

    % Go from cluster to tissue responsibilities
    Z_slice = cluster2template(Z_slice,part);    

    % Map slice into full responsibilities
    Z(:,:,z,:) = reshape(Z_slice,[dm(1:2) 1 K]);
end
%==========================================================================

%==========================================================================
function Z = get_resp_slice(slice,cluster,prop,part,miss,const,lnPzN,ix_tiny)
% FORMAT Z = get_resp_slice(slice,cluster,prop,part,miss,const,lnPzN,ix_tiny)
%
% Computes slice of subject responsibilities, filling in unobserved values 
% with the template.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Parameters
MU     = cluster{1}{1};
prec   = cluster{2};
lkp    = part.lkp;
mg     = part.mg;
lntiny = log(eps);

% Multiply bias field with observed data
BX = slice.bf.*slice.obs;

% Mask out where there are no observations and compute ln|bf|    
lnDetbf = log(prod(slice.bf,2));  

% Compute ln(p(x))
lnpX = spm_gmm_lib('Marginal', BX, [{MU} prec], const, {slice.code,miss.L}, slice.bin_var);

for k=1:numel(lkp)
    lnpX(slice.code == 0,k) = log(1e-3);
end

% Compute ln(p(z))
if isempty(slice.template)   
    lnPI         = reshape(prop,[],numel(prop));
    slice.labels = ones(1,size(lnPI,2));
else
    lnPI = log(spm_matcomp('softmax',slice.template,prop) + eps);
    lnPI = lnPI(:,lkp);
end

% Get log of label part
lnPl = log(slice.labels);
lnPl = lnPl(:,lkp);

% MRF part
lnPzN = lnPzN(:,lkp);

% Force responsibilties to zero for ix_tiny classes
lnpX(:,ix_tiny) = lntiny;

% Compute responsibilities
Z = spm_gmm_lib('Responsibility', lnpX,  lnPI,         lnDetbf,   lnPl,log(mg),lnPzN);   

% for k=1:numel(lkp)
%     Z(slice.code == 0,k) = 0;
% end
%==========================================================================

%==========================================================================
function write_dartel(Z,y,mat_s,dm_s,dm_a,mat_a,vx_a,vx_s,pth,nam,opt)
% FORMAT write_dartel(Z,y,mat_s,dm_s,dm_a,mat_a,vx_a,vx_s,pth,nam,opt)
%
% Write segmentations ready for DARTEL import.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

K       = size(Z,4);
[~,fg0] = get_par('bg_fg',opt.seg.bg,K);
tc      = opt.write.tc(:,5);

% Generate mm coordinates of where deformations map from
x = affind(rgrid(dm_s(1:3)),mat_s);

% Generate mm coordinates of where deformation maps to
y1 = affind(y,mat_a);
clear y

% Procrustes analysis to compute the closest rigid-body
% transformation to the deformation, weighted by the
% interesting tissue classes.
[dummy,R] = spm_get_closest_affine(x,y1,Z(:,:,:,fg0(1)));
clear x y1

% Voxel-to-world of original image space
mat0 = R\mat_a; 

fwhm = max(vx_a./vx_s - 1,0.01);
for k=find(tc' == true) 
    % Low pass filtering to reduce aliasing effects in downsampled images,
    % then reslice and write to disk
    tmp1 = decimate(Z(:,:,:,k),fwhm);
    
    % Create NIfTI
    Ni      = nifti;
    Ni.dat  = file_array(fullfile(pth,['rc', num2str(k), nam, '.nii']),...
                         dm_a(1:3),...
                         [spm_type('float32') spm_platform('bigend')],...
                         0,1,0);
    Ni.mat         = mat_a;
    Ni.mat_intent  = 'Aligned';
    Ni.mat0        = mat0;
    Ni.mat0_intent = 'Aligned';
    Ni.descrip     = ['Imported Tissue ' num2str(k)];
    create(Ni);

    for z=1:dm_a(3)
        tmp           = spm_slice_vol(tmp1,mat_s\mat0*spm_matrix([0 0 z]),dm_a(1:2),[1,NaN]);
        Ni.dat(:,:,z) = tmp;
    end
    clear tmp1
end
%==========================================================================

%==========================================================================
function dat = decimate(dat,fwhm)
% FORMAT dat = decimate(dat,fwhm)
%
% Helper function for write_dartel()
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Convolve the volume in memory (fwhm in voxels).
lim = ceil(2*fwhm);
x  = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x  = x/sum(x);
y  = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y  = y/sum(y);
z  = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(dat,dat,x,y,z,-[i j k]);
%==========================================================================

%==========================================================================
function y1 = affind(y0,M)
% FORMAT y1 = affind(y0,M)
%
% Helper function for write_dartel()
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

y1 = zeros(size(y0),'single');
for d=1:3
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
%==========================================================================

%==========================================================================
function x = rgrid(d)
% FORMAT x = rgrid(d)
%
% Helper function for write_dartel()
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3)
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
%==========================================================================