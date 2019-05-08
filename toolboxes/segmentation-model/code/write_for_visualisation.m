function dat = write_for_visualisation(dm_s,obs,bf,dat,Template,v,labels,scl,miss,nam,prm,opt)

K    = opt.template.K;
ix1  = round(dm_s(1)/2);
ix2  = round(dm_s(2)/2);
ix3  = round(dm_s(3)/2);
is2d = dm_s(3) == 1;

% Responsibilities
Z = get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt);   

if is2d    
    dat.pth.seg2d = {fullfile(opt.dir_seg2d,['seg2d_' nam '.nii'])};
    spm_misc('create_nii',dat.pth.seg2d{1},Z,eye(4),[spm_type('float32') spm_platform('bigend')],'seg2d');        
else
    dat.pth.seg2d    = cell(1,3);
    im               = Z(:,:,ix3,:);    
    dat.pth.seg2d{1} = fullfile(opt.dir_seg2d,['seg2d_3' nam '.nii']);
    spm_misc('create_nii',dat.pth.seg2d{1},im,eye(4),[spm_type('float32') spm_platform('bigend')],'seg2d');            
    im               = Z(:,ix2,:,:);
    dat.pth.seg2d{2} = fullfile(opt.dir_seg2d,['seg2d_2' nam '.nii']);
    spm_misc('create_nii',dat.pth.seg2d{2},im,eye(4),[spm_type('float32') spm_platform('bigend')],'seg2d');          
    im               = Z(ix1,:,:,:);
    dat.pth.seg2d{3} = fullfile(opt.dir_seg2d,['seg2d_1' nam '.nii']);
    spm_misc('create_nii',dat.pth.seg2d{3},im,eye(4),[spm_type('float32') spm_platform('bigend')],'seg2d');          
end
clear Z

% Image (only one channel)
[~,~,~,~,~,~,~,chn_names] = obs_info(dat); 
for i=1:numel(chn_names)
   if strcmpi(chn_names{i},'T1')
       break
   end
end
Im = reshape(obs(:,i),dm_s(1:3));

[x,y]             = hist(Im(:),100);
dat.verbose.bf.x0 = x;
dat.verbose.bf.y0 = y;

dat.pth.im2d = fullfile(opt.dir_seg2d,['im2d_' nam '.nii']);
spm_misc('create_nii',dat.pth.im2d,Im(:,:,ix3),eye(4),[spm_type('float32') spm_platform('bigend')],'im2d');    

% Bias field modulated image
if numel(bf) > 1
    bf = reshape(bf(:,i),dm_s(1:3));
    Im = bf.*Im;
    clear bf
end

[x,y]             = hist(Im(:),100);
dat.verbose.bf.x1 = x;
dat.verbose.bf.y1 = y;

if is2d  
    dat.pth.bfim2d = {fullfile(opt.dir_seg2d,['bfim2d_' nam '.nii'])};
    spm_misc('create_nii',dat.pth.bfim2d{1},Im(:,:,ix3),eye(4),[spm_type('float32') spm_platform('bigend')],'bfim2d');        
else
    dat.pth.bfim2d  = cell(1,3);
    im              = Im(:,:,ix3);    
    dat.pth.bfim2d{1} = fullfile(opt.dir_seg2d,['bfim2d_3' nam '.nii']);
    spm_misc('create_nii',dat.pth.bfim2d{1},im,eye(4),[spm_type('float32') spm_platform('bigend')],'bfim2d');          
    im              = Im(:,ix2,:,:);
    dat.pth.bfim2d{2} = fullfile(opt.dir_seg2d,['bfim2d_2' nam '.nii']);
    spm_misc('create_nii',dat.pth.bfim2d{2},im,eye(4),[spm_type('float32') spm_platform('bigend')],'bfim2d');          
    im              = Im(ix1,:,:,:);
    dat.pth.bfim2d{3} = fullfile(opt.dir_seg2d,['bfim2d_1' nam '.nii']);
    spm_misc('create_nii',dat.pth.bfim2d{3},im,eye(4),[spm_type('float32') spm_platform('bigend')],'bfim2d');       
end

% Jacobian determinants
[~,v] = spm_shoot3d(v,prm,[4 2 2]);
v     = spm_diffeo('det',v);
    
if is2d 
    dat.pth.v2d = {fullfile(opt.dir_seg2d,['v2d_' nam '.nii'])};
    spm_misc('create_nii',dat.pth.v2d{1},v,eye(4),[spm_type('float32') spm_platform('bigend')],'v2d');        
else
    dat.pth.v2d    = cell(1,3);
    im             = v(:,:,ix3,:);    
    dat.pth.v2d{1} = fullfile(opt.dir_seg2d,['v2d_3' nam '.nii']);
    spm_misc('create_nii',dat.pth.v2d{1},im,eye(4),[spm_type('float32') spm_platform('bigend')],'v2d');          
    im             = v(:,ix2,:,:);
    dat.pth.v2d{2} = fullfile(opt.dir_seg2d,['v2d_2' nam '.nii']);
    spm_misc('create_nii',dat.pth.v2d{2},im,eye(4),[spm_type('float32') spm_platform('bigend')],'v2d');      
    im             = v(ix1,:,:,:);
    dat.pth.v2d{3} = fullfile(opt.dir_seg2d,['v2d_1' nam '.nii']);
    spm_misc('create_nii',dat.pth.v2d{3},im,eye(4),[spm_type('float32') spm_platform('bigend')],'v2d');        
end
clear v

% Warped, scaled template    
Template = spm_matcomp('softmax',Template,dat.gmm.prop);
for k=1:K
    Template(miss.C == 0,k) = 0;
end
Template = reshape(Template,[dm_s(1:3) K]);

if is2d 
    dat.pth.temp2d = {fullfile(opt.dir_seg2d,['temp2d_' nam '.nii'])};
    spm_misc('create_nii',dat.pth.temp2d{1},Template,eye(4),[spm_type('float32') spm_platform('bigend')],'temp2d');        
else
    dat.pth.temp2d    = cell(1,3);
    im                = Template(:,:,ix3,:);    
    dat.pth.temp2d{1} = fullfile(opt.dir_seg2d,['temp2d_3' nam '.nii']);
    spm_misc('create_nii',dat.pth.temp2d{1},im,eye(4),[spm_type('float32') spm_platform('bigend')],'temp2d');          
    im                = Template(:,ix2,:,:);
    dat.pth.temp2d{2} = fullfile(opt.dir_seg2d,['temp2d_2' nam '.nii']);
    spm_misc('create_nii',dat.pth.temp2d{2},im,eye(4),[spm_type('float32') spm_platform('bigend')],'temp2d');      
    im                = Template(ix1,:,:,:);
    dat.pth.temp2d{3} = fullfile(opt.dir_seg2d,['temp2d_1' nam '.nii']);
    spm_misc('create_nii',dat.pth.temp2d{3},im,eye(4),[spm_type('float32') spm_platform('bigend')],'temp2d');        
end
%==========================================================================