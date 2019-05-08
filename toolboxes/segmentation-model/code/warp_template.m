function [Template,y] = warp_template(model,y,Affine)
nii  = model.template.nii;
bg   = model.template.bg;
dm_a = nii.dat.dim;     
K    = dm_a(4);
 
y    = spm_warps('transform',Affine,y);
dm_s = size(y);
if dm_a(3) == 1
    y(:,:,:,3) = 1;
end

% Template0 = single(nii.dat(:,:,:,:));
Template  = zeros([dm_s(1:3) K],'single');
for k=1:K
    % Warp template to subject
    Template(:,:,:,k) = spm_diffeo('pull',single(nii.dat(:,:,:,k)),y);
end        
% clear Template0

Template = reshape(Template,[prod(dm_s(1:3)) K]);

% For masking outside FOV
msk1 =   y(:,:,:,1) >= 1 & y(:,:,:,1) <= dm_a(1) ...
       & y(:,:,:,2) >= 1 & y(:,:,:,2) <= dm_a(2) ...
       & y(:,:,:,3) >= 1 & y(:,:,:,3) <= dm_a(3);
msk2 =   y(:,:,:,3) < 1;

% Set values of outside FOV voxels
for k=1:K
    Template(~msk1,k) = bg{2}(k); 
    Template(msk2,k)  = bg{1}(k); 
end
%==========================================================================