function a = rotate_template(nii_a,opt)
% Initial starting estimates for updating template
dm = nii_a.dat.dim;
a  = zeros([dm(1:3),dm(4) - 1],'single');
for z=1:dm(3) % Loop over planes
    sz = nii_a.dat(:,:,z,:);
    for j1=1:(dm(4) - 1)
        az = zeros(dm(1:2));
        for j2=1:dm(4)
            az = az + opt.template.R(j2,j1)*sz(:,:,j2); % Note the rotation
        end
        a(:,:,z,j1) = az;
    end
end
%==========================================================================