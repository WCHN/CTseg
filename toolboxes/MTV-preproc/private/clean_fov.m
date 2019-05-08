function img = clean_fov(img,dat)
% Zeros voxels outside of the field of view
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

mat1 = dat.mat;      
dm1  = dat.dm; 
msk  = cell(dat.N,3);
for n=1:dat.N
    mat0 = dat.A(n).mat;  
    dm0  = dat.A(n).dm;

    % Get the mapping from Mref to Mmod
    T = mat0\mat1;

    % Use ndgrid to give an array of voxel indices
    [x0,y0,z0] = ndgrid(single(1:dm1(1)),...
                        single(1:dm1(2)),...
                        single(1:dm1(3)));

    % Transform these indices to the indices that they point to in the reference image
    D = cat(4,T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4), ...
              T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4), ...
              T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4));

    if dm1(3) == 1
        D(:,:,:,3) = 1;
    end
    
    % Mask according to whether these are < 1 or > than the dimensions of the reference image.        
    for i=1:3
        msk{n,i} = D(:,:,:,i) >= 1 & D(:,:,:,i) <= dm0(i);
    end
end  

% Generate cleaned up image
for n=1:dat.N
    for i=1:3
        img = msk{n,i}.*img;
    end
end
%==========================================================================