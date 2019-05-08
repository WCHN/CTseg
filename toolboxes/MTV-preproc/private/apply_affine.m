function y = apply_affine(M,dm)
% Apply an affine matrix to a grid
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

[x0,y0,z0] = ndgrid(single(1:dm(1)),...
                    single(1:dm(2)),...
                    single(1:dm(3)));
y          = cat(4,M(1,1)*x0 + M(1,2)*y0 + M(1,3)*z0 + M(1,4), ...
                   M(2,1)*x0 + M(2,2)*y0 + M(2,3)*z0 + M(2,4), ...
                   M(3,1)*x0 + M(3,2)*y0 + M(3,3)*z0 + M(3,4));
if dm(3) == 1
    y(:,:,:,end) = 1;
end
%==========================================================================   