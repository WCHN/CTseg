function Y = At(X,dat,tau,n)  
% Adjoint of forward model (y=A'x)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 3, tau = ones(1,dat.N); end
    
if nargin < 4
    Y = single(0);   
    for n=1:dat.N      
        R = dat.A(n).R;
        T = dat.mat\R*dat.A(n).mat;
        y = apply_affine(T,dat.A(n).dm);

        if strcmp(dat.method,'superres')
            tmp = pushpull('pushc',X{n},y,single(dat.A(n).J),double(dat.A(n).win),double(dat.dm));     
        elseif strcmp(dat.method,'denoise')
            tmp = spm_diffeo('pushc',X{n},y,dat.dm);
        end
        clear y

        tmp(~isfinite(tmp)) = 0;

        Y = Y + tau(n).*tmp;           
    end
else
    R = dat.A(n).R;
    T = dat.mat\R*dat.A(n).mat;
    y = apply_affine(T,dat.A(n).dm);

    if strcmp(dat.method,'superres')
        tmp = pushpull('pushc',X,y,single(dat.A(n).J),double(dat.A(n).win),double(dat.dm));     
    elseif strcmp(dat.method,'denoise')
        tmp = spm_diffeo('pushc',X,y,dat.dm);
    end
    clear y

    tmp(~isfinite(tmp)) = 0;

    Y = tau(n).*tmp;      
end
%==========================================================================