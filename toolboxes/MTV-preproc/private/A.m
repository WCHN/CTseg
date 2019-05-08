function X = A(Y,dat,n)
% Forward model (x=Ay)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 3
    X = cell(1,dat.N);
    for n=1:dat.N
        R = dat.A(n).R;
        T = dat.mat\R*dat.A(n).mat;
        y = apply_affine(T,dat.A(n).dm);

        if strcmp(dat.method,'superres')
            X{n} = pushpull('pullc',Y,y,single(dat.A(n).J),double(dat.A(n).win));    
        elseif strcmp(dat.method,'denoise')
            X{n} = spm_diffeo('pullc',Y,y);
        end
        clear y    

        X{n}(~isfinite(X{n})) = 0;
    end
else
    R = dat.A(n).R;
    T = dat.mat\R*dat.A(n).mat;
    y = apply_affine(T,dat.A(n).dm);

    if strcmp(dat.method,'superres')
        X = pushpull('pullc',Y,y,single(dat.A(n).J),double(dat.A(n).win));    
    elseif strcmp(dat.method,'denoise')
        X = spm_diffeo('pullc',Y,y);
    end
    clear y    

    X(~isfinite(X)) = 0;
end
%==========================================================================  