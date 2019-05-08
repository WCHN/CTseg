function ll = get_ll1(use_projmat,inc_bf,y,Nii_x,Nii_b,tau,dat)
% Compute log of likelihood part
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
  
if use_projmat
    % We use the projection matrices (A, At)
        
    ll = 0;
    for n=1:dat.N
        Ay  = A(y,dat,n);
        x   = get_nii(Nii_x(n));
        msk = get_msk(x,Ay);
        ll  = ll - 0.5*tau(n)*sum(( double(x(msk)) - double(Ay(msk))).^2);
%         ll  = ll - 0.5*tau(n)*sum(double(Ay{n}(msk).^2 - 2*Ay{n}(msk).*x{n}(msk)));
    end
else
    % We do not use the projection matrices (A, At)
    ll = 0;
    for n=1:dat.N
        x   = get_nii(Nii_x(n));
        msk = get_msk(x);
        if inc_bf
            bf = exp(get_nii(Nii_b(n)));
            bf = bf(msk);
        else
            bf  = 1;
        end
        ll  = ll - 0.5*tau(n)*sum((double(x(msk)) - double(bf.*y(msk))).^2);
    end
%     ll  = -0.5*tau*sum(double(y(msk).^2 - 2*y(msk).*x{1}(msk)));
end   
%==========================================================================