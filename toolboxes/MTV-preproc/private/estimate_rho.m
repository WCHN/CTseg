function rho = estimate_rho(tau,lam)
% Estimate rho (this value seems to lead to reasonably good convergence)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging    

C       = numel(tau);
all_tau = [];
for c=1:C
    N = numel(tau{c});
    for n=1:N
        all_tau = [all_tau tau{c}(n)];
    end
end
rho = sqrt(mean(all_tau))/mean(lam);        
%==========================================================================    