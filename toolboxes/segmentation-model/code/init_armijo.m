function dat = init_armijo(dat)
% FORMAT dat = init_armijo(dat)
% dat   - Subjects data structure
% 
% Initialise armijo factors for all Gauss-Newton updated variables:
% * dat.armijo.bf
% * dat.armijo.prop
% * dat.armijo.nl
% * dat.armijo.prop
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
S0 = numel(dat);
for s=1:S0
    [~,~,~,C] = obs_info(dat{s});
    
    dat{s}.armijo.bf   = ones(1,C);
    dat{s}.armijo.aff  = 1;
    dat{s}.armijo.nl   = 1;
    dat{s}.armijo.prop = 1;
end
%==========================================================================