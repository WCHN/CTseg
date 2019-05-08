function dat = init_lkp(dat,opt)
% FORMAT dat = init_lkp(dat,opt)
% dat   - Subjects data structure
% opt   - Options structure
%
% Initialise mapping between GMM clusters and Template classes:
% * dat.gmm.part.lkp: GMM to Template mapping
% * dat.gmm.part.mg:  Within class weighting
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
K   = opt.template.K;
S0  = numel(dat);   
lkp = 1:K;
mg  = ones(1,K);

for s=1:S0
    dat{s}.gmm.part.lkp = lkp;
    dat{s}.gmm.part.mg  = mg;
end
%==========================================================================