function dat = init_mrf(dat,opt)
% FORMAT dat = init_mrf(dat,opt)
% dat   - Subjects data structure
% opt   - Options structure
%
% Initialise Markov Random Field related variables:
% * dat.mrf.G    - (Expected) "confusion" matrix
% * dat.mrf.ElnG - (Expected) log-"confusion" matrix
% * dat.mrf.w    - Additional weight accounting for voxel-size
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Parameters
S0 = numel(dat);
K  = opt.template.K;

% Inititalise for all subjects
for s=1:S0
    dat{s}.mrf.do = opt.do.mrf;
    dat{s}.mrf.G  = ones(K);
    
    if ~opt.do.mrf
        continue
    end

    dat{s}.mrf.ml = opt.seg.mrf.ml;
    
    [~,~,vs]     = obs_info(dat{s});    
    dat{s}.mrf.w = single(1./(vs.^2)); 
        
    if opt.seg.mrf.ml
        % Maximum-likelihood MRF
        val_diag = opt.seg.mrf.val_diag;

        G              = (1 - val_diag)/(K - 1)*ones(K);
        G(1 == eye(K)) = val_diag;        
        G              = bsxfun(@rdivide,G,sum(G,2));
        
        dat{s}.mrf.G = G;
    else
        % MRF with prior
        alpha = opt.seg.mrf.alpha;              
        G     = alpha*ones(K);

        ElnG = zeros(K);
        for k=1:K
            ElnG(k,:) = psi(G(k,:)) - psi(sum(G(k,:)));
        end
        
        dat{s}.mrf.G    = G;
        dat{s}.mrf.ElnG = ElnG;
    end
end
%==========================================================================