function [objval,dat] = init_objval(dat,P)
% FORMAT [objval,dat] = init_objval(dat,[P])
% dat    - Subjects data structure
% P      - Number of populations
% objval - 
%
% Initialise structure holding lower bound parts:
% dat.lb.{
%   sum      - Sum of all parts
%   last     - Previous value of sum
%   X        - Conditional log-likelihood E[lnp(X|...)]
%   lnDetbf  - Normalising factor due to the bias field
%   Z        - KL-divergence of labels
%   MU       - KL-divergence of GMM mean
%   A        - KL-divergence of GMM precision
%   bf_reg   - Prior term of bias field
%   aff_reg  - Prior term of affine reg
%   v_reg    - Prior term of velocity
%   lab      - 
%   prop_reg - Prior term of tissue proportions
%   mg       -
%   ZN       -
% }
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin<2, P = 1; end

objval          = struct;
objval.main     = [];
objval.template = [];
for p=1:P
    objval.hyperpar(p).KL_qVpV = [];  
    objval.hyperpar(p).ElnDetV = [];
end    

S0  = numel(dat);
for s=1:S0
    [~,~,~,C,~,~,~,~] = obs_info(dat{s});

    lb0 = struct('sum', -Inf, 'last', -Inf, ...
                 'X', 0, 'lnDetbf', 0, 'Z', ...
                 0, 'MU', 0, ...
                 'A', 0, 'bf_reg', zeros(1,C), 'aff_reg', 0, 'v_reg', 0, 'lab', 0, 'prop_reg', 0, 'mg', 0, 'ZN', 0);

    dat{s}.lb = lb0;
end
%========================================================================== 