function M = model_slice_gap(M,gap,vs)
% Modify M to model a (potential) slice-gap
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

S = sqrt(sum(M(1:3,1:3).^2));
R = M(1:3,1:3)/diag(S);
S = S - gap./vs;
M = R*diag(S);
%==========================================================================