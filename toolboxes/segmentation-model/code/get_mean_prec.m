function [MU,A] = get_mean_prec(cluster)
% FORMAT [MU,A] = get_mean_prec(cluster)
% cluster - GMM posterior parameters {MU,b,V,n}
% MU      - Corresponding expected mean
% A       - Corresponding expected precision
%
% Get expected values of mean and precision from posterior Gauss-Wishart 
% parameters.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
MU = cluster{1}{1};
K  = size(MU,2);
W0 = cluster{2}{1};
n0 = cluster{2}{2};
A  = bsxfun(@times,W0,reshape(n0,[1 1 K])); % E[Lambda]
%==========================================================================