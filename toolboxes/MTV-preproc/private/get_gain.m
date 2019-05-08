function gain = get_gain(ll)
% Compute log-likelihood gain, used to determine algorithm covergence
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

gain = abs((ll(end - 1)*(1 + 10*eps) - ll(end))/ll(end));
%==========================================================================