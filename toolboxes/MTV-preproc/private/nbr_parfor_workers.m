function nw = nbr_parfor_workers
% Get number of CPU cores
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

c  = parcluster('local');
nw = c.NumWorkers;
%==========================================================================