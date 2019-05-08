function set_boundary_conditions
% Set model boundary conditions
% 0: BOUND_CIRCULANT
% A point that disappears off the right side of field of view, will appear again on the left.
% 1: BOUND_NEUMANN
% Gradient is zero at the edge of the field of view
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

rng('default');
rng(1);

spm_field('boundary',1);
pushpull('boundary',1); 
spm_diffeo('boundary',1); 
%==========================================================================   