function set_globals

% Set random seed, in case random numbers are used
rng('default');
rng(1);

% Set boundary conditions 
% BOUND_CIRCULANT 0 - a point that disappears off the right side of field of view, will appear again on the left.
% BOUND_NEUMANN 1   - gradient is zero at the edge of the field of view
spm_diffeo('boundary',0); 
spm_field('boundary',1);
%==========================================================================