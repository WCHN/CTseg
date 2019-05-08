function varargout = pushpull(varargin)
% Mex function called for super-resolution
%__________________________________________________________________________
%
% FORMAT [sub,dx,dy,dz] = pushpull('pull', ref, y, J, win)
% ref      - Input image(s)   [n1 n2 n3 n]
% y        - Points to sample [m1 m2 m3 3]
% J        - Jacobian tensor  [1 1 1 3 3]   (J = R*diag(s))
% win(1:3) - Slice selection window (0=dirac,1=gauss,2=rect)
% sub      - Output image     [m1 m2 m3 n]
% dx,dy,dz - Sampled first derivatives [m1 m2 m3 n]
%
% Convolve and sample an image encoded by exponential basis functions.
% sub = (win * ref)(y)
% - ref, y, J and sub are single precision floating point.
% - win is double precision floating point.
% Values sampled outside the field of view of ref are assigned a value
% of NaN.
%
%__________________________________________________________________________
%
% FORMAT [sub,dx,dy,dz] = pushpull('pull', ref, y, J, win)
% ref      - Input image(s)   [n1 n2 n3 n]
% y        - Points to sample [m1 m2 m3 3]
% J        - Jacobian tensor  [1 1 1 3 3]   (J = R*diag(s))
% win(1:3) - Slice selection window (0=dirac,1=gauss,2=rect)
% sub      - Output image     [m1 m2 m3 n]
% dx,dy,dz - Sampled first derivatives [m1 m2 m3 n]
%
% Convolve and sample an image encoded by exponential basis functions.
% sub = (win * ref)(y)
% - ref, y, J and sub are single precision floating point.
% - win is double precision floating point.
% Uses boundary condiditions that wrap around (circulant).
%
%__________________________________________________________________________
%
% FORMAT ref = pushpull('push', sub, y, J, win)
% sub      - Input image     [m1 m2 m3 n]
% y        - Points to sample [m1 m2 m3 3]
% J        - Jacobian tensor  [1 1 1 3 3]   (J = R*diag(s))
% win(1:3) - Slice selection window (0=dirac,1=gauss,2=rect)
% ref      - Output image(s)   [n1 n2 n3 n]
%
% Push values of a function according to a deformation.  Note that the
% deformation should be the same as the one used with 'pull', and 
% therefore inverse of the one used with 'samp' or 'bsplins'. 
% - ref, y, J and sub are single precision floating point.
% - win is double precision floating point.
% Voxels in sub that would be pushed outside the field of view of ref 
% are ignored.
%
%__________________________________________________________________________
%
% FORMAT ref = pushpull('pushc', sub, y, J, win)
% sub      - Input image     [m1 m2 m3 n]
% y        - Points to sample [m1 m2 m3 3]
% J        - Jacobian tensor  [1 1 1 3 3]   (J = R*diag(s))
% win(1:3) - Slice selection window (0=dirac,1=gauss,2=rect)
% ref      - Output image(s)   [n1 n2 n3 n]
%
% Push values of a function according to a deformation.  Note that the
% deformation should be the same as the one used with 'pull', and 
% therefore inverse of the one used with 'samp' or 'bsplins'. 
% - ref, y, J and sub are single precision floating point.
% - win is double precision floating point.
% Uses boundary condiditions that wrap around (circulant).
%
%__________________________________________________________________________
%
% FORMAT b = pushpull('boundary')
% Get the current boundary condition.
% b - boundary condition
%     0 - field wraps around at the boundary, as if the field is on a
%         torus (circulant).  This is typically assumed when using
%         FFTs for convolution etc.
%     1 - Neumann boundary condition.
%     Note that after a `clear functions' in MATLAB, the boundary
%     condition is reset to 0.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging   


%-This is merely the help file for the compiled routine
error('pushpull.c not compiled - use compile_pushpull')
