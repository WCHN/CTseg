function [y,eul_its] = make_deformation(v,prm,int_args,Greens,cyc_its)
% FORMAT [y,eul_its] = make_deformation(v,prm,int_args,Greens,cyc_its)
% v        - Initial velocity field
% prm      - Differential operator parameters
% int_args - Number of integration time steps
% Greens   - Fourier transform of Green's function 
% cyc_its  - Number of Full Multigrid cycles and number of relaxation 
%            iterations per cyc
%
% y       - Forward deformation
% eul_its - Number of integration time steps (just for debugging
%           information)
%
% Make a deformation (y) from an initial velocity (v), either by shooting
% or small displacements.
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 5, cyc_its = [3 3]; end

dm = size(v);

if int_args <= 1
    % Small displacements
    id        = cell(3,1);
    [id{1:3}] = ndgrid(single(1:dm(1)),single(1:dm(2)),single(1:dm(3)));
    id        = cat(4,id{:});

    eul_its = [];
    y       = id + v;
else
    if int_args == Inf
        % Pick a suitable number of Euler iterations heuristically
        eul_its = max(double(floor(sqrt(max(max(max( sum(v.^2, 4) ))))) + 1),2);
        eul_its = min(eul_its,12);
    else
        eul_its = int_args;
    end
    
    if eul_its == 1
        % Small displacements
        id        = cell(3,1);
        [id{1:3}] = ndgrid(single(1:dm(1)),single(1:dm(2)),single(1:dm(3)));
        id        = cat(4,id{:});
        y         = id + v;
    else
        % Shoot
        y = spm_shoot3d(v,prm,[eul_its cyc_its],Greens);
    end
end

if dm(3) == 1
    % Data is 2D
    y(:,:,:,3) = 1;
end
%==========================================================================