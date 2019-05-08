function reg = get_reg(dat,fct,opt)
% FORMAT reg = get_reg(dat,fct,opt)
% dat - Subject's data structure (one subject)
% fct - Facts structure obtained from `get_facts`
% opt - Options structure
% reg - Registration structure with fields:
%       * prm    - Regularisation parameters [vs a m b le1 le2]
%       * Greens - Regularisation oerator's Greens functions
%       * v      - Loaded velocity field
%       * y      - Exponentiated warp
%       * Affine - Exponentiated affine matrix
%
% Create registration structure with:
% * loaded velocity
% * exponentiated warp
% * exponentiated affine matrix
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% FFT of Green's function
reg.prm                             = [fct.subj.vs fct.subj.ff*opt.reg.rparam*prod(fct.subj.vs)];
if opt.reg.int_args > 1, reg.Greens = spm_shoot_greens('kernel',fct.subj.dm(1:3),reg.prm);
else,                    reg.Greens = [];
end

% Get initial velocities
if isnumeric(dat.reg.v)
    % Initial velocity stored in array
    reg.v = dat.reg.v;                            
else
    % Initial velocity read from NIfTI
    reg.v = single(dat.reg.v.dat(:,:,:,:));       
end

% Make deformation
reg.y      = make_deformation(reg.v,reg.prm,opt.reg.int_args,reg.Greens);

% Affine matrix
E          = spm_dexpm(dat.reg.r,opt.reg.B);
reg.Affine = fct.templ.mat\E*fct.subj.mat;

%==========================================================================