function [opt,nscl,oscl] = modify_opt(opt,iter)
% FORMAT opt = modify_opt(opt,iter)
% opt       - Options structure
% iter      - Current EM iteration
% nscl      - Current registration scaling value
% oscl      - Previous registration scaling value
%
% Modify options based on the current iteration.
% This allows to have parameters change every few iterations:
% - decrease regularisation (template, deformation)
% - increase precision (GMM iterations, shooting integration steps)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Get scheduling struct
sched = opt.sched;

if isfield(opt,'template') && isfield(opt.template,'reg0')
    % Template regularisation decreases with iteration number
    opt.template.reg  = [opt.template.reg0(1:2) sched.a(min(numel(sched.a),max(iter - opt.reg.strt_nl + 1,1)))*opt.template.reg0(3)];        
end

if isfield(opt,'gmm')
    % Number of GMM iterations increases with iteration number
    opt.gmm.niter     = sched.gmm(min(numel(sched.gmm),iter));
end

nscl = 1;
oscl = 1;
if isfield(opt,'reg')
    % Non-linear registration regularisation decreases with iteration number        
    opt.reg.rparam = opt.reg.rparam0;
    
    ix             = min(numel(sched.reg),max(iter - opt.reg.strt_nl + 1,1));
    opt.reg.rparam = sched.reg(ix)*opt.reg.rparam;     
    
    nscl = sched.reg(ix);
    if ix == 1
        oscl = sched.reg(ix);
    else
        oscl = sched.reg(ix - 1);
    end
end

if isfield(opt,'reg')
    % Integration parameters increases with iteration number, i.e., starts out
    % doing small deformation, then going full diffeomorphic
    ix                = min(numel(sched.eul),max(iter - opt.reg.strt_nl + 1,1));
    opt.reg.int_args  = sched.eul(ix);
end
%==========================================================================