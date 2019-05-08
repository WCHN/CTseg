function plot_model_lb(dat,model,it,opt)
% FORMAT plot_model_lb(dat,model,it,opt)
% dat   - Subjects data structure
% model - Model structure
% it    - Current EM iteration
% opt   - Options structure
%
% Plot the lower bound (for each population) and its different parts:
% * Template posterior
% * Template prior
% * Template likelihood
% * GMM constraint
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Get figure (create if it does not exist)
figname = '(SPM) Model lower bound';
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

P = 0;
if opt.gmm.GaussPrior.constrained
    populations = spm_json_manager('get_populations',dat);
    P           = numel(populations);
end

% Plot
nfigs = 4 + P;
nrows = floor(sqrt(nfigs));
ncols = ceil(nfigs/nrows);             

if it>2
    x = 3;
else
    x = 1;
end

subplot(nrows,ncols,1)
plot(model.lb(x:end))
title(['Model (iter=' num2str(it) ')'])

subplot(nrows,ncols,2)
plot(model.template.objval.post(x:end))
title('-ln(p(a|.))')

subplot(nrows,ncols,3)
plot(model.template.objval.likel(x:end))
title('-ln(p(.|a))')

subplot(nrows,ncols,4)
plot(model.template.objval.pr(x:end))
title('-ln(p(a))')

if opt.gmm.GaussPrior.constrained
    for p=1:P
        name = populations{p}.name;
        pr   = model.GaussPrior(name);

        subplot(nrows,ncols,4 + p)
        plot(pr{6}.KL_qVpV(x:end))
        title(['KL(qV|pV) (' name ')'])
    end
end

drawnow;

deal_figs(model);
%==========================================================================