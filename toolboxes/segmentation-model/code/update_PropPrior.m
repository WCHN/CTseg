function model = update_PropPrior(dat,model,opt,it_mod)
% FORMAT model = update_PropPrior(dat,model,opt,it_mode)
% dat    - Subjects data structure
% model  - Model structure
% opt    - Options structure
% it_mod - Current EM iteration
%
% Update the parameters of the Dirichlet prior over tissue weights.
% This is done iteratively using Gauss-Newton optimisation.
%
% Uses:
% model.PropPrior.alpha
% model.PropPrior.norm
% dat.gmm.prop
% opt.model.PropPrior.do
% opt.start_it.do_prop
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if ~opt.model.PropPrior.do
    % If Dirichlet parameters should not be optimised, return
    return;
end

if it_mod < opt.start_it.do_prop
    % If it's too soon for that, return
    return;
end

S0    = numel(dat);            % Number of subjects
alpha = model.PropPrior.alpha; % Previous value of the Dirichlet parameters

% meanLogX: [Kx1] Mean of the log of the observations (mean(log(X)))
meanLogX = 0;
for s=1:S0
    meanLogX = meanLogX + log(spm_matcomp('softmax',dat{s}.gmm.prop) + eps);
end
meanLogX = meanLogX./S0;

alpha    = double(alpha(:));    % Dirichlet parammeters
logalpha = log(alpha);          % Log of Dirichlet parameters
meanLogX = double(meanLogX(:)); % Mean of the log of the observations
K        = numel(alpha);        % Number of classes

% Gauss-Newton iterations
E  = NaN;
for gn=1:100000
    
    % Compute objective function (negtive log-likelihood)
    Eo = E;
    E  = sum(gammaln(alpha)) ...
         - gammaln(sum(alpha)) ...
         - sum((alpha-1).*meanLogX);
%     fprintf('E = %10.6g\n', E);
    if E < Eo && abs(Eo - E) < 1E-7
        % It sometimes overshoots during the first iterations but gets back
        % on track on its own.
        break
    end
    
    % Compute grad/hess
    g = alpha .* ( psi(alpha) - psi(sum(alpha)) - meanLogX );
    H = (alpha * alpha') .* (diag(psi(1,alpha)) - psi(1,sum(alpha)) * ones(K));
    
    H = spm_matcomp('LoadDiag', H);
    
    % update
    logalpha = logalpha - H\g;
    alpha    = exp(logalpha);
end

model.PropPrior.alpha = alpha;
model.PropPrior.norm  = S0*(gammaln(sum(alpha)) - sum(gammaln(alpha)));

% % Show results
% show_PropPrior(dat,model,opt);

% Save updated PropPrior
PropPrior = model.PropPrior;
fname     = fullfile(opt.dir_model,'PropPrior.mat');
save(fname,'PropPrior');
%==========================================================================