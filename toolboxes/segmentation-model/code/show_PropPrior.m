function show_PropPrior(dat,model,opt,figname)

if nargin < 4
    figname = '(SPM) Tissue proportions';
end

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);
 
ticklabels = opt.model.nam_cls;
if isempty(ticklabels)
    K          = opt.template.K;
    ticklabels = 1:K;
else
    K          = numel(ticklabels);    
end

% Get mean of softmaxed weights
S0        = numel(dat);
mean_Prop = 0;
for s=1:S0
    mean_Prop = mean_Prop + spm_matcomp('softmax',dat{s}.gmm.prop);
end
mean_Prop = mean_Prop./S0;

% subplot(121)
% b = bar([mean_Prop(:)'; nan(1,numel(mean_Prop))], 'Stacked');
% set(gca,'xtick',1);
% xlim([0.99 1.01]);
% ylim([0 1]);
% title('mean(w)')
% 
% % Set color of bars
% colors = hsv(numel(b));
% for i=1:numel(b)
%     b(i).FaceColor = colors(i,:);
% end
% 
% colormap(colors);
% cb =  colorbar;
% set(gca, 'clim', [0.5 K+0.5]);
% set(cb, 'ticks', 1:K, 'ticklabels', ticklabels); 

% Compute mean of Dirichlet distribution
alpha      = model.PropPrior.alpha;
mean_alpha = bsxfun(@rdivide,alpha,sum(alpha));

% Draw stacked bars
% subplot(122)
b = bar([mean_alpha(:)'; nan(1,numel(mean_alpha))], 'Stacked');
set(gca,'xtick',1);
xlim([0.99 1.01]);
ylim([0 1]);
title('mean(\alpha)')

% Set color of bars
colors = hsv(numel(b));
for i=1:numel(b)
    b(i).FaceColor = colors(i,:);
end

colormap(colors);
cb =  colorbar;
set(gca, 'clim', [0.5 K+0.5]);
set(cb, 'ticks', 1:K, 'ticklabels', ticklabels); 

drawnow;

deal_figs(model);
%==========================================================================