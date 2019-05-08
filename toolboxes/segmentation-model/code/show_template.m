function show_template(model,opt,S0,figname)

if nargin < 4
    figname = '(SPM) Template';
end

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

mu = single(model.template.nii.dat(:,:,:,:));
mu = spm_matcomp('softmax',mu);

title_name = 'Template';
if nargin >= 3
    title_name = [title_name ' (S=' num2str(S0) ')'];
end

alpha = model.PropPrior.alpha;

show_cat_img(mu,alpha,{title_name},opt.model.nam_cls);

deal_figs(model);
%==========================================================================

%==========================================================================
function show_cat_img(img,alpha,title_nam,ticklabels)

if isnumeric(img)
    img = {img};
end 

dm0    = size(img{1});
is2d   = dm0(3) == 1;
K      = dm0(4);
colors = hsv(K);
            
if nargin < 4 || isempty(ticklabels), ticklabels = 1:K; end

if numel(ticklabels) ~= K
    ticklabels = 1:K;
end

if is2d
    N = 2;
else
    N = 4;
end

% Compute mean of Dirichlet distribution
mean_alpha = bsxfun(@rdivide,alpha,sum(alpha));
    
if dm0(3)==1
    % 2d 
    %----------------------------------------------------------------------
              
    subplot(1,N,1)

    slice = img{1}(:,:,1,:);
    slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
    slice = squeeze(slice(:,:,:,:));
    slice = permute(slice,[2 1 3]);
    imagesc(slice); axis off image xy;      

    if ~isempty(title_nam{1})
        title(title_nam{1})
    end
    
    colormap(colors);
%     cb = colorbar;
%     set(gca, 'clim', [0.5 K+0.5]);
%     set(cb, 'ticks', 1:K, 'ticklabels', ticklabels);    
else
    % 3d
    %----------------------------------------------------------------------
    
    subplot(1,N,1)

    slice = img{1}(:,:,floor(dm0(3)/2) + 1,:);
    slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
    slice = squeeze(slice(:,:,:,:));
    slice = permute(slice,[2 1 3]);
    imagesc(slice); axis off image xy;  

    subplot(1,N,2)

    slice = permute(img{1}(:,floor(dm0(2)/2) + 1,:,:),[3 1 2 4]);
    slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
    slice = squeeze(slice(:,:,:,:));    
    imagesc(slice); axis off image xy;  

    if ~isempty(title_nam{1})
        title(title_nam{1})
    end

    subplot(1,N,3)

    slice = permute(img{1}(floor(dm0(1)/2) + 1,:,:,:),[2 3 1 4]);
    slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
    slice = squeeze(slice(:,:,:,:));
    slice = permute(slice,[2 1 3]);
    imagesc(slice); axis off image xy;  

    colormap(colors);
%         cb = colorbar;
%         set(gca, 'clim', [0.5 K+0.5]);
%         set(cb, 'ticks', 1:K, 'ticklabels', ticklabels);      
end

% Draw stacked bars
subplot(1,N,N)
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
    
drawnow
%==========================================================================