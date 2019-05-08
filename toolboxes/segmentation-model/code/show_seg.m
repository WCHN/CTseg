function show_seg(obs,Template,prop,Z,dm,modality,title_nam,ticklabels,figname)

if nargin < 9
    figname = '(SPM) Observed, template and responsibilities';
end

C       = size(obs,2);
C0      = C + 2;
K       = size(Template,2);
colors  = hsv(K);

if nargin < 7,                        title_nam = cell(1,C + 1); end
if nargin < 8 || isempty(ticklabels), ticklabels = 1:K; end

if numel(ticklabels) ~= K
    ticklabels = 1:K;
end

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   

obs      = reshape(obs,[dm C]);
if ~isempty(Z)
    Z    = reshape(Z,[dm K]);
end

Template = spm_matcomp('softmax',Template,prop);  
Template = reshape(Template,[dm K]);

if dm(3)==1    
    %----------------------------------------------------------------------
    % 2d 
    %----------------------------------------------------------------------
            
    % Show observed data
    %----------------------------------------------------------------------
    for c=1:C    
        subplot(1,C0,c)
        
        slice = obs(:,:,1,c);
        slice = permute(slice,[2 1 3]);
        if strcmpi(modality,'CT')
            imagesc(slice,[0 100]); axis off xy;
        else
            imagesc(slice); axis off xy;
        end
        colormap(gca,'gray')
        
        if ~isempty(title_nam{c})
            title(title_nam{c})
        end
    end
    
    % Show template
    %----------------------------------------------------------------------
    c = C + 1;
    subplot(1,C0,c)
    
    slice = Template(:,:,1,:);
    slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
    slice = squeeze(slice(:,:,:,:));
    slice = permute(slice,[2 1 3]);
    imagesc(slice); axis off xy;      

    colormap(gca,colors)
    cb = colorbar;
    set(gca, 'clim', [0.5 K+0.5]);
    set(cb, 'ticks', 1:K, 'ticklabels', ticklabels); 
    
    if numel(title_nam) == c && ~isempty(title_nam{c})
        title(title_nam{c})
    end
    
    if ~isempty(Z)
        % Show Z
        %------------------------------------------------------------------
        c = C + 2;
        subplot(1,C0,c)

        slice = Z(:,:,1,:);
        slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
        slice = squeeze(slice(:,:,:,:));
        slice = permute(slice,[2 1 3]);
        imagesc(slice); axis off xy;      

        colormap(gca,colors)
        cb = colorbar;
        set(gca, 'clim', [0.5 K+0.5]);
        set(cb, 'ticks', 1:K, 'ticklabels', ticklabels); 

        if numel(title_nam) == c && ~isempty(title_nam{c})
            title(title_nam{c})
        end
    end    
    
else
    %----------------------------------------------------------------------
    % 3d
    %----------------------------------------------------------------------        
    
    % Show observed data
    %----------------------------------------------------------------------    
    for c=1:C  
        subplot(C0,3,3*(c - 1) + 1)
        
        slice = obs(:,:,floor(dm(3)/2) + 1,c);
        slice = permute(slice,[2 1 3]);
        if strcmpi(modality,'CT')
            imagesc(slice,[0 100]); axis off xy;
        else
            imagesc(slice); axis off xy;
        end
        colormap(gca,'gray')
        
        if ~isempty(title_nam{c})
            title(title_nam{c})
        end
        
        subplot(C0,3,3*(c - 1) + 2)
        
        slice = permute(obs(:,floor(dm(2)/2) + 1,:,c),[3 1 2 4]);  
        if strcmpi(modality,'CT')
            imagesc(slice,[0 100]); axis off xy;
        else
            imagesc(slice); axis off xy;
        end
        colormap(gca,'gray')

        if ~isempty(title_nam{c})
            title(title_nam{c})
        end
        
        subplot(C0,3,3*(c - 1) + 3)
        
        slice = permute(obs(floor(dm(1)/2) + 1,:,:,c),[2 3 1 4]);
        slice = permute(slice,[2 1 3]);
        if strcmpi(modality,'CT')
            imagesc(slice,[0 100]); axis off xy;
        else
            imagesc(slice); axis off xy;
        end
        colormap(gca,'gray')
        
        if ~isempty(title_nam{c})
            title(title_nam{c})
        end     
    end    
    
    % Show template
    %----------------------------------------------------------------------
    c = C + 1;
    subplot(C0,3,3*(c - 1) + 1)
    
    slice = Template(:,:,floor(dm(3)/2) + 1,:);
    slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
    slice = squeeze(slice(:,:,:,:));
    slice = permute(slice,[2 1 3]);
    imagesc(slice); axis off xy;  
    colormap(gca,colors)
    
    if ~isempty(title_nam{c})
        title(title_nam{c})
    end

    subplot(C0,3,3*(c - 1) + 2)

    slice = permute(Template(:,floor(dm(2)/2) + 1,:,:),[3 1 2 4]);
    slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
    slice = squeeze(slice(:,:,:,:));    
    imagesc(slice); axis off xy;  
    colormap(gca,colors)

    if ~isempty(title_nam{c})
        title(title_nam{c})
    end

    subplot(C0,3,3*(c - 1) + 3)

    slice = permute(Template(floor(dm(1)/2) + 1,:,:,:),[2 3 1 4]);
    slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
    slice = squeeze(slice(:,:,:,:));
    slice = permute(slice,[2 1 3]);
    imagesc(slice); axis off xy;  
    colormap(gca,colors)

    if ~isempty(title_nam{c})
        title(title_nam{c})
    end

    cb = colorbar;
    set(gca, 'clim', [0.5 K+0.5]);
    set(cb, 'ticks', 1:K, 'ticklabels', ticklabels);   
    
    if ~isempty(Z)
        % Show Z
        %------------------------------------------------------------------
        c = C + 2;
        subplot(C0,3,3*(c - 1) + 1)

        slice = Z(:,:,floor(dm(3)/2) + 1,:);
        slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
        slice = squeeze(slice(:,:,:,:));
        slice = permute(slice,[2 1 3]);
        imagesc(slice); axis off xy;  
        colormap(gca,colors)

        if numel(title_nam) == c && ~isempty(title_nam{c})
            title(title_nam{c})
        end

        subplot(C0,3,3*(c - 1) + 2)

        slice = permute(Z(:,floor(dm(2)/2) + 1,:,:),[3 1 2 4]);
        slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
        slice = squeeze(slice(:,:,:,:));    
        imagesc(slice); axis off xy;  
        colormap(gca,colors)

        if numel(title_nam) == c && ~isempty(title_nam{c})
            title(title_nam{c})
        end

        subplot(C0,3,3*(c - 1) + 3)

        slice = permute(Z(floor(dm(1)/2) + 1,:,:,:),[2 3 1 4]);
        slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
        slice = squeeze(slice(:,:,:,:));
        slice = permute(slice,[2 1 3]);
        imagesc(slice); axis off xy;  
        colormap(gca,colors)

        if  numel(title_nam) == c && ~isempty(title_nam{c})
            title(title_nam{c})
        end

        cb = colorbar;
        set(gca, 'clim', [0.5 K+0.5]);
        set(cb, 'ticks', 1:K, 'ticklabels', ticklabels);   
    end    
end

drawnow
%==========================================================================