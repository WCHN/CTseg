function show_bf_and_ivel(obs,dm,varargin)

% Get figure (create if it does not exist)
figname = '(SPM) Bias field and initial velocity';
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   

C      = size(obs,2);
z      = floor(dm(3)/2) + 1;
ix     = ix_slice(z,prod(dm(1:2)));
colors = 'hsv';
C0     = C + 1;

if numel(size(varargin{1})) == 2
    % Show bias field
    %----------------------------------------------------------------------        

    bf = varargin{1};
    if bf == 1
        bf = ones(size(obs));
    end

    mn1 = min(obs,[],1);
    mn2 = min(bf.*obs,[],1);
    mn  = [mn1; mn2];
    mn  = min(mn,[],1);

    mx1 = max(obs,[],1);
    mx2 = max(bf.*obs,[],1);
    mx  = [mx1; mx2];
    mx  = max(mx,[],1);

    zobs = reshape(obs(ix,:),[dm(1:2) C]);
    zbf  = reshape(bf(ix,:),[dm(1:2) C]);   

    for c=1:C
        subplot(C0,3,1 + (c-1)*3);                    
        imagesc((zobs(:,:,c))'); axis xy off; colorbar
        colormap(gca,colors)
%         caxis([mn(c) mx(c)]);
        title('obs')

        subplot(C0,3,2 + (c-1)*3);                    
        imagesc((zbf(:,:,c))'); axis xy off; colorbar
        colormap(gca,colors)
        title('bf')

        subplot(C0,3,3 + (c-1)*3);                    
        imagesc((zbf(:,:,c).*zobs(:,:,c))'); axis xy off;  colorbar
        colormap(gca,colors)
%         caxis([mn(c) mx(c)]);
        title('bf*obs')
    end
elseif numel(size(varargin{1})) == 4
    % Show initial velocities
    %----------------------------------------------------------------------
    
    v = varargin{1};
    c = C + 1;
    
    subplot(C0,3,1 + (c-1)*3);                    
    imagesc(v(:,:,z,1)'); axis xy off; colorbar
    colormap(gca,'gray')
    title('v1')

    subplot(C0,3,2 + (c-1)*3);                    
    imagesc(v(:,:,z,2)'); axis xy off; colorbar
    colormap(gca,'gray')
    title('v2')

    subplot(C0,3,3 + (c-1)*3);                    
    imagesc(v(:,:,z,3)'); axis xy off;  colorbar
    colormap(gca,'gray')
    title('v3')    
end

drawnow
%==========================================================================