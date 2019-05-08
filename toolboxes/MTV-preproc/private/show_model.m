function varargout = show_model(varargin)
% Show various parts of the model
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help show_model
    error('Not enough argument. Type ''help spm_prob'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'solution'
        [varargout{1:nargout}] = show_solution(varargin{:});
    case 'll'
        [varargout{1:nargout}] = show_ll(varargin{:});
    case 'mtv'
        [varargout{1:nargout}] = show_mtv(varargin{:});
    case 'rgb'
        [varargout{1:nargout}] = show_rgb(varargin{:});        
    case 'bf'
        [varargout{1:nargout}] = show_bf(varargin{:});               
    otherwise
        help show_model
        error('Unknown function %s. Type ''help show_model'' for help.', id)
end
%==========================================================================

%==========================================================================
function show_solution(use_projmat,modality,Nii,dat)

figname          = '(SPM) MTV solution';
fig              = findobj('Type', 'Figure', 'Name', figname);
if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', fig);  
clf(fig);

C = numel(Nii.x);

if use_projmat
     for c=1:C
        img = single(Nii.y(c).dat(:,:,:));        
        dm  = size(img);
        dm  = [dm 1];
        ix  = floor(dm./2);
        
        if ix(3) > 1
            % 3D
            tmp = squeeze(img(:,:,ix(3)));
            subplot(3,C,(c-1)*C + 1);        
            show_img(tmp,modality);
            colormap(gray); 
            
            tmp = squeeze(img(:,ix(2),:));
            subplot(3,C,(c-1)*C + 2)
            show_img(tmp,modality);
            colormap(gray); 
            
            tmp = squeeze(img(ix(1),:,:));
            subplot(3,C,(c-1)*C + 3)
            show_img(tmp,modality);
            colormap(gray);            
        else
            % 2D
            subplot(1,C,c);        
            show_img(img,modality);
            colormap(gray);             
        end
    end   
else
    nrows = get_maxN(dat) + 1;
    cnt   = 1;
    
    for rw=1:nrows
        for c=1:C
            if rw < nrows
                % Observations
                if rw <= dat(c).N
                    dm  = dat(c).A(rw).dm;
                    z   = round(dm(3)/2);
                    img = single(Nii.x{c}(rw).dat(:,:,z));       
                else
                    img = [];
                end                                
            else
                % Recovered
                dm  = dat(c).dm;
                z   = round(dm(3)/2);
                img = single(Nii.y(c).dat(:,:,z));      
            end

            % 2D
            if ~isempty(img)   
                subplot(nrows,C,cnt)
                show_img(img,modality);
                colormap(gray);             
            end
            
            cnt = cnt + 1;
        end
    end
end

drawnow;
%==========================================================================

%==========================================================================
function show_ll(ll,llpart)

figname          = '(SPM) MTV log-likelihood';
fig              = findobj('Type', 'Figure', 'Name', figname);
if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', fig);  

idx = 0:(numel(ll)-1);
color = {'r' 'b'};

if nargin > 1 && ~isempty(llpart)
    for i=unique(llpart)
        subll = ll;
        subll(llpart ~= i) = NaN;
        prev = find(llpart == i) - 1;
        prev = prev(prev>0);
        subll(prev) = ll(prev);
        plot(idx,subll,[color{i} '-'],'LineWidth',2);
        hold on
    end
    hold off
else
    plot(ll,'r-','LineWidth',2);
end
grid on
title('log-likelihood')

drawnow;
%==========================================================================

%==========================================================================
function show_mtv(mtv)

figname          = '(SPM) MTV prior';
fig              = findobj('Type', 'Figure', 'Name', figname);
if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', fig);  

imagesc3d(mtv); axis off image xy; colormap(parula); colorbar;
title('Scaling')

drawnow;
%==========================================================================

%==========================================================================
function show_rgb(Nii)

figname          = '(SPM) RGB';
fig              = findobj('Type', 'Figure', 'Name', figname);
if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', fig);  

C    = numel(Nii);
dm   = size(Nii(1).dat(:,:,:));
is3d = numel(dm) == 3;

Gmag = imgradient3(single(Nii(1).dat(:,:,:)));
for c=2:C
    if is3d
        Gmag = cat(4,Gmag,imgradient3(single(Nii(c).dat(:,:,:))));
    else
        Gmag = cat(3,Gmag,imgradient3(single(Nii(c).dat(:,:,:))));
    end
end

% Make max value for each channel <= 1
for c=1:C
    if is3d
        Gmag(:,:,:,c) = Gmag(:,:,:,c)/max(reshape(Gmag(:,:,:,c),[],1));
    else
        Gmag(:,:,c) = Gmag(:,:,c)/max(reshape(Gmag(:,:,c),[],1));
    end
end

imagesc3d(Gmag); axis off image xy;
title('RGB')

drawnow;
%==========================================================================

%==========================================================================
function show_bf(dat,Nii,modality)

figname          = '(SPM) Bias field';
fig              = findobj('Type', 'Figure', 'Name', figname);
if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', fig);  
clf(fig)

C      = numel(Nii.y);
C0     = 5;
nrows  = C;
colors = 'hsv';
n      = 1;

for c=1:C

    dm  = dat(c).dm;    
    z   = round(dm(3)/2);
    y   = get_nii(Nii.y(c),z); 

    bfc = get_nii(Nii.b{c}(n),z);
    bf  = exp(bfc);    
        
    x   = get_nii(Nii.x{c}(n),z);
    
    subplot(nrows,C0,(c - 1)*C0 + 1)
    show_img(x,modality);
    colormap(gca,colors)
    colorbar
    if c == 1
        title('x')
    end
    
    subplot(nrows,C0,(c - 1)*C0 + 2)
    show_img(y,modality);
    colormap(gca,colors)
    colorbar
    if c == 1
        title('y')
    end
    
    subplot(nrows,C0,(c - 1)*C0 + 3)
    show_img(bf,'');
    colormap(gca,colors)
    colorbar
    if c == 1
        title('bf')
    end
    
    subplot(nrows,C0,(c - 1)*C0 + 4)
    show_img(bf.*y,modality);
    colormap(gca,colors)
    colorbar
    if c == 1
        title('bf*y')
    end    
    
    subplot(nrows,C0,(c - 1)*C0 + 5)
    show_img(bf.*y - x,modality);
    colormap(gca,colors)
    colorbar
    if c == 1
        title('bf*y - x')
    end    
end

drawnow;
%==========================================================================

%==========================================================================
% Helper functions
%==========================================================================

%==========================================================================
function show_img(img,modality)
if strcmpi(modality,'CT')
    imagesc(img',[0 100]);
else
    imagesc(img');
end
axis xy image off;
%==========================================================================

%==========================================================================
function N = get_maxN(dat)
N = 0;
for c=1:numel(dat)
    if dat(c).N > N
        N = dat(c).N;
    end
end
%==========================================================================