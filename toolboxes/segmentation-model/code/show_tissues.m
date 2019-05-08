function show_tissues(dat,model,opt)

% Get figure (create if does not exist)
%--------------------------------------------------------------------------
figname = '(SPM) Tissue classes';
f       = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);  
clf(f);

% Generate figures
%--------------------------------------------------------------------------
populations = spm_json_manager('get_populations',dat);
P           = numel(populations);
S0          = numel(dat);
is2d        = model.template.nii.dat.dim(3) == 1;

if is2d
    % image, segmentations, warped template, proportions
    ncols = 4;
else
    % image x 3, segmentations x 3, warped template x 3, proportions
    ncols = 10;
end

nrows         = min(S0,opt.verbose.mx_rows); 
nrows_per_pop = floor(nrows/P);

rng('default') 
rng(1);    
K      = opt.template.K;
colors = hsv(K);
% colors = colors(randperm(K),:);

titles0   = {'Image','Segmentation','Template','Proportions'};
titles1   = {'','','',''};
add_title = true;

cnt_plots = 1;
for p=1:P
    population0 = populations{p}.name;
    modality    = populations{p}.type;
    
    cnt = 1;
    for s=1:S0
        population = dat{s}.population;

        if strcmp(population0,population)
 
            if add_title
                titles    = titles0;
                add_title = false;
            else
                titles    = titles1;
            end
            
            % Image
            img = get_img(dat{s}.pth.bfim2d);                        
            show_img(img,1,nrows,ncols,cnt_plots,is2d,titles{1},modality);
            
            % Segmentations            
            img = get_img(dat{s}.pth.seg2d,colors);  
            show_img(img,2,nrows,ncols,cnt_plots,is2d,titles{2},modality);
            
            % Warped and scaled template                        
            img = get_img(dat{s}.pth.temp2d,colors);  
            show_img(img,3,nrows,ncols,cnt_plots,is2d,titles{3},modality);
            
            % Proportions
            subplot(nrows,ncols,ncols + (cnt_plots - 1)*ncols);
            PI  = dat{s}.gmm.prop;
            PI = spm_matcomp('softmax',PI);
            mg  = dat{s}.gmm.part.mg;            
            lkp = dat{s}.gmm.part.lkp;
            colors = hsv(max(lkp));
            PI  = mg.*PI(:,lkp);
            hold on
            for k=1:numel(lkp)
                bar(k, PI(k), 'FaceColor', colors(lkp(k),:));
            end
            box on
            hold off            
            axis off
            if ~isempty(titles{4}), title(titles{4},'FontSize',16); end
            
            cnt_plots = cnt_plots + 1;
                        
            if (p < P  && cnt == nrows_per_pop) || ...
               (p == P && cnt_plots - 1 == nrows)
                break
            end
            
            cnt = cnt + 1;
        end
    end
end

drawnow;
%==========================================================================

%==========================================================================
function img = get_img(pth,colors)
img  = cell(1,numel(pth));
perm = {[2 1],[2 1],[2 1]};
for i=1:numel(pth)
    nii    = nifti(pth{i});    
    img{i} = single(nii.dat(:,:,:,:));   
    img{i} = squeeze(img{i});    
    if nargin > 1
        dm     = size(img{i});
        img{i} = reshape(img{i},[dm(1:2) 1 dm(3)]);
        img{i} = spm_gmm_lib('plot', 'cat2rgb', img{i}, colors);   
        img{i} = squeeze(img{i});      
        
        img{i} = permute(img{i},[perm{i} 3]);
    else
        img{i} = permute(img{i},perm{i});
    end
end
%==========================================================================

%==========================================================================
function show_img(img,k,nrows,ncols,cnt_plots,is2d,title_nam,modality)
if nargin < 8, modality = ''; end

if is2d
    make_subplot(img{1},k,nrows,ncols,cnt_plots,title_nam,modality); 
else
    cnt = 1;
    for k=3*(k - 1) + 1:3*k
        if cnt == 2
            make_subplot(img{cnt},k,nrows,ncols,cnt_plots,title_nam,modality);   
        else
            make_subplot(img{cnt},k,nrows,ncols,cnt_plots,'',modality);   
        end
        cnt = cnt + 1;
    end
end
%==========================================================================

%==========================================================================
function make_subplot(img,k,nrows,ncols,cnt_plots,title_nam,modality)
if nargin < 7, modality = ''; end

sb = subplot(nrows,ncols,k + (cnt_plots - 1)*ncols);  
if strcmpi(modality,'CT')
    imagesc(img,[0 100]); axis off image xy;
else
    imagesc(img); axis off image xy;
end
if ~isempty(title_nam), title(title_nam,'FontSize',16); end
colormap(sb,gray)   

% pos = get(gca, 'Position');
% pos(1) = 0.055;
% pos(3) = 0.9;
% set(gca, 'Position', pos)
%==========================================================================