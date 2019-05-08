function show_segmentations(dat,opt)
% FORMAT show_segmentations(dat,opt)
% dat - Subjects data structure
% opt - Options structure
%
% Plot a selection of segmentations.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

figname = '(SPM) Sample: Images, segmentations';

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);  
clf(f);

populations = spm_json_manager('get_populations',dat);
P           = numel(populations);
S0          = numel(dat);

nrows         = min(S0,opt.verbose.mx_rows); 
nrows_per_pop = floor(nrows/P);
K             = opt.template.K;
ncols         = K + 1;

cnt_plots = 1;
for p=1:P
    population0 = populations{p}.name;
    modality    = populations{p}.type;
    
    cnt = 1;
    for s=1:S0
        population = dat{s}.population;

        if strcmp(population0,population)
 
            % img
            nii  = nifti(dat{s}.pth.im2d);    
            img  = single(nii.dat(:,:,:,:));                  
            img  = img';
            img1 = img;
            
            % bf*img
            if isfield(dat{s}.pth,'bfim2d')
                nii  = nifti(dat{s}.pth.bfim2d);    
                img1 = single(nii.dat(:,:,:));
                img1 = img1';
            end       
            
            sb = subplot(nrows,ncols,1 + (cnt_plots - 1)*ncols);
            if strcmpi(modality,'CT')
                imagesc(img1,[0 100]); axis off xy;
            else
                imagesc(img1); axis off xy;
            end
            colormap(sb,gray)        
            
            % Z
            nii = nifti(dat{s}.pth.seg2d);    
            Z   = single(nii.dat(:,:,:,:));
            
            img = [];
            for k=1:size(Z,4)
                img = [img Z(:,:,:,k)'];
            end

            sb = subplot(nrows,ncols,[2:ncols - 1] + (cnt_plots - 1)*ncols);
            imagesc(img,[0 1]); axis off xy;             
            colormap(sb,gray)
                   
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