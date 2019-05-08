function show_registration(dat,opt)
% FORMAT show_registration(dat,opt)
% dat - Subjects data structure
% opt - Options structure
%
% Plot a selection of segmentations.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

figname = '(SPM) Sample: Images, warped template, velocities, affine matrices';

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
ncols         = K + 5;

cnt_plots = 1;
for p=1:P
    population0 = populations{p}.name;
    modality    = populations{p}.type;
    
    cnt = 1;
    for s=1:S0
        population = dat{s}.population;

        if strcmp(population0,population)
            
            % img
            nii = nifti(dat{s}.pth.im2d);    
            img = single(nii.dat(:,:,:,:));
            img = img';       
            
            sb = subplot(nrows,ncols,1 + (cnt_plots - 1)*ncols);
            if strcmpi(modality,'CT')
                imagesc(img,[0 100]); axis off xy;
            else
                imagesc(img); axis off xy;
            end
            colormap(sb,gray)

            % warped and scaled template
            nii = nifti(dat{s}.pth.temp2d);    
            Z   = single(nii.dat(:,:,:,:));
            
            img = [];
            for k=1:size(Z,4)
                img = [img Z(:,:,:,k)'];
            end

            sb = subplot(nrows,ncols,[2:ncols - 4] + (cnt_plots - 1)*ncols);
            imagesc(img,[0 1]); axis off xy;             
            colormap(sb,gray)
                                    
            % v       
            nii = nifti(dat{s}.pth.v2d);    
            img = single(nii.dat(:,:,:,:));
            dm  = size(img);
            img = reshape(img,[dm(1) 3*dm(2)]);
            
            sb = subplot(nrows,ncols,[ncols - 3:ncols - 1] + (cnt_plots - 1)*ncols);            
            imagesc(img); axis off xy;
            colormap(sb,gray)
            
            % affine
            img = spm_dexpm(dat{s}.reg.r,opt.reg.B);
            
            sb = subplot(nrows,ncols,ncols + (cnt_plots - 1)*ncols);            
            imagesc(img); axis off;
            colormap(sb,jet)            
            
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