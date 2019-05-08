function show_params(dat,model,opt)

% Get figure (create if does not exist)
%--------------------------------------------------------------------------
figname = '(SPM) Parameters';
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
K           = opt.template.K;
is2d        = model.template.nii.dat.dim(3) == 1;

nrows         = min(S0,opt.verbose.mx_rows); 
nrows_per_pop = floor(nrows/P);
if is2d
    ncols     = 7;
else
    ncols     = 9;
end

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
            
            sb = subplot(nrows,ncols,1 + (cnt_plots - 1)*ncols);
            if strcmpi(modality,'CT')
                imagesc(img1,[0 100]); axis off image xy;
            else
                imagesc(img1); axis off image xy;
            end
            colormap(sb,gray) 
            
            % Histogram (native)      
            subplot(nrows,ncols,2 + (cnt_plots - 1)*ncols);
            bar(dat{s}.verbose.bf.y0,dat{s}.verbose.bf.x0)            
            set(gca,'ytick',[])
            
            % bf*img
            if isfield(dat{s}.pth,'bfim2d')
                nii  = nifti(dat{s}.pth.bfim2d{1});    
                img1 = single(nii.dat(:,:,:));
                img1 = img1';
            end       
            
            sb = subplot(nrows,ncols,3 + (cnt_plots - 1)*ncols);
            if strcmpi(modality,'CT')
                imagesc(img1,[0 100]); axis off image xy;
            else
                imagesc(img1); axis off image xy;
            end
            colormap(sb,gray)                             
                       
            % Histogram (corrected)      
            subplot(nrows,ncols,4 + (cnt_plots - 1)*ncols);
            bar(dat{s}.verbose.bf.y1,dat{s}.verbose.bf.x1)            
            set(gca,'ytick',[])   
                              
            % lower bound            
            sb = subplot(nrows,ncols,5 + (cnt_plots - 1)*ncols);
            plot(dat{s}.lb.sum); axis off;
            
            % v       
            if is2d
                nii = nifti(dat{s}.pth.v2d{1});    
                img = single(nii.dat(:,:));

                sb = subplot(nrows,ncols,6 + (cnt_plots - 1)*ncols);            
                imagesc(img'); axis off image xy;
                colormap(sb,gray)
%                 colorbar('westoutside');
            else
                for i=1:3
                    nii = nifti(dat{s}.pth.v2d{i});    
                    img = squeeze(single(nii.dat(:,:,:)));

                    sb = subplot(nrows,ncols,(6 + i - 1) + (cnt_plots - 1)*ncols);            
                    imagesc(img'); axis off image xy;
                    colormap(sb,gray)
                    
%                     if i == 1
%                         colorbar('westoutside');
%                     end
                end
            end
            
            % affine
            r = dat{s}.reg.r;
            B = opt.reg.B;
            M = spm_dexpm(r,B);
            q = spm_imatrix(M);
            
            q(4:6)  = 180/pi*q(4:6);
            q       = abs(q);
%             q(q==0) = [];
            
            subplot(nrows,ncols,ncols + (cnt_plots - 1)*ncols);
            hold on
            for k=1:numel(q)
                bar(k, q(k));
            end
            box on
            hold off            
            set(gca,'xtick',[])  
            
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