function dat = segment_subject(dat,model,opt)

%--------------------------------------------------------------------------
% Initialise model parameters
%--------------------------------------------------------------------------

set_globals;

% Population parameters
K          = opt.template.K;                   % Number of tissue classes
mat_a      = model.template.nii.mat;           % Template orientation matrix                  
dm_a       = model.template.nii.dat.dim;       % Template orientation matrix   
PropPrior  = model.PropPrior;                  % Proportion prior parameter
GaussPrior = model.GaussPrior(dat.population); % Intensity prior parameters
it_mod     = opt.model.it;                     % Model iteration

% Subject parameters
[obs,dm_s,mat_s,vs_s,scl,~,~,~,~,nam,subsmp,grd] = get_obs(dat,'mskonlynan',opt.seg.mskonlynan,'samp',opt.seg.samp);
labels                                           = get_labels(dat,opt,opt.seg.samp,subsmp,grd); % Get labels
clear grd

modality = dat.modality{1}.name;             % Imaging modality
I        = size(obs,1);                      % Total number of voxels
ff       = get_ff(vs_s);                     % Fudge factor (see spm_preproc8)   
do_reg   = opt.reg.do_aff || opt.reg.do_nl ; % Update registration?

miss = get_par('missing_struct',obs);

% Bias field
do_bf = opt.bf.do && strcmpi(modality,'MRI'); % Update bias-field?
if do_bf
    % Compute bias-field
    bf = get_bf(dat.bf.chan,dm_s);             
else  
    % Bias field not estimated
    bf = 1;                                    
end

% Adjust lower bound
lblnDetbf                   = bf;
lblnDetbf(isnan(lblnDetbf)) = 1;        
lblnDetbf                   = log(prod(lblnDetbf,2));  
lblnDetbf                   = sum(lblnDetbf);
dat.lb.lnDetbf(end + 1)     = lblnDetbf; % Bias field ln|Bf|
dat.lb.prop_reg(end + 1)    = dot(PropPrior.alpha - 1,log(spm_matcomp('softmax',dat.gmm.prop) + eps)); % Tissue prop prior

% Get initial velocities
%--------------------------------------------------------------------------
if isnumeric(dat.reg.v), v = dat.reg.v;                            
else,                    v = single(dat.reg.v.dat(:,:,:,:));       
end

% FFT of Green's function
prm_v                           = [subsmp.sk.*vs_s ff*opt.reg.rparam];
if opt.reg.int_args > 1, Greens = spm_shoot_greens('kernel',dm_s(1:3),prm_v);
else,                    Greens = [];
end

% Deformation
y = make_deformation(v,prm_v,opt.reg.int_args,Greens);

if opt.verbose.gmm >= 3
    % Verbose stuff
    [~,~,~,~,~,~,~,chn_names] = obs_info(dat);    
    chn_names{end + 1}        = 'Template';
    chn_names{end + 1}        = 'Z';
end

if ~opt.template.do
    % Set initial posteriors to same values as GaussPrior
    dat.gmm.cluster{1}{1} = GaussPrior{1};
    dat.gmm.cluster{1}{2} = GaussPrior{2};
    dat.gmm.cluster{2}{1} = GaussPrior{3};
    dat.gmm.cluster{2}{2} = GaussPrior{4};
end

%--------------------------------------------------------------------------
% Initial alignment of template
%--------------------------------------------------------------------------

if ~opt.template.do && dm_s(3) > 1
    % Align template by mutual information registration
    
    dat = maffreg_template2subject(dat,model,opt);
    
    E      = spm_dexpm(dat.reg.r,opt.reg.B);
    Affine = (model.template.nii.mat\E*mat_s)*subsmp.MT;

    % Warp template to subject      
    Template = warp_template(model,y,Affine);  
elseif it_mod == 1 && opt.reg.do_aff            
    % Align the template by updating the affine parameters until convergence    

    % Affine matrix    
    E      = spm_dexpm(dat.reg.r,opt.reg.B);
    Affine = (mat_a\E*mat_s)*subsmp.MT;

    % Warp template to subject
    Template = warp_template(model,y,Affine);  

    % Disable MRF just here
    do_mrf0    = dat.mrf.do;
    dat.mrf.do = false;
    
    for it_reg=1:opt.reg.nit_init_aff
        
        if opt.verbose.reg >= 3 
            % Show observed data, warped template and current responsibilities
            Z = get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt);
            show_seg(obs,Template,dat.gmm.prop,Z,dm_s,modality,chn_names,opt.model.nam_cls);
            clear Z
        end
        
        [dat,Affine,Template,gain] = update_affine(dat,model,obs,Template,bf,labels,mat_a,mat_s,y,dm_s,scl,miss,opt,subsmp,true);
        
        if gain < opt.reg.init_aff_tol
            % Finished updating registration
            break;
        end
    end
    
    dat.mrf.do = do_mrf0;
else
    % Template already aligned
    
    % Affine matrix    
    E      = spm_dexpm(dat.reg.r,opt.reg.B);
    Affine = (mat_a\E*mat_s)*subsmp.MT;

    % Warp template to subject
    Template = warp_template(model,y,Affine);  
end

if ~strcmpi(modality,'ct') && ~opt.template.do
    % Get posteriors using the aligned template
    dat.gmm.cluster = get_cluster(obs,bf,dm_s,GaussPrior,miss,{Template,dat.gmm.prop,labels,dat.gmm.part},'sort_pars',false);
end

if opt.do.mrf
    % When using an MRF, store previous responsibilities, etc., in the mrf
    % field of dat
    dat.mrf.oZ = reshape(get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt,'use_uint',true),[dm_s K]); 
    % figure; imshow3D(squeeze(dat.mrf.oZ))
    
    if opt.verbose.mrf >= 3
        show_mrf(dat.mrf);
    end
end

%--------------------------------------------------------------------------
% Verbose
%--------------------------------------------------------------------------
    
if opt.verbose.gmm >= 3
    plot_subject_lowerbound(dat.lb);
end
if opt.verbose.gmm >= 3
    % Show observed data, warped template and current responsibilities
    Z = get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt);
    show_seg(obs,Template,dat.gmm.prop,Z,dm_s,modality,chn_names,opt.model.nam_cls);
    clear Z
end
if opt.verbose.gmm >= 4
    % Show GMM
    show_gmm(dat,obs)
end
if opt.verbose.bf >= 3
    % Show bias field
    show_bf_and_ivel(obs,dm_s,bf);
end
if opt.verbose.reg >= 3 
     % Show initial velocities
    show_bf_and_ivel(obs,dm_s,v);
end
if opt.verbose.reg >= 3 || opt.verbose.gmm >= 3  || opt.verbose.bf >= 3
    % Distribute figures
    deal_figs(model);
end

%--------------------------------------------------------------------------
% Start EM segmentation algorithm
%--------------------------------------------------------------------------

if opt.verbose.gmm, tic; end % Start timer

if opt.template.do
    % Only do one within subject iteration
    nit_main = 1;
else
    nit_main = opt.seg.niter;
end

gain = 0;
for it_seg=1:nit_main   
    
    do_prop = opt.prop.do   && (it_mod >= opt.start_it.do_prop || it_seg >= opt.start_it.do_prop);    
    do_mg   = opt.do.mg     && (it_mod >= opt.start_it.do_mg   || it_seg >= opt.start_it.do_mg);
    do_nl   = opt.reg.do_nl && (it_mod >= opt.reg.strt_nl      || it_seg >= opt.reg.strt_nl);
    
    if it_seg == 1 && do_mg && ~opt.template.do && ~isequal(GaussPrior{7},get_par('lkp',modality,opt))
        % Introduce multiple Gaussians per tissue
        [dat,GaussPrior] = introduce_lkp(dat,model,opt,GaussPrior);                        
    end
    
    %----------------------------------------------------------------------
    % UPDATE BIAS FIELD
    %----------------------------------------------------------------------
            
    updt_mg = do_mg & ...
              (it_mod >= opt.start_it.upd_mg || it_seg >= opt.start_it.upd_mg);
    
    for it_bf=1:opt.bf.niter       
        
        % Start with updating GMM parameters
        dat = update_gmm(obs,bf,dat,Template,labels,scl,dm_s,GaussPrior,miss,updt_mg,opt);                   
                              
        if do_bf 
            % Update bias-field parameters    
            [dat,bf] = update_bf(dat,obs,bf,Template,dm_s,labels,miss,opt);             
        end               
    
        if it_bf > 1 && ~((dat.lb.last - ooll) > 2*opt.bf.tol*I)
            break; 
        end
        ooll = dat.lb.last;
        
        if opt.verbose.gmm >= 3                        
            show_seg(obs,Template,dat.gmm.prop, ...
                     get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt), ...
                     dm_s,modality,chn_names,opt.model.nam_cls);            
        end                
    end   
       
    %----------------------------------------------------------------------
    % UPDATE TISSUE PROPORTIONS
    %----------------------------------------------------------------------
    
    if do_prop         
        if opt.template.do
            prop_niter = opt.prop.niter(min(numel(opt.prop.niter),it_mod));
        else
            prop_niter = opt.prop.niter(min(numel(opt.prop.niter),it_seg));            
        end

        for it_prop=1:prop_niter 
            
            % Start with updating GMM parameters
            dat = update_gmm(obs,bf,dat,Template,labels,scl,dm_s,GaussPrior,miss,updt_mg,opt);                               

            % Update tissue proportions
            dat = update_prop(obs,bf,dat,Template,labels,scl,dm_s,PropPrior,miss,opt);          

            if it_prop > 1 && ~((dat.lb.last - oll)  > 2*opt.bf.tol*I)
                break; 
            end
            oll = dat.lb.last;
            
            if opt.verbose.gmm >= 3                        
                show_seg(obs,Template,dat.gmm.prop, ...
                         get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt), ...
                         dm_s,modality,chn_names,opt.model.nam_cls);            
            end
        end   
    end
    
    %----------------------------------------------------------------------
    % UPDATE MRF WEIGHTS
    %----------------------------------------------------------------------
    
    if opt.do.update_mrf && opt.do.mrf && (it_mod >= opt.start_it.do_upd_mrf || it_seg >= opt.start_it.do_upd_mrf)
        for it_mrf=1:opt.seg.mrf.niter  
            
            % Start with updating GMM parameters
            dat = update_gmm(obs,bf,dat,Template,labels,scl,dm_s,GaussPrior,miss,updt_mg,opt);                   

            % Update MRF weights        
            dat = update_mrf(dat,obs,bf,Template,dm_s,labels,scl,miss,opt);                    

            if opt.verbose.gmm >= 3                        
                show_seg(obs,Template,dat.gmm.prop, ...
                         get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt), ...
                         dm_s,modality,chn_names,opt.model.nam_cls);            
            end
        end  
    end
    
    %----------------------------------------------------------------------
    % UPDATE REGISTRATION
    %----------------------------------------------------------------------
        
    if do_reg              
        for it_reg=1:opt.reg.niter 
        
            if opt.reg.do_aff      
                
                % Start with updating GMM parameters
                dat = update_gmm(obs,bf,dat,Template,labels,scl,dm_s,GaussPrior,miss,updt_mg,opt); 
                
                % Update affine parameters                                
                [dat,Affine,Template,gain] = update_affine(dat,model,obs,Template,bf,labels,mat_a,mat_s,y,dm_s,scl,miss,opt,subsmp);               
            
                if opt.verbose.gmm >= 3
                    show_seg(obs,Template,dat.gmm.prop, ...
                             get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt), ...
                             dm_s,modality,chn_names,opt.model.nam_cls);        
                end   
            end
            
            if do_nl
                
                % Start with updating GMM parameters
                dat = update_gmm(obs,bf,dat,Template,labels,scl,dm_s,GaussPrior,miss,updt_mg,opt); 
                
                % Update initial velocities                                                  
                [dat,y,v,Template,Greens,gain,opt] = update_nonlin(dat,model,obs,Template,bf,labels,v,y,Affine,Greens,it_seg,it_mod,scl,miss,dm_s,vs_s,subsmp,opt);                

                if opt.verbose.gmm >= 3
                    show_seg(obs,Template,dat.gmm.prop, ...
                             get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt), ...
                             dm_s,modality,chn_names,opt.model.nam_cls);        
                end    
            end        
            
            if gain < opt.reg.tol
               % Finished updating registration
               break;
            end
        end
    end  
        
    if opt.template.do
        dat = update_gmm(obs,bf,dat,Template,labels,scl,dm_s,GaussPrior,miss,updt_mg,opt);
    end
    
    %----------------------------------------------------------------------
    % CONVERGED?
    %----------------------------------------------------------------------
    
    if opt.verbose.gmm                
        % Print overall convergence        
        fprintf('%3s | %3d | lb = %10.6f | diff = %1.5f | tol = %1.5f\n', 'seg', it_seg, dat.lb.last, gain, 1e-4);
    end
    
    if it_seg >= 10 && ~((dat.lb.last - ooll) > 2*opt.seg.tol*I)
        if opt.verbose.gmm
            fprintf('Segmentation converged in %i iterations and %0.0f seconds.\n',it_seg,toc);
        end
        
        break
    end
end
clear Greens

if opt.verbose.model >= 3
    % Write some results to disk (for visualisation in SegModel.m)
    dat = write_for_visualisation(dm_s,obs,bf,dat,Template,v,labels,scl,miss,nam,prm_v,opt);
end

if do_nl
    % Save initial velocities
    if isnumeric(dat.reg.v), dat.reg.v              = v;
    else,                    dat.reg.v.dat(:,:,:,:) = v;
    end
end
clear v

if opt.template.do && opt.template.load_a_der
    %----------------------------------------------------------------------
    % Compute template derivatives (will be loaded later in
    % update_template)
    %----------------------------------------------------------------------

    if opt.seg.samp > 0
        % Compute responsibilities on original image size, not the
        % sub-sampled version
        clear obs bf Template miss labels
        
        % Get original image(s) and dimensions        
        [obs,dm_s] = get_obs(dat,'samp',0); 
        
        if do_bf
            % Resize bias field
            [x1,y1] = ndgrid(1:dm_s(1),1:dm_s(2),1);
            z1      = 1:dm_s(3);
            dat_tmp = dat;
            for c=1:numel(dat_tmp.bf.chan)
                d3                    = [size(dat_tmp.bf.chan(c).T) 1];
                dat_tmp.bf.chan(c).B3 = spm_dctmtx(dm_s(3),d3(3),z1);
                dat_tmp.bf.chan(c).B2 = spm_dctmtx(dm_s(2),d3(2),y1(1,:)');
                dat_tmp.bf.chan(c).B1 = spm_dctmtx(dm_s(1),d3(1),x1(:,1));
            end        
            clear x1 y1 z1
            bf = get_bf(dat_tmp.bf.chan,dm_s);
            clear dat_tmp
        else
            bf = 1;
        end
        
        % Resize deformation
        y = resize_def(y,dm_s,subsmp);
        
        % Warp template
        E        = spm_dexpm(dat.reg.r,opt.reg.B); 
        Affine   = model.template.nii.mat\E*mat_s;               
        Template = warp_template(model,y,Affine);
        
        % Missing data struct
        miss = get_par('missing_struct',obs);
        
        % Labels
        labels = get_labels(dat,opt);
    end   
    
    % Get responsibilities    
    Z = get_resp(obs,bf,dat,Template,labels,scl,miss,dm_s,opt);   
    
%     % Some verbose
%     if opt.verbose.gmm >= 3
%         % Show observed data, warped template and current responsibilities
%         show_seg(obs,Template,dat.gmm.prop,Z,dm_s,modality,chn_names,opt.model.nam_cls);
%     end
%     if opt.verbose.reg >= 3 
%          % Show initial velocities
%         show_bf_and_ivel(obs,dm_s,v);
%     end
%     if opt.verbose.bf >= 3
%         % Show bias field
%         show_bf_and_ivel(obs,dm_s,bf);
%     end
    clear Template bf labels

    % Push responsibilities from subject space to template space
    y = spm_warps('transform',Affine,y);
    if dm_s(3) == 1
        y(:,:,:,3) = 1;
    end
    
    % Push responsibilities from subject space into template space
    Z = push_responsibilities(Z,y,dm_a(1:3));                      
    clear y
    
    % Compute template gradients and Hessian (in template space)
    a               = rotate_template(model.template.nii,opt);
    [gr,H,dll]      = diff_template(a,Z,dat.gmm.prop,opt);         
    dat.template.ll = dll;                        
    clear a Z
    
    % Write derivatives to disk                                    
    spm_misc('create_nii',dat.template.pth_gr,gr,mat_a,[spm_type('float32') spm_platform('bigend')],'a_gr'); 
    spm_misc('create_nii',dat.template.pth_H,H,mat_a,[spm_type('float32') spm_platform('bigend')],'a_H');  
    clear gr H
end

if isfield(dat.mrf,'oZ') && opt.template.do
    % Remove oZ field from dat.mrf
    dat.mrf = rmfield(dat.mrf,'oZ');
end   
%==========================================================================
