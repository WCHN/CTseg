function dat = gmm_from_priors(dat,model,labels,opt)
[obs,dm_s,mat_s] = get_obs(dat,'do_scl',true,'mskonlynan',opt.seg.mskonlynan); % Do not subsample!
GaussPrior       = model.GaussPrior(dat.population);
miss             = get_par('missing_struct',obs);

% Get affine matrix parameters by mutual information registration
if dm_s(3) > 1                
    [~,P] = maffreg_template2subject(dat,model);

    % Do a least squares fit on the matrix logs to map from 12-parameter
    % maffreg representation to number of basis functions in opt.reg.B
    M1 = P2M1(P);
    e  = eig(M1);
    if isreal(e) && any(e<=0), disp('Possible problem!'); disp(eig(M1)); end
    B  = reshape(opt.reg.B,[16 size(opt.reg.B,3)]);
    P  = B\reshape(real(logm(M1)),[16 1]);
    
    dat.reg.r = P; % Set affine parameter
end

E      = spm_dexpm(dat.reg.r,opt.reg.B);
Affine = model.template.nii.mat\E*mat_s;

% Warp template to subject      
template = warp_template(model,spm_warps('identity',dm_s),Affine);

if 0     
    [~,~,~,~,~,~,~,chn_names] = obs_info(dat);    
    chn_names{end + 1}        = 'Template';
    show_seg(obs,template,dat.gmm.prop,[],dm_s,'',chn_names,opt.model.nam_cls);
end

% Compute suffstats from template and GaussPrior -> then GMM parameters
prop            = dat.gmm.prop;
part            = dat.gmm.part;
dat.gmm.cluster = get_cluster(obs,dm_s,GaussPrior,miss,{template,prop,labels,part},'sort_pars',false);
%==========================================================================

