function dat = update_gmm(obs,bf,dat,Template,labels,scl,dm,pr,miss,do_mg,opt)

% Parse input
IterMax   = opt.gmm.niter;
Tolerance = opt.gmm.tol;
Verbose   = opt.verbose.gmm;
part      = dat.gmm.part;
lb        = dat.lb;
cluster   = dat.gmm.cluster;
prop      = dat.gmm.prop;
ix_tiny   = get_par('ix_tiny',dat.population,part.lkp,opt);

% Update GMM parameters
%----------------------------------------------------------------------
[cluster,prop,lb,mom,dat.mrf,part] = gmm_img('update',obs,bf,cluster,prop,Template,lb,dm,miss,part,dat.mrf,ix_tiny,do_mg, ...
                                        'GaussPrior',pr(1:4),'Labels',labels, ...
                                        'BinWidth',scl, ...
                                        'IterMax',IterMax,'Verbose',Verbose,...
                                        'Constrained',pr{6},...
                                        'Tolerance',Tolerance);        

dat.gmm.cluster = cluster;
dat.gmm.prop    = prop;
dat.gmm.mom     = mom;
dat.lb          = lb;
dat.gmm.part    = part;

if opt.verbose.gmm >= 3    
    [MU,A] = get_mean_prec(dat.gmm.cluster);
    spm_gmm_lib('Plot', 'GMM', obs, {MU,A}, spm_matcomp('softmax',dat.gmm.prop), dat.gmm.part);
end
%==========================================================================