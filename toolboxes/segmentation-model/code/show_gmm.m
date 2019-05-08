function show_gmm(dat,obs)
cluster = dat.gmm.cluster;
prop    = dat.gmm.prop;
part    = dat.gmm.part;

[MU,A] = get_mean_prec(cluster);

spm_gmm_lib('Plot', 'GMM', obs, {MU,A}, spm_matcomp('softmax',prop), part);
%==========================================================================