function dat = update_mg(dat,obs,bf,template,dm,labels,scl,miss,opt)

% Parse input
cluster = dat.gmm.cluster;
prop    = dat.gmm.prop;
part    = dat.gmm.part;

% Parameters
const    = spm_gmm_lib('Const', cluster{1}, cluster{2}, miss.L);
mg       = part.mg;
lkp      = part.lkp;
K        = max(lkp);
ix_tiny  = get_par('ix_tiny',dat.population,part.lkp,opt);
 
% Init
[dlb,mom] = gmm_img('init_lb_and_mom',miss);

% Neighborhood part
lnPzN = gmm_mrf('apply',dat.mrf);

for z=1:dm(3) % Loop over slices

    % Get slices
    [slice,ix] = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels,scl);

    if dat.mrf.do && numel(lnPzN) > K           
        lnPzNz = double(lnPzN(ix,:));
    else
        lnPzNz = lnPzN;
    end
    
    % Compute responsibilities and lb
    [Z,dlb,BX] = gmm_img('slice_resp_and_lb',slice,cluster{1},cluster{2},prop,part,miss,const,lnPzNz,ix_tiny,dlb);

    if dat.mrf.do        
        dat.mrf.oZ(:,:,z,:) = reshape(cluster2template(uint8((2^8)*Z),part),[dm(1:2) 1 K]);
    end
        
    % Compute sufficient statistics 
    mom = gmm_img('slice_mom',mom,Z,slice,miss,BX);

end % <-- Loop over slices

lbX = spm_gmm_lib('MarginalSum', mom.SS0, mom.SS1, mom.SS2, cluster{1}, cluster{2}, miss.L, mom.SS2b);    

dat.lb.X(end + 1)       = lbX;     
dat.lb.Z(end + 1)       = dlb.Z;     
if numel(part.mg) > numel(prop)
    dat.lb.mg(end + 1)  = dlb.mg;   
end
if ~isempty(labels)
    dat.lb.lab(end + 1) = dlb.lab;
end
if dat.mrf.do   
    dat.lb.ZN(end + 1)  = gmm_mrf('lowerbound',dat.mrf);
end
    
dat.lb = check_convergence('mg',dat.lb,1,opt.verbose.prop);

MU = cluster{1}{1};
A  = cluster{2}{1};
    
SS0 = spm_gmm_lib('SuffStat', 'infer', mom.SS0, mom.SS1, mom.SS2, {MU,A}, miss.L);

for k=1:numel(lkp)
    tmp   = SS0(lkp == lkp(k));
    mg(k) = (SS0(k) + eps)/sum(tmp + eps); 
end

part.mg = mg;

if opt.verbose.prop >= 3
    % Plot GMM
    [MU,A] = get_mean_prec(cluster);
    spm_gmm_lib('Plot', 'GMM', obs, {MU,A}, spm_matcomp('softmax',prop), part);
end

% Init
dlb = gmm_img('init_lb_and_mom',miss);

% Neighborhood part
lnPzN = gmm_mrf('apply',dat.mrf);

for z=1:dm(3) % Loop over slices

    % Get slices
    [slice,ix] = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels,scl);

    if dat.mrf.do && numel(lnPzN) > K           
        lnPzNz = double(lnPzN(ix,:));
    else
        lnPzNz = lnPzN;
    end
    
    % Compute responsibilities and lb
    [Z,dlb] = gmm_img('slice_resp_and_lb',slice,cluster{1},cluster{2},prop,part,miss,const,lnPzNz,ix_tiny,dlb);
    
    if dat.mrf.do        
        dat.mrf.oZ(:,:,z,:) = reshape(cluster2template(uint8((2^8)*Z),part),[dm(1:2) 1 K]);
    end    
    
%     % Compute sufficient statistics 
%     mom = gmm_img('slice_mom',mom,Z,slice,miss,BX);
end % <-- Loop over slices
   
% lbX = spm_gmm_lib('MarginalSum', mom.SS0, mom.SS1, mom.SS2, cluster{1}, cluster{2}, miss.L, mom.SS2b);    

% dat.lb.X(end + 1)  = lbX;     
% dat.lb.Z(end + 1)  = dlb.Z;   
dat.lb.mg(end + 1) = dlb.mg;   

dat.lb = check_convergence('mg',dat.lb,1,opt.verbose.prop);

dat.gmm.part.mg = mg;
%==========================================================================