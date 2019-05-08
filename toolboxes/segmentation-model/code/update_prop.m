function dat = update_prop(obs,bf,dat,template,labels,scl,dm,PropPrior,miss,opt)

% Parse input
cluster = dat.gmm.cluster;
armijo  = dat.armijo.prop;
part    = dat.gmm.part;
prop    = double(dat.gmm.prop(:))';
alpha   = double(PropPrior.alpha(:))';
% tiny    = get_par('tiny');

% Parameters
const    = spm_gmm_lib('Const', cluster{1}, cluster{2}, miss.L);
K        = size(template,2);
ix_tiny  = get_par('ix_tiny',dat.population,part.lkp,opt);

if numel(alpha) < numel(prop)
    alpha = padarray(alpha, [0 numel(prop) - numel(alpha)], 'post', 'replicate');
end

%--------------------------------------------------------------------------
% GN loop
%--------------------------------------------------------------------------

for gn=1:opt.prop.gnniter 
    
    % Neighborhood part
    lnPzN = gmm_mrf('apply',dat.mrf);

    sumZ = 0; sumPi = 0; sumPi2 = 0;
    for z=1:dm(3)
        
        % Get slices
        [slice,ix] = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels,scl);

        if dat.mrf.do && numel(lnPzN) > K           
            lnPzNz = double(lnPzN(ix,:));
        else
            lnPzNz = lnPzN;
        end
    
        % Compute responsibilities and lb
        Z = gmm_img('slice_resp_and_lb',slice,cluster{1},cluster{2},prop,part,miss,const,lnPzNz,ix_tiny);      
        
        % Map cluster responsibilities to tissue responsibilities
        Z = cluster2template(Z,part);   
        
        if dat.mrf.do        
            dat.mrf.oZ(:,:,z,:) = reshape(uint8((2^8)*Z),[dm(1:2) 1 K]);
        end
    
        % Mask
        msk = sum(Z,2) > 0;        

        [dsumZ, dsumPi, dsumPi2] = suffstat(prop, Z, slice.template, msk);

        sumZ   = sumZ + dsumZ;
        sumPi  = sumPi + dsumPi;
        sumPi2 = sumPi2 + dsumPi2;
    end    
    
    p  = spm_matcomp('softmax',prop);    
    p2 = p'*p;

    % ---------------------------------------------------------------------
    % Compute grad/hess
    g = sumPi - sumZ - (alpha - 1).*(1 - p);
    H = diag(sumPi) - sumPi2 + diag(alpha - 1) .* (diag(p) - p2);
    H = spm_matcomp('LoadDiag', H);

    % ---------------------------------------------------------------------
    % GN step
    Update = H\g';
    Update = Update';

    % ---------------------------------------------------------------------
    % Line search    
    oprop     = prop;   
    oprop_reg = sum((alpha - 1) .* log(spm_matcomp('softmax',oprop) + eps));
    for line_search=1:opt.nline_search.prop

        % GN step
        prop = oprop - armijo * Update;
        
        % Prior part of lower bound
        prop_reg = sum((alpha - 1) .* log(spm_matcomp('softmax',prop) + eps));

        % Z part of lower bound
        [dlb,~,dat.mrf] = gmm_img('img_lb_and_mom',obs,bf,scl,template,labels,prop,cluster{1},cluster{2},miss,part,dm,dat.mrf,ix_tiny,{'prop',oprop});                               

        if (dlb.Z + prop_reg) > (dat.lb.Z(end) + oprop_reg)
            armijo = min(armijo*1.25,1);
            
            dat.lb.Z(end + 1)        = dlb.Z;  
            dat.lb.prop_reg(end + 1) = prop_reg;

            break;
        else
            armijo = armijo*0.5;

            if line_search == opt.nline_search.prop                        
                prop = oprop;            
            end        
        end

    end

    [dat.lb,gain] = check_convergence('prp',dat.lb,gn,opt.verbose.prop,armijo);
    
    if gain < opt.prop.tol
       % Finished updating registration
       break;
    end
end

if opt.verbose.prop >= 3
    % Plot GMM
    [MU,A] = get_mean_prec(cluster);
    spm_gmm_lib('Plot', 'GMM', obs, {MU,A}, spm_matcomp('softmax',prop), part);
end

dat.armijo.prop = armijo;
dat.gmm.prop    = prop;
%==========================================================================

%==========================================================================
function [sumZ, sumPi, sumPi2] = suffstat(w, Z, A, msk)

K    = size(Z,2);
sumZ = sum(double(Z),1);
Pi   = spm_matcomp('softmax',bsxfun(@plus, double(A), double(w)));
for k=1:size(Pi,2)
    Pi(~msk,k) = 0;
end
sumPi2 = zeros(K);
for k=1:K
    sumPi2(k,k) = sum(Pi(:,k).^2);
    for l = k+1:K
        sumPi2(k,l) = sum(Pi(:,k).*Pi(:,l));
        sumPi2(l,k) = sumPi2(k,l);
    end
end
sumPi = sum(Pi,1);
%==========================================================================