function [dat,Affine,template,gain] = update_affine(dat,model,obs,template,bf,labels,mat_a,mat_s,y,dm,scl,miss,opt,subsamp,comp_lb)

if nargin < 15, comp_lb = false; end

% Parse input
r       = dat.reg.r;
prop    = dat.gmm.prop;
cluster = dat.gmm.cluster;
armijo  = dat.armijo.aff;
B       = opt.reg.B;
verbose = opt.verbose.reg;
part    = dat.gmm.part;

% Parameters
K       = size(template,2);
Nr      = numel(r);
[E,dE]  = spm_dexpm(r,B);
Affine  = (mat_a\E*mat_s)*subsamp.MT;
const   = spm_gmm_lib('Const', cluster{1}, cluster{2}, miss.L);
ix_tiny = get_par('ix_tiny',dat.population,part.lkp,opt);

%--------------------------------------------------------------------------
% Compute objective function and its first and second derivatives
%--------------------------------------------------------------------------

y1 = spm_warps('transform',Affine,y);
if dm(3) == 1
    y1(:,:,:,3) = 1;
end
dA = zeros(4,4,Nr);
for i1=1:Nr
    dA(:,:,i1) = double(mat_a\dE(:,:,i1)*mat_s);
end
lkp    = [1 4 5; 4 2 6; 5 6 3];
ga     = zeros([Nr, 1]);
Ha     = zeros([Nr,Nr]);
   
[dlb,mom] = gmm_img('init_lb_and_mom',miss);

% Neighborhood part
lnPzN = gmm_mrf('apply',dat.mrf);

Template0 = single(model.template.nii.dat(:,:,:,:));
for z=1:dm(3)

    % Get slices
    [slice,ix] = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels,scl);
    
    if dat.mrf.do && numel(lnPzN) > K           
        lnPzNz = double(lnPzN(ix,:));
    else
        lnPzNz = lnPzN;
    end
    
    if comp_lb
        % Compute responsibilities and lb
        [Z,dlb,BX] = gmm_img('slice_resp_and_lb',slice,cluster{1},cluster{2},prop,part,miss,const,lnPzNz,ix_tiny,dlb);     

        % Compute sufficient statistics 
        mom = gmm_img('slice_mom',mom,Z,slice,miss,BX);
    else
        % Compute responsibilities and lb
        [Z,dlb] = gmm_img('slice_resp_and_lb',slice,cluster{1},cluster{2},prop,part,miss,const,lnPzNz,ix_tiny,dlb);  
    end
    
    % Map cluster responsibilities to tissue responsibilities
    Z = cluster2template(Z,part);       

    % Get gradient and Hessian of multinomial objective function
    Z = reshape(Z,[dm(1:2) 1 K]);
    
    if dat.mrf.do        
        dat.mrf.oZ(:,:,z,:) = uint8((2^8)*Z);
    end
    
    [gz,Hz] = mnom_objfun_slice(Z,Template0,y1,z,prop);
    gz      = double(gz);
    Hz      = double(Hz);
    
    % Compute affine gradient and Hessian
    dAff = cell(Nr,3);
    for i1=1:Nr
        for d1=1:3
            tmp = dA(d1,1,i1)*double(y(:,:,z,1)) + dA(d1,2,i1)*double(y(:,:,z,2)) + dA(d1,3,i1)*double(y(:,:,z,3)) + dA(d1,4,i1);
            dAff{i1,d1} = tmp(:);
        end
    end
    for d1=1:3
        tmp = gz(:,:,d1);
        tmp = tmp(:)';
        for i1=1:Nr
            ga(i1) = ga(i1) + tmp*dAff{i1,d1};
        end
    end

    % Could probably be speeded up by computing Hessian for changes to elements of Affine,
    % and then transforming this Hessian to obtain changes w.tr.t. rigid parameters.
    for d1=1:3
        for d2=1:3
            tmp1 = Hz(:,:,lkp(d1,d2));
            tmp1 = tmp1(:);
            for i1=1:Nr
                tmp2 = (tmp1.*dAff{i1,d1})';
                for i2=i1:Nr % Only compute upper/lower triangle of symmetric matrix
                    Ha(i1,i2) = Ha(i1,i2) + tmp2*dAff{i2,d2};
                end
            end
        end
    end
end
clear y1 Z gZ Hz Template0

if comp_lb
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
end

% Fill in missing triangle
for i1=1:Nr
    for i2=(i1+1):Nr
        Ha(i2,i1) = Ha(i1,i2);
    end
end

Ha = spm_matcomp('LoadDiag', Ha);

% Compute GN step
ICO    = opt.reg.aff_reg_ICO;
Update = reshape((Ha + ICO)\(ga + ICO*r(:)),size(r));

% Linesearch
%--------------------------------------------------------------------------
oTemplate    = template;
or           = r;
nline_search = opt.nline_search.aff;
for line_search=1:nline_search
    
    % Update affine matrix
    r      = r - armijo*Update;
    E      = spm_dexpm(r,B);
    Affine = (mat_a\E*mat_s)*subsamp.MT;
    
    % Compute new lower bound
    aff_reg = double(-0.5*r(:)'*ICO*r(:));    

    template = warp_template(model,y,Affine);        

    [dlb,~,dat.mrf] = gmm_img('img_lb_and_mom',obs,bf,scl,template,labels,prop,cluster{1},cluster{2},miss,part,dm,dat.mrf,ix_tiny,{'Template',oTemplate});                      
    
    % Check new lower bound
    if (dlb.Z + aff_reg) > (dat.lb.Z(end) + dat.lb.aff_reg(end))
        armijo = min(armijo*1.25,1);
        
        dat.lb.Z(end + 1)       = dlb.Z;                                                      
        dat.lb.aff_reg(end + 1) = aff_reg;   

        break;
    else
        armijo = armijo*0.5;
        r      = or;
        
        if line_search == nline_search                  
              E      = spm_dexpm(r,B);
              Affine = (mat_a\E*mat_s)*subsamp.MT;
    
              template = oTemplate; 
        end
    end
end
clear oTemplate

[dat.lb,gain] = check_convergence('aff',dat.lb,1,verbose,armijo);

% Set output
dat.reg.r      = r;
dat.armijo.aff = armijo;
%==========================================================================