function [model,dat] = update_template(dat,model,opt,is_init)
% FORMAT [model,dat] = update_template(dat,model,opt,is_init)
% dat     - Subject data structure
% model   - Model structure
% opt     - Options structure
% is_init - 
%
% Update the template by Gauss-Newton.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 4, is_init = false; end
if is_init
    opt.template.niter = 1;             % Number of GN iterations
end

% Options
reg        = opt.template.reg;          % [a m b] Template regularisation
verbose    = opt.template.verbose;      % Verbosity level
shrink     = opt.template.shrink;       % Crop template so it fits pushed data
B          = opt.reg.B;                 % Affine Lie basis
prm_reg    = opt.reg.rparam;            % [a m b le1 le1] Velocity regularisation
int_args   = opt.reg.int_args;          % Number of integration steps
load_a_der = opt.template.load_a_der;   % Precomputed template derivatives
R          = opt.template.R;            % Template null-space rotation
nii_a      = model.template.nii;        % Nifti holding the template

% Parameters
S0     = numel(dat);                    % Number of subjects
dm     = [nii_a.dat.dim,1,1,1];        
K      = dm(4);                         % Number of template classes
rits   = [3 3];                         % FMG cycles and relaxation iterations
mat_a  = nii_a.mat;                     % Template vox2world matrix
dm_a   = nii_a.dat.dim;                 % Template dimensions
vs_a   = spm_misc('vxsize',mat_a);      % Template voxel size
prm_a  = [vs_a prod(vs_a)*reg];         % Template reg param vector

% Initial starting estimates
a = rotate_template(nii_a,opt);         % Send previous template in null space

%--------------------------------------------------------------------------
% Compute/load pushed responsibilities then update template using Gauss-Newton
%--------------------------------------------------------------------------

if load_a_der && ~is_init
    % Template derivatives are loaded from disk
    %----------------------------------------------------------------------
                    
    % Add to global derivatives using bounding boxes
    H  = zeros([dm(1:3) round(((K-1)*K)/2)],'single'); % 2nd derivatives
    gr = zeros([dm(1:3),K-1],'single');                % 1st derivatives
    ll = 0;
    for s=1:S0
        % Add up gradients
        nii = nifti(dat{s}.template.pth_gr);
        gr  = gr + nii.dat(:,:,:,:);
        
        % Add up Hessians
        nii = nifti(dat{s}.template.pth_H);
        H   = H + nii.dat(:,:,:,:);
        
        % Add up ll:s
        ll = ll + dat{s}.template.ll;
    end            
    
    % Do Gauss-Newton step to update template
    [a,~,ll1,ss2] = solve_template(a,gr,H,prm_a,rits);
    
    if verbose>=2
        fprintf('update_template | %2d | ll = %6.6f    ll_a = %6.6f    ll + ll_a = %6.6f    ss2 = %6.6f\n',1,ll/prod(dm(1:3)), ll1/prod(dm(1:3)), (ll+ll1)/prod(dm(1:3)), (ss2)/prod(dm(1:3)));
    end
else
    % Template derivatives are computed on the fly.
    % This allows for iterating the Gauss-Newton solver.
    %----------------------------------------------------------------------
    
    for iter=1:opt.template.niter
        
        H  = zeros([dm(1:3) round(((K-1)*K)/2)],'single'); % 2nd derivatives
        gr = zeros([dm(1:3),K-1],'single');                % 1st derivatives    
        ll = 0;                                            % Template log-likelihood
        
        parfor s=1:S0        
%         for s=1:S0, fprintf('obs! for s=1:S0\n')
                         
            if is_init
                samp = 0;
            else
                samp = opt.seg.samp;
            end

            % Subject parameters
            [obs,dm_s,mat_s,vs_s,scl,~,~,~,~,nam,subsmp,grd] = get_obs(dat{s},'mskonlynan',opt.seg.mskonlynan,'samp',samp);                            
            labels                                           = get_labels(dat{s},opt,samp,subsmp,grd);
            miss                                             = get_par('missing_struct',obs);
            grd                                              = [];            
            
            % Misc parameters            
            ff        = get_ff(vs_s);     
            prop      = dat{s}.gmm.prop;
            if is_init
                prm_v = [vs_s ff*prm_reg];
            else
                prm_v = [subsmp.sk.*vs_s ff*prm_reg];                
            end
            modality  = dat{s}.modality{1}.name; 
            
            % Get initial velocity field
            if is_init
                v = zeros([dm_s 3],'single');
            else
                if isnumeric(dat{s}.reg.v)
                    % Initial velocity stored in array
                    v = dat{s}.reg.v;              
                else
                    % Initial velocity read from NIfTI
                    v = single(dat{s}.reg.v.dat(:,:,:,:));       
                end 
            end            
            
            % Get bias field
            do_bf = opt.bf.do;% && strcmpi(modality,'MRI');
            if do_bf
                if is_init
                    [x1,y1] = ndgrid(1:dm_s(1),1:dm_s(2),1);
                    z1      = 1:dm_s(3);
                    dat_tmp = dat{s};
                    for c=1:numel(dat_tmp.bf.chan)
                        d3                    = [size(dat_tmp.bf.chan(c).T) 1];
                        dat_tmp.bf.chan(c).B3 = spm_dctmtx(dm_s(3),d3(3),z1);
                        dat_tmp.bf.chan(c).B2 = spm_dctmtx(dm_s(2),d3(2),y1(1,:)');
                        dat_tmp.bf.chan(c).B1 = spm_dctmtx(dm_s(1),d3(1),x1(:,1));
                    end        
                    bf      = get_bf(dat_tmp.bf.chan,dm_s);
                    dat_tmp = []; x1 = []; y1 = []; z1 = [];                    
                else
                    bf = get_bf(dat{s}.bf.chan,dm_s);
                end
            else     
                bf = 1;
            end                        

            % Build and apply FFT of Green's function to map from momentum
            % to velocity
            if int_args > 1, Greens = spm_shoot_greens('kernel',dm_s(1:3),prm_v);
            else             Greens = [];
            end

            % Generates deformations from initial velocity fields by
            % gedesic shooting (if opt.reg.int_args > 1)
            y = make_deformation(v,prm_v,int_args,Greens);
            
            % Compute affine transformation matrix
            E      = spm_dexpm(dat{s}.reg.r,B); % Compute matrix exponential
            Affine = (mat_a\E*mat_s);   

            % Warp template to subject and softmax       
            [Template,y] = warp_template(model,y,Affine);              

            % Get responsibilities
            Z        = get_resp(obs,bf,dat{s},Template,labels,scl,miss,dm_s,opt);   
            % figure; imshow3D(squeeze(reshape(Z,[dm_s K])))                                          

            if opt.verbose.model >= 3
                % Write some results to disk (for visualisation in SegModel.m)
                dat{s} = write_for_visualisation(dm_s,obs,bf,dat{s},Template,v,labels,scl,miss,nam,prm_v,opt);
            end      
            Greens   = [];
            v        = [];
            labels   = []; 
            Template = []; 
            bf       = [];
            miss     = [];
            
            % Push responsibilities in subject space to template space
            Z = push_responsibilities(Z,y,dm_a(1:3)); 
            y = []; 
            % figure; imshow3D(squeeze(Z))
    
            % Compute gradients and Hessian
            [gr_s,H_s,ll_s] = diff_template(a,Z,prop,opt,is_init); 
            Z               = [];            

            % Add to global derivatives using bounding box         
            gr = gr + gr_s;
            H  = H  + H_s;            
            ll = ll + ll_s;                        
        end    

        % Do Gauss-Newton step to update template
        [a,ss1,ll1,ss2] = solve_template(a,gr,H,prm_a,rits);

        if verbose>=2
            fprintf('update_template | %2d | ll = %6.6f    ll_a = %6.6f    ll + ll_a = %6.6f    ss2 = %6.6f\n', iter,ll/prod(dm(1:3)), ll1/prod(dm(1:3)), (ll+ll1)/prod(dm(1:3)), (ss2)/prod(dm(1:3)));
        end

        if ss2/ss1<1e-4 
            % Converged?
            break; 
        end        
    end
end
ll = -(ll + ll1);

if ~isfinite(ll)
    % Just a simple sanity check
    error('~isfinite(ll)');
end

model.template.objval.post(end + 1)  = -(ll + ll1);
model.template.objval.likel(end + 1) = -ll;
model.template.objval.pr(end + 1)    = -ll1;

% Rotate back null-space (soft-maxed version is just used for
% visualisation)
[~,a] = softmax_template(a,R);

% Update the NIfTI that stores the template
nii_a.dat(:,:,:,:) = a;

if shrink
    % Use all the computed bounding-boxes from the push operations to
    % shrink the template.
    nii_a = spm_impreproc('mult_bb_crop',nii_a,BB);
end

model.template.nii = nii_a;

% Update values to be used for FOV voxels when warping 
model = init_template_bg(model,opt);

if verbose>=1
    % Show soft-maxed template
    show_template(model,opt,S0);
end
%==========================================================================

%==========================================================================
function [a,ss1,ll1,ss2] = solve_template(a,gr,H,prm,its)
% FORMAT [a,ss1,ll1,ss2] = solve_template(a,gr,H,prm,its)
% a   - log-template in null space
% gr  - Template gradient
% H   - Template Hessian
% prm - Regularisation parameters  [vs abs memb bend]
% its - Full multi grid parameters [cycles relax] 
%
% Perform one Gauss-Newton update of the template.
% ss1 and ss2 are for examining how close the 1st derivatives are to zero.
% At convergence, the derivatives from the likelihood term should match those
% from the prior (regularisation) term.
K   = size(a,4);
ss1 = sum(sum(sum(sum(gr.^2))));
gr1 = spm_field('vel2mom',a,prm);     % 1st derivative of the prior term
ll1 = 0;
for k=1:K
    ll1 = ll1 + 0.5*sum(sum(sum(sum(double(gr1(:,:,:,k)).*double(a(:,:,:,k)))))); % -ve log probability of the prior term
end
gr  = gr + gr1;                       % Combine the derivatives of the two terms
clear gr1

if ~isfinite(ll1)
    warning('~isfinite(ll1)');
end
            
ss2 = sum(sum(sum(sum(gr.^2))));      % This should approach zero at convergence

a = a - spm_field(H,gr,[prm(1:3) prm(4) prm(5:6) its]); % Gauss-Newton update  
%==========================================================================