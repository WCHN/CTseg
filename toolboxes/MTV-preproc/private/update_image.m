function [Nii,ll1,ll2,mtv_scale]= update_image(Nii,dat,tau,rho,lam,num_workers,p)
% Update Nii.y, Nii.u, Nii.w by an ADMM algorithm
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Some parameters from options struct
modality      = p.Results.Modality;
method        = p.Results.Method;
nitgn         = p.Results.IterGaussNewtonImage; 
speak         = p.Results.Verbose; 
EstimateRigid = p.Results.EstimateRigid;
EstimateBias  = p.Results.EstimateBias;
ApplyBias     = p.Results.ApplyBias;

% Flag saying if we solve using projection matrices (A, At), or not
use_projmat = ~(strcmpi(method,'denoise') && ~EstimateRigid);
inc_bf      = EstimateBias || ApplyBias;

% Get data from Nii struct (otherwise we get parfor errors)
Nii_x = Nii.x;
Nii_y = Nii.y;
Nii_H = Nii.H;
Nii_w = Nii.w;
Nii_u = Nii.u;
Nii_b = Nii.b;

C  = numel(Nii_x);
vx = sqrt(sum(dat(1).mat(1:3,1:3).^2));
dm = dat(1).dm;

%--------------------------------------------------------------------------
% First update y
%--------------------------------------------------------------------------

ll1 = zeros(1,C);
% for c=1:C, fprintf('OBS! for c=1:C\n')
parfor (c=1:C,num_workers) % Loop over channels

    set_boundary_conditions;

    u = get_nii(Nii_u(c));   
    w = get_nii(Nii_w(c));   
    y = get_nii(Nii_y(c));   
    
    if use_projmat
        % We use the projection matrices (A, At)
        
        for gnit=1:nitgn % Iterate Gauss-Newton

            % Gradient      
            rhs = w/rho - u; 
            rhs = lam(c)*imdiv(rhs,vx);            
            for n=1:dat(c).N
                % Here we discard missing data, for MRI these are
                % assumed to be zeros and NaNs.
                x         = get_nii(Nii_x{c}(n));
                Ayx       = A(y,dat(c),n);
                msk       = get_msk(x,Ayx);
                Ayx       = Ayx - x;
                Ayx(~msk) = 0;
                
                rhs = rhs + At(Ayx,dat(c),tau{c},n)*(1/rho); 
            end                  
            msk = [];
            x   = [];                        
            Ayx = [];
            
            rhs = rhs + spm_field('vel2mom',y,[vx 0 lam(c)^2 0]);

            % Hessian
            H   = get_nii(Nii_H(c));
            lhs = H*sum(tau{c})/rho;
            H   = [];

            % Compute GN step
            y   = y - spm_field(lhs,rhs,[vx 0 lam(c)^2 0 2 2]);
            lhs = [];
            rhs = [];
        end                                
    else
        % We do not use the projection matrices (A, At)
        
        % RHS
        rhs = u - w/rho; 
        rhs = lam(c)*imdiv(rhs,vx);
        for n=1:dat(c).N
            
            if inc_bf
                bf = exp(get_nii(Nii_b{c}(n)));
            else
                bf = 1;
            end
    
            x         = get_nii(Nii_x{c}(n));
            msk       = get_msk(x);
            tmp       = (bf.*x)*(tau{c}(n)/rho);
            tmp(~msk) = 0;
            rhs       = rhs + tmp;
        end
        x   = [];
        tmp = [];
        msk = [];
        bf  = [];
                
        if inc_bf && ~use_projmat
            bf = zeros([dm dat(c).N],'single');
            for n=1:dat(c).N
                bf(:,:,:,n) = exp(get_nii(Nii_b{c}(n)));
            end
        else
            bf = ones(dm,'single');
        end
            
        % LHS
        lhs = sum(reshape(tau{c},[1 1 1 dat(c).N]).*bf.^2,4)/rho;

        % Compute new y
        y   = spm_field(lhs,rhs,[vx 0 lam(c)^2 0 2 2]);
        lhs = [];
        rhs = [];        
    end    
    
    if strcmpi(modality,'MRI')
        % Ensure non-negativity (ad-hoc)
        y(y < 0) = 0;
    end    
    
    Nii_y(c) = put_nii(Nii_y(c),y);
    
    % Compute log of likelihood    
    ll1(c) = get_ll1(use_projmat,inc_bf,y,Nii_x{c},Nii_b{c},tau{c},dat(c));
    y      = [];    
    x      = [];
    
end % End loop over channels     

%--------------------------------------------------------------------------
% Then compute MTV
%--------------------------------------------------------------------------

unorm = 0;    
% for c=1:C, fprintf('OBS! for c=1:C\n')
parfor (c=1:C,num_workers) % Loop over channelsi

    set_boundary_conditions;
    
    y = get_nii(Nii_y(c));        
    G = lam(c)*imgrad(y,vx);
    y = [];              
    
    w = get_nii(Nii_w(c));        
    u = G + w/rho;    
    w = [];
    G = [];
    
    Nii_u(c) = put_nii(Nii_u(c),u);

    unorm = unorm + sum(sum(u.^2,4),5);
    u     = [];
    
end % End loop over channels

unorm     = sqrt(unorm);
mtv_scale = max(unorm - 1/rho,0)./(unorm + eps);
clear unorm

%--------------------------------------------------------------------------
% Update u and w
%--------------------------------------------------------------------------

ll2 = 0;
u2  = single(0);
% for c=1:C, fprintf('OBS! for c=1:C\n')
parfor (c=1:C,num_workers) % Loop over channels

    set_boundary_conditions;

    u = get_nii(Nii_u(c));   
    w = get_nii(Nii_w(c));       
    
    %------------------------------------------------------------------
    % Update proximal operator for u
    % Here we multiply each contrast image with the same scaling
    % matrix, this is a key addition of using MTV
    %------------------------------------------------------------------

    u        = bsxfun(@times,u,mtv_scale);
    Nii_u(c) = put_nii(Nii_u(c),u);
    u2       = u2 + sum(sum(u.^2,4),5);
    
    %------------------------------------------------------------------
    % Solve for w
    % Here we update the Lagrange variable
    %------------------------------------------------------------------
    
    y = get_nii(Nii_y(c));
    G = lam(c)*imgrad(y,vx);
    a = u - G;
    u = [];
    G = [];
    
    w        = w - rho*a;  
    Nii_w(c) = put_nii(Nii_w(c),w);
                   
    %------------------------------------------------------------------
    % Update augmented lagrangian part of the log-likelihood
    %------------------------------------------------------------------

    ll2  = ll2 - a(:)'*(0.5*rho*a(:) - w(:));
    a    = [];
    w    = [];
    
end % End loop over channels     

% Update augmented lagrangian part of the prior
ll2 = ll2 - sum(sum(sum(sqrt(double(u2)))));
clear u2

Nii.y = Nii_y;
Nii.w = Nii_w;
Nii.u = Nii_u;

if speak >= 3
    % Show MTV prior
    show_model('solution',use_projmat,modality,Nii,dat);
    show_model('mtv',mtv_scale);
%     show_model('rgb',Nii_y);
end
%==========================================================================