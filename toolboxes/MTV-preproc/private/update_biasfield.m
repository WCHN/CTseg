function [Nii,ll1,ll3] = update_biasfield(Nii,dat,tau,num_workers,p)
% Update bias field coefficients (stored in Nii.b) using a Gauss-Newton
% step.
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Some parameters from options struct
modality      = p.Results.Modality;
speak         = p.Results.Verbose; 
method        = p.Results.Method;
EstimateRigid = p.Results.EstimateRigid;
bfreg         = p.Results.BiasFieldReg;
use_projmat   = ~(strcmpi(method,'denoise') && ~EstimateRigid);
C             = numel(dat); 

%--------------------------------------------------------------------------
% Start updating
%--------------------------------------------------------------------------

Nii_x = Nii.x;
Nii_y = Nii.y;
Nii_b = Nii.b;

ll1   = zeros(1,C);
ll3   = zeros(1,C);
% for c=1:C, fprintf('OBS! for c=1:C\n')
parfor (c=1:C,num_workers) % Loop over channels
            
    set_boundary_conditions;
    
    [Nii_b{c},ll1(c),ll3(c)] = update_channel(Nii_x{c},Nii_y(c),Nii_b{c},dat(c),tau{c},use_projmat,bfreg);    
end % End loop over channels
clear Nii_x Nii_y

Nii.b = Nii_b;
clear Nii_b

if speak >= 2
    show_model('bf',dat,Nii,modality);   
end
%==========================================================================

%==========================================================================
function [Nii_b,ll1,ll3] = update_channel(Nii_x,Nii_y,Nii_b,dat,tau,use_projmat,reg) 

N = numel(Nii_x);
y = get_nii(Nii_y); 

ll3 = zeros(1,N);
for n=1:N % Loop over number of observations (N) of channel c
    
    % Get voxel size
    vx = sqrt(sum(dat.A(n).mat(1:3,1:3).^2));
    
    % Get image and bias field coefficients
    x         = get_nii(Nii_x(n));
    msk       = get_msk(x);
    bfc       = get_nii(Nii_b(n));
    bfy       = exp(bfc);
    bfy(~msk) = 1;
    bfy       = bfy.*y;
    
    % Compute gradient
    gr       = tau(n)*bfy.*(bfy - x);
    gr(~msk) = 0;
    gr       = gr + spm_field('vel2mom', bfc, [vx 0 0 reg]);
        
    % Compute Hessian
    H        = tau(n)*(bfy.^2);
    H        = H + eps('single')*max(H(:));
    H(~msk)  = 0;
    clear bfy
        
    % Compute GN step
    Update = spm_field(H,gr,[vx 0 0 reg 2 2]);                          
    clear H gr
    
    % Do GN step
    bfc = bfc - Update;
    clear Update
    
    % Compute prior log-likelihood of one observation
    tmp    = spm_field('vel2mom', bfc, [vx 0 0 reg]);    
    ll3(n) = -0.5*double(bfc(:))'*double(tmp(:));

    % Save new bfc
    Nii_b(n) = put_nii(Nii_b(n),bfc);
    clear bfc
end

% Compute new log-likelihood
ll1 = get_ll1(use_projmat,true,y,Nii_x,Nii_b,tau,dat);
ll3 = sum(ll3);
%==========================================================================