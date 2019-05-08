function [gr,H,ll] = diff_template(a,Z,prop,opt,is_init)
% FORMAT [gr,H,ll] = diff_template(a,Z,prop,opt)
% a    - Current log-template warped in subject-space [Nx Ny Nz K]
% Z    - Class responsibilities [Nx Ny Nz K]
% prop - Tissue proportions
% opt  - Options structure
% gr   - Gradient [Nx Ny Nz K]
% H    - Hessian [Nx Ny Nz K(K+1)/2]
% ll   - Log-likelihood
%
% Compute derivatives of the log-template for a given subject.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 5, is_init = false; end

% Parameters
dm      = [size(Z),1,1,1];
K       = dm(4);
ln_prop = reshape(prop,1,1,K);
R       = null(ones(1,K));

if opt.template.sym
    % Make all classes left-right symmetric
    Z = Z + Z(end:-1:1,:,:,:);
end

% Preproc sufficient statistics (here Z)
smZ = sum(Z,4);
msk = smZ > 0;

if ~isempty(opt.template.tc_miss) && ~is_init
    % Set missing data as being from a specific class    
    valZ                            = mean(smZ(:));
    Z                               = reshape(Z,[prod(dm(1:3)) K]);
    Z(~msk(:),opt.template.tc_miss) = valZ;    
    Z                               = reshape(Z,[dm(1:3) K]);
    
    smZ = sum(Z,4);
    msk = msk >= 0; % Do not mask anymore..
end

msk = reshape(msk,[prod(dm(1:2)) dm(3)]);
Z   = bsxfun(@rdivide,Z,smZ + eps('single'));
mxZ = max(smZ(:));

%--------------------------------------------------------------------------
% Compute gradients and Hessian
%--------------------------------------------------------------------------

gr  = zeros([dm(1:3),K-1],'single');                % 1st derivatives
H   = zeros([dm(1:3) round(((K-1)*K)/2)],'single'); % 2nd derivatives
dgr = size(gr);
ll  = 0;
for z=1:size(Z,3) % Loop over planes

    % Compute softmax for this plane
    mu = softmax_template(a(:,:,z,:),R,ln_prop);
    
    % Mask where sum(Z) == 0
    mu = reshape(mu,[prod(dm(1:2)) K]);
    for k=1:K        
        mu(~msk(:,z),k) = NaN;        
    end
    
    % Reshape and double
    mu = double(reshape(mu,[dm(1:2) K]));
    
    % Compupte log-likelihood
    ll = ll - nansum(nansum(nansum(log(mu + eps).*reshape(Z(:,:,z,:),[dm(1:2),dm(4)]),3).*smZ(:,:,z)));
            
    if ~isfinite(ll)
        warning('~isfinite(ll)');
    end
            
    % Compute first derivatives (d(4)-1) x 1 
    grz = mu - double(reshape(Z(:,:,z,:),[dm(1:2),K]));
    for j1=1:(K-1)
        gr1 = zeros([dm(1:2) 1 dgr(4)],'single');
        for j2=1:K
            gr1(:,:,1,j1) = gr1(:,:,1,j1) + R(j2,j1)*grz(:,:,j2); % Note the rotation
        end
        gr(:,:,z,j1) = gr(:,:,z,j1) + gr1(:,:,1,j1).*smZ(:,:,z);
    end

    % Compute d(4) x d(4) matrix of second derivatives at each voxel.
    % These should be positive definate, but rounding errors may prevent this.
    % Regularisation is included to enforce +ve definateness.
    hz = zeros([dm(1:2),K,K]);
    for j1=1:K
        hz(:,:,j1,j1) =   (1 - mu(:,:,j1)).*mu(:,:,j1).*smZ(:,:,z);
        for j2=1:(j1-1)
            hz(:,:,j1,j2) = -mu(:,:,j1) .*mu(:,:,j2).*smZ(:,:,z);
            hz(:,:,j2,j1) = hz(:,:,j1,j2);
        end
    end

    % First step of rotating 2nd derivatives to (d(4)-1) x (d(4)-1)
    % by R'*W*R
    hz1 = zeros([dm(1:2),K,K-1]);
    for j1=1:K
        for j2=1:(K-1)
            tmp = zeros(dm(1:2));
            for j3=1:K
                tmp = tmp + hz(:,:,j1,j3)*R(j3,j2);
            end
            hz1(:,:,j1,j2) = tmp;
        end
    end

    % Second step of rotating 2nd derivatives to (d(4)-1) x (d(4)-1)
    % by R'*W*R
    hz = zeros([dm(1:2),K-1,K-1]);
    for j1=1:(K-1)
        for j2=1:(K-1)
            tmp = zeros(dm(1:2));
            for j3=1:K
                tmp = tmp + R(j3,j1)*hz1(:,:,j3,j2);
            end
            hz(:,:,j1,j2) = tmp;
        end
    end

    % First pull out the diagonal of the 2nd derivs
    for j1=1:K-1
        H(:,:,z,j1) = H(:,:,z,j1) + hz(:,:,j1,j1) + mxZ*sqrt(eps('single'))*K^2;
    end

    % Then pull out the off diagonal parts (note that matrices are symmetric)
    jj = K;
    for j1=1:K-1
       for j2=(j1+1):(K-1)
           H(:,:,z,jj) = H(:,:,z,jj) + hz(:,:,j2,j1);
           jj = jj+1;
       end
    end
end

gr(~isfinite(gr)) = 0;
H(~isfinite(H))   = 0;
%==========================================================================