function varargout = gmm_mrf(varargin)
%__________________________________________________________________________
% MRF functions for images (2D and 3D).
%
% FORMAT help gmm_mrf>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help gmm_mrf
    error('Not enough argument. Type ''help gmm_mrf'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'apply'
        [varargout{1:nargout}] = mrf_apply(varargin{:});            
    case 'lowerbound'
        [varargout{1:nargout}] = mrf_lb(varargin{:});               
    case 'update'
        [varargout{1:nargout}] = mrf_update(varargin{:});                       
    otherwise
        help gmm_mrf
        error('Unknown function %s. Type ''help gmm_mrf'' for help.', id)
end
%==========================================================================

%==========================================================================
function lnPzN = mrf_apply(mrf,z)
if nargin < 2, z = 0; end

K = size(mrf.G,1);
if mrf.do && isfield(mrf,'oZ')
    % Compute log-MRF part of responsibilities   
    lnG = single(log(mrf.G));
    w   = single(mrf.w);
    dm  = size(mrf.oZ);
    % figure; imshow3D(squeeze(mrf.oZ))
%     lnPzN  = compute_mrf(mrf.oZ,lnG,w);

    if z > 0 && z <= dm(3) 
        lnPzN = mex_mrf('apply',mrf.oZ(:,:,z,:),lnG,w);
                
        lnPzN = reshape(lnPzN,[prod(dm(1:2)) K]); % Reshape
    else
        lnPzN = mex_mrf('apply',mrf.oZ,lnG,w);
        
        lnPzN = reshape(lnPzN,[prod(dm(1:3)) K]); % Reshape
    end
    % figure; imshow3D(squeeze(lnPzN))        
else
    lnPzN = zeros([1 K]);
end
%==========================================================================

%==========================================================================
function lbZN = mrf_lb(mrf)
if mrf.do && isfield(mrf,'oZ')
    lnG  = single(log(mrf.G));
    w    = single(mrf.w);
    
%     lbZN = compute_mrf_lb(mrf.oZ,lnG,w);

    lbZN = mex_mrf('lowerbound',mrf.oZ,lnG,w);
else
    lbZN = 0;
end

return
%==========================================================================

%==========================================================================
function mrf = mrf_update(mrf)
if mrf.do && isfield(mrf,'oZ')
    w    = double(mrf.w);
    tiny = eps('single');
    
%     G = compute_mrf_update(mrf.oZ,w);
%     G = max(G,eps('single'));
%     G = bsxfun(@rdivide,G,sum(G,2));
    
    G = mex_mrf('update',mrf.oZ,single(w));
    
    % Normalise
    G  = max(G,tiny);
    sG = sum(G,2);
    G  = bsxfun(@rdivide,G,sG(:)');    
    
    mrf.G = G; 
end
%==========================================================================

%==========================================================================
% CORE FUNCTIONS
%==========================================================================

%==========================================================================
function lnPzN = compute_mrf(Z,lnG,w)
dm = size(Z);
K  = dm(4);

lnPzN = zeros(dm,'single');
for z=1:dm(3) % loop over dm(3)
    for y=1:dm(2) % loop over dm(2)
        for x=1:dm(1) % loop over dm(1)
            
            N = get_neighborhood(Z,x,y,z,dm,w); % Get neighborhood voxels

            for k=1:K % loop over tissue classes
                sN = 0;
                for l=1:K % loop over tissue classes
                    for j=1:6 % loop over neighborhood
                        sN = sN + N(j,l)*lnG(k,l);
                    end
                end

                lnPzN(x,y,z,k) = sN;
            end
        end
    end
end
%==========================================================================

%==========================================================================
function lbZN = compute_mrf_lb(Z,lnG,w)
dm = size(Z);
K  = dm(4);

% Checkerboard loop (only over one color)
lnPzN = zeros(K,1);
lbZN  = 0;   
for z=1:dm(3) % loop over dm(3)
    y1 = (1 == mod(z,2));
    for y=1:dm(2) % loop over dm(2)
        x1 = (y1 == mod(y,2)) + 1;
        for x=x1:2:dm(1) % loop over dm(1)                

            N = get_neighborhood(Z,x,y,z,dm,w); % Get neighborhood voxels

            for k=1:K % loop over tissue classes
                sN = 0;
                for l=1:K % loop over tissue classes        
                    for j=1:6 % loop over neighborhood
                        sN = sN + N(j,l)*lnG(k,l);
                    end
                end

                lnPzN(k) = sN;
            end

            Z1   = Z(x,y,z,:);            % Get responsibility for voxel Z(x,y,z)
            Z1   = single(Z1(:))./255;    % Convert to single
            lbZN = lbZN + sum(Z1.*lnPzN); % Add to lower bound

        end
    end
end
%==========================================================================

%==========================================================================
function G = compute_mrf_update(Z,w)
dm = size(Z);
K  = dm(4);
G  = zeros(K,'double');
for k=1:K % loop over tissue classes
    for l=1:K % loop over tissue classes
                
        % Checkerboard loop (only over one color)
        g = 0;
        for z=1:dm(3) % loop over dm(3)
            y1 = (1 == mod(z,2));
            for y=1:dm(2) % loop over dm(2)
                x1 = (y1 == mod(y,2)) + 1;
                for x=x1:2:dm(1) % loop over dm(1)                

                    N = get_neighborhood(Z,x,y,z,dm,w); % Get neighborhood voxels

                    Z1 = double(Z(x,y,z,k))./255; % Convert to single
                    g  = g + Z1*sum(double(N(:,l)));
                end
            end
        end

        G(k,l) = g; % Update (k,l) element of G
    end
end
%==========================================================================

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

%==========================================================================
function N = get_neighborhood(Z,x,y,z,dm,w)
N = zeros([6 dm(4)],'single');
if x > 1,     N(1,:) = single(Z(x - 1,y,z,:))./255; end
if x < dm(1), N(2,:) = single(Z(x + 1,y,z,:))./255; end
if y > 1,     N(3,:) = single(Z(x,y - 1,z,:))./255; end
if y < dm(2), N(4,:) = single(Z(x,y + 1,z,:))./255; end
if z > 1,     N(5,:) = single(Z(x,y,z - 1,:))./255; end
if z < dm(3), N(6,:) = single(Z(x,y,z + 1,:))./255; end
%==========================================================================