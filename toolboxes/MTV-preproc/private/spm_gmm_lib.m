function varargout = spm_gmm_lib(varargin)
%__________________________________________________________________________
%
% Library of functions for Gaussian Mixture modelling
%
%--------------------------------------------------------------------------
% Convention
% ----------
% N  - Number of observations
% Nm - Number of observations that correspond to a given "missing code"
% P  - Dimension of observations
% Po - Number of non-missing observations in a given "missing code"
% Pm - Number of     missing observations in a given "missing code"
% K  - Number of clusters
% M  - Number of unique "missing codes"
%
% X   -  N  x P        Observations                    (+ inferred missing)
% Z   -  N  x K        Clusters' responsibility
% PI  - [N] x K        Clusters' proportions
% a   -  1  x K        Proportions concentration       (Dirichlet prior)
% W   - [N] x 1        Observations weights            (Histogram GMM) 
% MU  -  P  x K        Clusters' [expected] mean
% A   -  P  x P  x K   Clusters' [expected] precision matrix
% b   -  1  x K        Mean degrees of freedom         (Gauss prior)
% V   -  P  x P  x K   Precision scale matrix          (Wishart prior)
% n   -  1  x K        Precision degrees of freedom    (Wishart prior)
% E   - [N] x P        Bin variance                    (Histogram GMM) 
% SS0 -  1  x K        Zeroth order sufficient statistics
% SS1 -  P  x K        First  order sufficient statistics
% SS2 -  P  x P  x K   Second order sufficient statistics
% C   -  N  x 1        Code image: one code / missing combination
% L   -  1  x M        List of unique codes (saves time)
% 
%--------------------------------------------------------------------------
% Update functions
% ----------------
%
% logp = spm_gmm_lib('Marginal', X, {MU,A},   const, {C, L}, E)
% logp = spm_gmm_lib('Marginal', X, {MU,V,n}, const, {C, L}, E)
% > Observation's marginal log probability within each cluster
%
% const = spm_gmm_lib('Const',  MU,     A,    (L))
% const = spm_gmm_lib('Const', {MU,b}, {V,n}, (L))
% > Constant term (normalisation) of a Gaussian log-distribution
%   If L is provided -> marginal distributions
%
% Z = spm_gmm_lib('Responsibility', logpX, logPI)
% > Compute & normalise responsibilities (safe softmax)
%
% FORMAT [SS0,SS1,SS2] = spm_gmm_lib('SuffStat', X, Z, W)
% FORMAT [SS0,SS1,SS2] = spm_gmm_lib('SuffStat', 'base',  X, Z, W)
% FORMAT [SS0,SS1,SS2] = spm_gmm_lib('SuffStat', 'infer', SS0, SS1, SS2, {MU,A}, L)
% FORMAT         [SS2] = spm_gmm_lib('SuffStat', 'bin',   E, Z, W, {C,L})
% > Compute sufficient statistics (0th, 1st, 2nd order)
%   default: use only non-missing data => E[z], E[z]*x, E[z]*xx'
%   'base':  statistics / each config  => E[z], E[z]*g, E[z]*gg'
%   'infer': base to full statistics   => E[z], E[z*x], E[z*xx']
%   'bin':   binning uncertainty       => Tr(S\cov[gg'])
%
% [MU,A]       = spm_gmm_lib('UpdateClusters', SS0, SS1, SS2)
% [MU,A,b,V,n] = spm_gmm_lib('UpdateClusters', SS0, SS1, SS2, {MU0,b0,V0,n0})
% > Update GMM parameters (ML or Bayesian posterior)
%
% [PI,logPI,a] = spm_gmm_lib('UpdateProportions', SS0, a0)
% > Update cluster proportions
%
% [GaussPrior,extras] = spm_gmm_lib('updatehyperpars',cluster,GaussPrior,varargin)
% > Update VB-GMM hyper-parameters (MU,b,V,n)
%
% X = spm_gmm_lib('InferMissing', X, Z, {MU,A}, {C,L})
% > Infer missing values (this should only be done once, at the end)
%
%--------------------------------------------------------------------------
% Lower bound functions
% ---------------------
%
% [lb,const] = spm_gmm_lib('MarginalSum', SS0, SS1, SS2,  MU,     A,    L, SS2b)
% [lb,const] = spm_gmm_lib('MarginalSum', SS0, SS1, SS2, {MU,b}, {V,n}, L, SS2b)
% > Compute conditional datasum: E[ln p(g|MU,A,Z)]
%   Also returns the result of spm_gmm_lib('const')
%
% [klMU,klA] = spm_gmm_lib('KL', 'GaussWishart', {MU,b}, {V,n}, {MU0,b0}, {V0,n0})
% > KL divergence between two Gauss-Wishart distributions
%
% klP = spm_gmm_lib('KL', 'Dirichlet', a, a0)
% > KL divergence between two Dirichlet distributions
%
% klZ = spm_gmm_lib('KL', 'Categorical', Z, W, logPI)
% > KL divergence between two Categorical distributions
%
%--------------------------------------------------------------------------
% "Missing code" image 
% --------------------
%
% C = spm_gmm_lib('obs2code', X)
% > Transform a vector of observations (with NaNs) into a code image
%
% B = spm_gmm_lib('code2bin', C, length)
% > Convert a given code to a binary mask of "non-missing" dimensions
%
% T = spm_gmm_lib('range2int', maxval, minval)
% > Compute the smallest integer type needed to store a given range
%
% L = spm_gmm_lib('double2int', L)
% > Convert an array of floats/doubles into an array of integers
%   (chooses the integer type of minimal length)
%
%--------------------------------------------------------------------------
% Visualisation 
% -------------
%
% spm_gmm_lib('Plot', 'LB', lb)
% > Plot lower bound
%
% spm_gmm_lib('Plot', 'GMM', {X,W}, {MU,A}, PI)
% > Plot mixture fit
%
% spm_gmm_lib('plot', 'cat', dm, Z, Template, (wintitle))
% > Plot (categorical) responsibilities and template (if available)
%
% spm_gmm_lib('plot', 'gaussprior', GaussPrior, (wintitle))
% > Plot VB-GMM hyper-parameters
%
% c = spm_gmm_lib('plot', 'cat2rgb', f, pal)
% > Generate an RGB volume from a categorical (e.g. responsibilities) volume.
%
%--------------------------------------------------------------------------
% Extras
% -------------
%
% gmm = spm_gmm_lib('extras', 'more_gmms', gmm, part)
% > A crude heuristic to replace a single Gaussian by a bunch of Gaussians.
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin == 0
    help spm_gmm_lib
    error('Not enough argument. Type ''help spm_gmm_lib'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'infermissing'
        [varargout{1:nargout}] = infermissing(varargin{:});
    case 'marginal'
        [varargout{1:nargout}] = marginal(varargin{:});
    case 'responsibility'
        [varargout{1:nargout}] = responsibility(varargin{:});
    case 'suffstat'
        [varargout{1:nargout}] = suffstat(varargin{:});
    case 'const'
        [varargout{1:nargout}] = const(varargin{:});
    case 'updateclusters'
        [varargout{1:nargout}] = updateclusters(varargin{:});   
    case 'updateproportions'
        [varargout{1:nargout}] = updateproportions(varargin{:});  
    case 'updatehyperpars'
        [varargout{1:nargout}] = updatehyperpars(varargin{:});           
    case 'marginalsum'
        [varargout{1:nargout}] = marginalsum(varargin{:});
    case 'kl'
        [varargout{1:nargout}] = kl(varargin{:}); 
    case 'obs2code'
        [varargout{1:nargout}] = obs2code(varargin{:}); 
    case 'code2bin'
        [varargout{1:nargout}] = code2bin(varargin{:}); 
    case 'range2int'
        [varargout{1:nargout}] = range2int(varargin{:}); 
    case 'double2int'
        [varargout{1:nargout}] = double2int(varargin{:});  
    case 'plot'
        [varargout{1:nargout}] = gmmplot(varargin{:});       
    case 'extras'
        [varargout{1:nargout}] = gmm_extras(varargin{:});             
    otherwise
        help spm_gmm_lib
        error('Unknown function %s. Type ''help spm_gmm_lib'' for help.', id)
end

%--------------------------------------------------------------------------
% Update functions
%--------------------------------------------------------------------------

% =========================================================================
function X = infermissing(X, Z, cluster, codes)
% FORMAT X = spm_gmm_lib('missing', X, Z, {MU,A}, {C,L})
% X  - NxP   observations
% Z  - NxK   responsibilities
% MU - PxK   (expected) means
% A  - PxPxK (expected) precision matrices
% C  - Nx1   "missing value" code image
% L  -       list of existing codes
%
% X - NxP    observations with inferred values
%
% Compute the mean expected value of missing voxels.

MU = [];
A  = [];
C  = [];
L  = [];

%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(cluster)
    MU = cluster;
else
    if numel(cluster) >= 1
        MU = cluster{1};
        if numel(cluster) >= 2
            A = cluster{2};
        end
    end
end
if nargin >= 4
    if ~iscell(codes)
        C = codes;
    else
        if numel(codes) >= 1
            C = codes{1};
            if numel(codes) >= 2
                L = codes{2};
            end
        end
    end
    if isempty(L)
        L = unique(C);
    end
end

%--------------------------------------------------------------------------
% Dimensions
P = size(X, 2);
K = size(Z, 2);
if isempty(L)
    L = 2^P - 1; % None missing
end

% -------------------------------------------------------------------------
% For each missing combination
for i=1:numel(L)
    
    % ---------------------------------------------------------------------
    % Get code / mask / missing modalities
    c        = L(i);
    observed = code2bin(c, P);
    missing  = ~observed;
    Pm       = sum(missing);
    if Pm == 0, continue; end
    msk      = (C == c);
    Nm       = sum(msk);
    if Nm == 0, continue; end
    
    % ---------------------------------------------------------------------
    % Initialise
    X(msk,missing) = 0;

    % ---------------------------------------------------------------------
    % Compute posterior mean (expected value)
    % 1) t = sum_k {z * ( mu[m] + A[m]/A[m,o]*(mu[o]-g) ) }
    for k=1:K
        X1k = zeros(1, 'like', X);
        X1k = bsxfun(@plus,X1k,MU(missing,k).');
        X1k = bsxfun(@plus,X1k,bsxfun(@minus, MU(observed,k).', X(msk,observed)) * (A(observed,missing,k) / A(missing,missing,k)));
        X(msk,missing) = X(msk,missing) + bsxfun(@times, X1k, Z(msk,k));
    end
                    
end

% =========================================================================
function logpX = marginal(X, cluster, const, codes, E)
% logp = spm_gmm_lib('marginal', X, {MU,A}, const, {C, L}, E)
% logp = spm_gmm_lib('marginal', X, {MU,V,n}, const, {C, L}, E)
% 
% X         - NxP   Observed values
% MU        - PxK   (Expected) means
% A         - PxPxK (Expected) precision matrices
% const     - MxK   Constant terms. If M > 1, marginal distributions.
% C         - Nx1   Image of "missing" codes
% L         - Mx1   List of existing codes
% E         - 1xP   Binning uncertainty
%
% logpX     - NxK   (Expected) log-likelihood of belonging to each class
%
% Compute the expected log-likelihood of each observation belonging to each
% cluster: logpx(i,k) = E[ln p(g(i) | MU_k,A_k)]


MU = [];
A  = [];
V  = [];
n  = [];
C  = [];
L  = [];

%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(cluster)
    MU = cluster;
else
    if numel(cluster) >= 1
        MU = cluster{1};
        if numel(cluster) >= 2
            A = cluster{2};
            if numel(cluster) >= 3
                V = A;
                n = cluster{3};
            end
        end
    end
end
if ~iscell(codes)
    C = codes;
else
    if numel(codes) >= 1
        C = codes{1};
        if numel(codes) >= 2
            L = codes{2};
        end
    end
end
if isempty(L)
    L = unique(C);
end
if nargin < 5
    E = [];
end

% -------------------------------------------------------------------------
% Use double precision
X  = double(X);
MU = double(MU);
A  = double(A);
V  = double(V);
n  = double(n);
E  = double(E);
const = double(const);

% -------------------------------------------------------------------------
% Dimensions
N  = size(X,1);
P  = size(X,2);
K  = size(MU,2);
if isempty(L), L = 2^P - 1;    end % None missing
if isempty(E), E = zeros(1,P, 'like', X); end % No uncertainty
if size(const, 1) == 1
    const = repmat(const, [numel(L) 1]);
end

logpX  = zeros([N K], 'like', X);

% -------------------------------------------------------------------------
% For each combination of missing voxels
for i=1:numel(L)
    
    % ---------------------------------------------------------------------
    % Get mask of missing values and modalities (with this particular code)
    c        = L(i);
    observed = code2bin(c, P);
    Po       = sum(observed);
    if Po == 0, continue; end
    if isempty(C) && (c == 2^P-1),  msk = ones(N,1, 'logical');
    else,                           msk = (C == c);   end
    missing  = ~observed;
    Pm       = P-Po;
    Nm       = sum(msk(:));
    if Nm == 0, continue; end
    
    % ---------------------------------------------------------------------
    % Initialise with constant term
    logpX(msk,:) = repmat(const(i,:), [Nm 1]);
    
    % ---------------------------------------------------------------------
    % Non constant terms
    X1 = X(msk,observed)';
    for k=1:K
        
        % /!\ Sub-covariance is different from the inverse sub-precision
        % ML case:
        %   inv(S(o,o)) = A(o,o) - A(o,m)*A(m,m)\A(m,o)
        % Bayesian case:
        %   inv(S(o,o)) ~ W(V(o,o) - V(o,m)*V(m,m)\V(m,o), n - Pm)
        if sum(n) > 0
%             Ao = V(observed,observed,k) - V(observed,missing,k)*(V(missing,missing,k)\V(missing,observed,k));
            Ao = V(observed,observed,k) - V(observed,missing,k)*(spm_matcomp('inv',V(missing,missing,k))*V(missing,observed,k));
            Ao = (n(k)-Pm) * Ao;
        else
%             Ao = A(observed,observed,k) - A(observed,missing,k)*(A(missing,missing,k)\A(missing,observed,k));
            Ao = A(observed,observed,k) - A(observed,missing,k)*(spm_matcomp('inv',A(missing,missing,k))*A(missing,observed,k));
        end
        
        % Quadratic term in observed values: (obs-mean) x (obs-mean)
        l = Ao * bsxfun(@minus, X1, 2*MU(observed,k));
        l = -0.5 * dot(l, X1, 1);
        
        % Binning uncertainty
        if any(any(E))
            if size(E,1)==1
                l = l - 0.5 * trace(Ao * diag(E(observed)));
            else
                l = l - 0.5 * sum(bsxfun(@times,diag(Ao)',E(msk,observed)),2)';
            end
        end
        
        % Reshape as a column vector
        logpX(msk,k) = logpX(msk,k) + l';
    end
        
end

% Set NaN value for voxels without observed dimensions
logpX(all(isnan(X),2),:) = NaN;

% =========================================================================
function Z = responsibility(logpX, logPI, varargin)
% FORMAT Z = spm_gmm_lib('responsibility', logpX, logPI)
%
% Compute responsibilities.
% Responsibilities are the posterior expected value of class-indexing 
% vectors z_n.
% The posterior is computed as:
%   r_nk = exp(E[log Pi_k] + E[log p(x_n | Theta_k)]) / sum_k {r_nk}
%
% Extra terms can be added to the responsibilities prior to softmax by the
% varargin argument. These arguments need to be compatible (w.r.t. size) with 
% the following function call: bsxfun(@plus, Z, varargin{i}).

% Use double precision
logpX = double(logpX);
logPI = double(logPI);


% Omit NaN
logpX(isnan(logpX)) = 0; 

% Add terms: E[log Pi_k] + E[log p(X | Theta_k)]
Z = bsxfun(@plus, logpX, logPI);

for i=1:numel(varargin)
    Z = bsxfun(@plus, Z, double(varargin{i}));
end

% Exponentiate and normalise
Z = bsxfun(@minus, Z, max(Z, [], 2));
Z = exp(Z);
Z = bsxfun(@rdivide, Z, sum(Z, 2));

% =========================================================================
function varargout = suffstat(varargin)
% FORMAT [SS0,SS1,SS2] = spm_gmm_lib('suffstat', X, Z, W)
% FORMAT [SS0,SS1,SS2] = spm_gmm_lib('suffstat', 'base',  X, Z, W, {C,L})
% FORMAT [SS0,SS1,SS2] = spm_gmm_lib('suffstat', 'infer', SS0, SS1, SS2, {MU,A}, L)
% FORMAT         [SS2] = spm_gmm_lib('suffstat', 'bin',   E, Z, W, {C,L})
%
% X    - NxP Observed + Inferred values
% S    - Posterior uncertainty for each code Cx{Nmx(Km(Km+1)/2)}
% W    - Nx1 Observation weights
% E    - 1xP Binning uncertainty
% Z    - NxK Responsibilities
% C    - Nx1 Missing values "code"
% L    - List of codes present (saves a tiny bit of time if provided)
%
% SS0 - 1xK   0th order suff stat (sum of resp)
% SS1 - PxK   1st order suff stat (weighted sum of intensities)
% SS2 - PxPxK 2nd order suff stat (weighted sum of squared intensities)
%
% Compute sufficient statistics up to 2nd order, taking into account
% inferred values and their uncertainty.

if nargin == 0
    help spm_gmm_lib>suffstat
    error('Not enough argument. Type ''help spm_gmm_lib>suffstat'' for help.');
end
if ~ischar(varargin{1})
    [varargout{1:nargout}] = suffstat_default(varargin{:});
    return
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case {'base'}
        [varargout{1:nargout}] = suffstat_base(varargin{:});
    case {'infer'}
        [varargout{1:nargout}] = suffstat_infer(varargin{:});
    case {'bin'}
        [varargout{1:nargout}] = suffstat_bin(varargin{:});              
    otherwise
        help spm_gmm_lib>suffstat
        error('Unknown function %s. Type ''help spm_gmm_lib>suffstat'' for help.', id)
end

% =========================================================================
function [SS0,SS1,SS2] = suffstat_default(X, Z, W)
% FORMAT [SS0,SS1,SS2] = suffstat_default(X, Z, W)
%
% Compute sufficient statistics (up to 2nd order)

%--------------------------------------------------------------------------
% Dimensions
N = size(X, 1);
P = size(X, 2);
K = size(Z, 2);

%--------------------------------------------------------------------------
% Use double precision
X = double(X);
Z = double(Z);
W = double(W);

%--------------------------------------------------------------------------
% Weight responsibilities
Z = bsxfun(@times, Z, W);

%--------------------------------------------------------------------------
% Oth order
SS0 = sum(Z, 1, 'omitnan');

%--------------------------------------------------------------------------
% 1st order
SS1 = sum(bsxfun(@times, X, reshape(Z, [N 1 K])), 1, 'omitnan');
SS1 = reshape(SS1, [P K]);

%--------------------------------------------------------------------------
% 1nd order
SS2 = zeros(P,P,K, 'like', Z);
for i=1:P
    SS2(i,i,:) = reshape(sum(bsxfun(@times, Z, X(:,i).^2),1,'omitnan'), [1 1 K]);
    for j=i+1:P
        SS2(i,j,:) = reshape(sum(bsxfun(@times, Z, X(:,i).*X(:,j)),1,'omitnan'), [1 1 K]);
        SS2(j,i,:) = SS2(i,j,:);
    end
end

% =========================================================================
function [SS0,SS1,SS2] = suffstat_base(X, Z, W, codes)
% FORMAT [{SS0},{SS1},{SS2}] = suffstat_base(X, Z, W, {C,L})
%
% Compute sufficient statistics (up to 2nd order)

if nargin < 3
    W = 1;
end


C  = [];
L  = [];

%--------------------------------------------------------------------------
% Read input arguments
if nargin >= 4
    if ~iscell(codes)
        C = codes;
    else
        if numel(codes) >= 1
            C = codes{1};
            if numel(codes) >= 2
                L = codes{2};
            end
        end
    end
    if isempty(L)
        L = unique(C);
    end
end


%--------------------------------------------------------------------------
% Use double precision
X = double(X);
Z = double(Z);
W = double(W);

%--------------------------------------------------------------------------
% Dimensions
N = size(X,1);
P = size(X,2);
K = size(Z,2);
Z = bsxfun(@times, Z, W); % Multiply resp with observation count
if isempty(L)
    L = 2^P-1;
end

SS0 = cell(1,numel(L));
SS1 = cell(1,numel(L));
SS2 = cell(1,numel(L));

%--------------------------------------------------------------------------
% Sum missing data
for i=1:numel(L)
    c       = L(i);
    missing = ~code2bin(c, P);
    if isempty(C),  msk = ones(N, 1, 'logical');
    else,           msk     = (C == c);
    end
    Pm      = sum(missing);
    Po      = P-Pm;
    Nm      = sum(msk);
    if Nm == 0, continue; end

    X1 = X(msk,~missing);
    Z1 = Z(msk,:);
    
    % ---------------------------------------------------------------------
    % Oth order moment
    SS0{i} = sum(Z1, 1);
    if nargout == 1, return; end

    % ---------------------------------------------------------------------
    % 1st order moment
    SS1{i} = sum(bsxfun(@times, X1, reshape(Z1, [Nm 1 K])), 1);
    SS1{i} = reshape(SS1{i}, [Po K]);
    if nargout == 2, return; end

    % ---------------------------------------------------------------------
    % 2nd order moment
    SS2{i} = zeros(Po,Po,K, 'like', Z);
    for k=1:Po
        SS2{i}(k,k,:) = reshape(sum(bsxfun(@times, Z1, X1(:,k).^2),1), [1 1 K]);
        for j=k+1:Po
            SS2{i}(k,j,:) = reshape(sum(bsxfun(@times, Z1, X1(:,k).*X1(:,j)),1), [1 1 K]);
            SS2{i}(j,k,:) = SS2{i}(k,j,:);
        end
    end
end


% =========================================================================
function [SS0,SS1,SS2] = suffstat_infer(lSS0, lSS1, lSS2, cluster, L)
% FORMAT [SS0,SS1,SS2] = suffstat_infer(SS0, SS1, SS2, {MU,A}, L)
%
% lSS0 - List of sufficient statistics for each pattern of missing data
% lSS1 - List of sufficient statistics for each pattern of missing data
% lSS2 - List of sufficient statistics for each pattern of missing data
% MU   - Clusters' mean
% A    - Clusters' precision matrix
% L    - List of unique codes
% 
% Compute "missing" 1st/2nd order statistics based on a simplified
% inferrence of missing values.
% simplified = inference is performed cluster-wise

MU = [];
A  = [];

%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(cluster)
    MU = cluster;
else
    if numel(cluster) >= 1
        MU = cluster{1};
        if numel(cluster) >= 2
            A = cluster{2};
        end
    end
end
if nargin < 5
    L = [];
end

%--------------------------------------------------------------------------
% Use double precision
MU = double(MU);
A  = double(A);

%--------------------------------------------------------------------------
% Dimensions
P = size(MU,1);
K = size(MU,2);

SS0 = zeros(1,K, 'like', lSS0{1});
if nargout > 1
    SS1 = zeros(P,K, 'like', lSS1{1});
    if nargout > 2
        SS2 = zeros(P,P,K, 'like', lSS2{1});
    end
end

for i=1:numel(L)
    c        = L(i);
    observed = code2bin(c, P);
    missing  = ~observed;
    
    for k=1:K
        % -----------------------------------------------------------------
        % 0th order moment
        SS0k   = lSS0{i}(k);
        SS0(k) = SS0(k) + SS0k;
        
        if nargout > 1
        % -----------------------------------------------------------------
        % 1st order moment
            SS1k = lSS1{i}(:,k);
            MUo  = MU(observed,k);
            MUm  = MU(missing,k);
%             SA   = A(missing,missing,k) \ A(missing,observed,k);
            SA   = spm_matcomp('inv',A(missing,missing,k)) * A(missing,observed,k);
            
            % 1) observed
            SS1(observed,k) = SS1(observed,k) + SS1k;
        
            % 2) missing
            % > t = mu(m) + A(m,m) \ A(m,o) * (mu(o) - g)
            SS1(missing,k) = SS1(missing,k) + SS0k * MUm;
            SS1(missing,k) = SS1(missing,k) + SA * (SS0k * MUo - SS1k);
        end
        
        if nargout > 2
        % -----------------------------------------------------------------
        % 2nd order moment: quadratic terms
            SS2k = lSS2{i}(:,:,k);
        
            % 0) precompute stuff
            MUMUm = SS0k * (MUm * MUm.');
            MUMUo = SS0k * (MUo * MUo.');
            GMUo  = SS1k * MUo.';
            GMUm  = SS1k * MUm.';
        
            % 1) observed x observed
            SS2(observed,observed,k) = SS2(observed,observed,k) + SS2k;
            
            % 2) missing x observed
            tmp = GMUm.' + SA * (GMUo.' - SS2k);
            SS2(missing,observed,k) = SS2(missing,observed,k) + tmp;
            SS2(observed,missing,k) = SS2(observed,missing,k) + tmp.';
            
            % 3) missing x missing
            SS2(missing,missing,k) = SS2(missing,missing,k) + MUMUm;
            tmp = SA * (SS0k * MUo - SS1k) * MUm.';
            SS2(missing,missing,k) = SS2(missing,missing,k) + tmp + tmp';
            SS2(missing,missing,k) = SS2(missing,missing,k) ...
                + SA * (SS2k + MUMUo - GMUo.' - GMUo) * SA.';
    
            % 4) uncertainty ~ missing
            SS2(missing,missing,k) = SS2(missing,missing,k) ...
                + SS0k * spm_matcomp('inv', A(missing,missing,k));
%                 + SS0k * inv(A(missing,missing,k));
        end
    end
end


% =========================================================================
function SS2 = suffstat_bin(E, Z, W, codes)
% FORMAT SS2 = suffstat_bin(E, Z, W, {C,L})
%
% E - Variance in each modality due to binning
% Z - Responisbilities
% W - Observation weights
% C - Missing code image
% L - List of unique codes
% P - Observed space dimension
% 
% Compute "uncertainty" 2nd order statistics based on the posterior
% precision matrix about inferred values.

C  = [];
L  = [];
if nargin < 4
    if nargin < 3
        W = 1;
    end
end

%--------------------------------------------------------------------------
% Read input arguments
if nargin >= 5
    if ~iscell(codes)
        C = codes;
    else
        if numel(codes) >= 1
            C = codes{1};
            if numel(codes) >= 2
                L = codes{2};
            end
        end
    end
    if isempty(L)
        L = unique(C);
    end
end

%--------------------------------------------------------------------------
% Use double precision
E = double(E);
Z = double(Z);
W = double(W);

%--------------------------------------------------------------------------
% Dimensions
N = size(Z,1);
K = size(Z,2);
P = size(E,2);
if sum(E) == 0
    SS2 = zeros(P,P,size(Z,2), 'like', Z);
    return
end
if isempty(L)
    L = 2^P - 1; % None missing
end
Z   = bsxfun(@times, Z, W); % Multiply resp with observation count
SS2 = zeros(P,P,K, 'like', Z);

% -------------------------------------------------------------------------
% 2nd order moment: uncertainty ~ binning
for i=1:numel(L)
    c        = L(i);
    observed = code2bin(c, P);
    Pp       = sum(observed);
    if Pp == 0, continue; end
    
    if isempty(C) && (c == 2^P-1),  msk = ones(N,1,'logical');
    else,                           msk = (C == c);  end
    Nm = sum(msk);
    if Nm == 0, continue; end
    
    
    list_p = 1:P;
    list_p = list_p(observed);
    for p=list_p
        if size(E,1)==1
            SS2(p,p,:) = SS2(p,p,:) ...
                + bsxfun(@times, reshape(sum(Z(msk,:), 1), [1 1 K]), E(p));
        else
            SS2(p,p,:) = SS2(p,p,:) ...
                + reshape(sum(bsxfun(@times, Z(msk,:), E(msk,p)), 1), [1 1 K]);
        end
    end
end

% =========================================================================
function c = const(mean,prec,L)
% FORMAT c = spm_gmm_lib('const', {MU,b}, {V,n})
% FORMAT c = spm_gmm_lib('const', {MU,b}, {A})
% FORMAT c = spm_gmm_lib('const', {MU},   {A})
% FORMAT c = spm_gmm_lib('const', ..., L)
% MU - (Expected) mean
% b  - Mean df (if isempty or 0 -> no Bayesian prior)
% V  - Scale matrix     (if not n isempty or 0) 
% A  - Precision matrix (if n isempty or 0)
% n  - Precision df (if isempty or 0 -> no Bayesian prior)
%
% L - If provided, compute one term for each combination of 
%             missing data
%
% Compute the constant term (w.r.t. voxels) of each Gaussian 
% (expected) log-likelihood.

MU = [];
b  = [];
V  = []; % It can actually be A (when n == 0)
n  = [];

%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(mean)
    MU = mean;
else
    if numel(mean) >= 1
        MU = mean{1};
        if numel(mean) >= 2
            b = mean{2};
        end
    end
end
if ~iscell(prec)
    V = prec;
else
    if numel(prec) >= 1
        V = prec{1};
        if numel(prec) >= 2
            n = prec{2};
        end
    end
end


%--------------------------------------------------------------------------
% Use double precision
MU = double(MU);
b  = double(b);
V  = double(V);
n  = double(n);

%--------------------------------------------------------------------------
% Dimensions
P = size(MU,1);
K = size(MU,2);
    
if nargin == 3 && ~isempty(L)
%--------------------------------------------------------------------------
% Marginal distribution (do not use missing dimensions)

% Assume that:
% . A is a KxK positive-definite matrix
% . A ~ W_K(V,n)
% . S = inv(A)
% . A is partitioned as [A11 A12: A12' A22] with A11 PxP
% . V is partitioned as [V11 V12: V12' V22] with V11 PxP
% . S is partitioned as [S11 S12: S12' S22] with S11 PxP
% Then
% . A11 ~ W_P(V11,n)
% . inv(S11) = A11 - A12*A22\A12'
% . inv(S11) ~ W_P(V11 - V12*V22\V12', n - K + P)
% This allows us to compute E[inv(S11)] and E[ln|S11]], which are needed to
% compute the expected marginal distribution within each cluster.

    c = zeros(numel(L), K, 'like', MU);
    for i=1:numel(L)
        code = L(i);
        observed = code2bin(code, P);
        missing  = ~observed;
        Pm = sum(missing);
        Po = sum(observed);
        for k=1:K
%             Vo = V(observed,observed,k) - V(observed,missing,k)*(V(missing,missing,k)\V(missing,observed,k));
            Vo = V(observed,observed,k) - V(observed,missing,k)*(spm_matcomp('inv',V(missing,missing,k))*V(missing,observed,k));
            c(i,k) = - 0.5 * Po * log(2*pi);
            if sum(n) > 0
                no = n(k) - Pm;
                c(i,k) = c(i,k) + 0.5 * spm_prob('W','ELogDet',Vo,no) ...
                                - 0.5 * no * MU(observed,k).' * Vo * MU(observed,k);
            else
                c(i,k) = c(i,k) + 0.5 * spm_matcomp('LogDet',Vo) ...
                                - 0.5 * MU(observed,k).' * Vo * MU(observed,k);
            end
            if sum(b) > 0
                c(i,k) = c(i,k) - 0.5 * Po / b(k);
            end
        end
        
    end
    
else
%--------------------------------------------------------------------------
% No missing dimensions
    c = zeros(1,K, 'like', MU);
    for k=1:K
        c(k) = - 0.5 * P * log(2*pi);
        if sum(n) > 0
            c(k) = c(k) + 0.5 * spm_prob('W','ELogDet',V(:,:,k),n(k)) ...
                        - 0.5 * n(k) * MU(:,k)' * V(:,:,k) * MU(:,k);
        else
            c(k) = c(k) + 0.5 * spm_matcomp('LogDet',V(:,:,k)) ...
                        - 0.5 * MU(:,k)' * V(:,:,k) * MU(:,k);
        end
        if sum(b) > 0
            c(k) = c(k) - 0.5 * P / b(k);
        end
    end
end

% =========================================================================
function [MU,A,b,V,n] = updateclusters(SS0,SS1,SS2,pr)
% FORMAT [MU,A,b,V,n] = updateclusters(SS0,SS1,SS2,{MU0,b0,V0,n0})
% SS0 - 0th order sufficient statistics (sum Z_i)
% SS1 - 1st order sufficient statistics (sum Z_i * X_i)
% SS2 - 2nd order sufficient statistics (sum Z_i * (X_i * X_i'))
% pr  - List of prior Gauss-Wishart parameters.
%
% Compute posterior GMM parameters from suff stats.

if nargin<4, pr =[]; end

K  = numel(SS0);
MU0 = [];
b0  = [];
V0  = [];
n0  = [];
if numel(pr) >= 1
    MU0 = pr{1};
    if numel(pr) >= 2
        b0 = pr{2};
        if numel(pr) >= 3
            V0 = pr{3};
            if numel(pr) >= 4
                n0 = pr{4};
            end
        end
    end
end

%--------------------------------------------------------------------------
% Use double precision
MU0 = double(MU0);
b0  = double(b0);
V0  = double(V0);
n0  = double(n0);

% -------------------------------------------------------------------------
% Mean
if sum(b0) == 0
    % ---------------------------------------------------------------------
    % Without prior
    b  = [];
    MU = bsxfun(@rdivide, SS1, SS0);
else
    % ---------------------------------------------------------------------
    % With prior
    b  = b0 + SS0;
    MU = bsxfun(@rdivide, SS1 + bsxfun(@times,b0,MU0), b);
end

% -------------------------------------------------------------------------
% Scale/Precision
if sum(n0) == 0
    % ---------------------------------------------------------------------
    % Without prior
    n   = [];
    for k=1:K
        SS2(:,:,k) = SS2(:,:,k) / SS0(k) - (MU(:,k) * MU(:,k).');
    end
else
    % ---------------------------------------------------------------------
    % With prior
    n = n0 + SS0;
    for k=1:K
        SS2(:,:,k) = SS2(:,:,k) +   b0(k) * MU0(:,k) * MU0(:,k).' ...
                                -    b(k) * MU(:,k)  * MU(:,k).' ...
                                + spm_matcomp('Inv',V0(:,:,k));
    end
end
V = SS2;
for k=1:K
    V(:,:,k) = spm_matcomp('Inv', V(:,:,k));
end
if sum(n) > 0
    A = bsxfun(@times, V, reshape(n, [1 1 K]));
else
    A = V;
    V = [];
end

% =========================================================================
function [PI,logPI,a] = updateproportions(SS0, a0)
% FORMAT [PI,logPI,a,ll] = updateproportions(SS0, a0)
%
% SS0 - 1xK 0th order sufficient statistics (sum of responsibilities)
% a0  - 1xK Dirichlet prior (can be 0)
%
% PI    - 1xK Cluster proportion posterior expected value
% logPI - 1xK ln(PI) or E[ln(PI)] (if Bayesian)
% a     - Dirichlet posterior (if Bayesian)
% ll    - 1x1 Lower bound: E[ln p(PI|a)] - E[ln q(PI)]
%
% Bayesian or ML update of cluster proportions.

a = a0 + SS0;
if sum(a0(:))
% Bayesian
    % expected values
    logPI = psi(a) - psi(sum(a));
    PI    = a ./ sum(a(:));
else
% Maximum Likelihood
    a     = max(a, eps);
    PI    = a ./ sum(a(:));
    logPI = log(PI);
end
% =========================================================================

% =========================================================================
function [GaussPrior,extras] = updatehyperpars(cluster,GaussPrior,varargin)
% FORMAT [GaussPrior,extras] = updatehyperpars(cluster,GaussPrior,varargin)
%
% REQUIRED
% cluster    - 1xS cell array where cluster{s} = {{MU,b},{V,n}}
% GaussPrior - {MU0,b0,V0,n0}
%
% OPTIONAL
% constrained - Optimise hierarchical prior on V [false]
% figname     - Postfix added to figure name ['']
% verbose     - Verbosity level: [false]=quiet, true=display
%
% OUTPUT
% GaussPrior - New {MU0,b0,V0,n0}
% extras     - Struct with lower bound information, etc.
%
% Update of VB-GMM hyper-parameters (m,b,W,n).

% Parse optional arguments
%--------------------------------------------------------------------------
p = inputParser;
p.FunctionName = 'updatehyperpars';
p.addParameter('constrained',0,@islogical);
p.addParameter('figname','',@ischar);
p.addParameter('verbose',0,@islogical);
p.addParameter('lkp',[],@isnumeric);
p.parse(varargin{:});
constrained = p.Results.constrained;
figname     = p.Results.figname;
verbose     = p.Results.verbose;
lkp         = p.Results.lkp;

% Parameters
S = numel(cluster); % Number of posteriors

m0 = GaussPrior{1};
b0 = GaussPrior{2};
W0 = GaussPrior{3};
n0 = GaussPrior{4};

N = size(m0,1);
K = size(m0,2);

% pre-allocate
LogDetW0  = zeros(size(n0));
V         = zeros(size(W0));
p         = zeros(size(n0));
p0        = 0;

% -------------------------------------------------------------------------
%   Gauss-Wishart "mean" parameters
% -------------------------------------------------------------------------

for k=1:K
    
    % ---------------------------------------------------------------------
    % Update m0 (mode, closed-form)
    
    Lambda   = 0;
    LambdaMu = 0;
    for s=1:S
        [m,~,W,n] = get_posteriors(cluster,s);
        Lambda    = Lambda   + n(k)*W(:,:,k);
        LambdaMu  = LambdaMu + n(k)*W(:,:,k)*m(:,k);
    end
    m0(:,k) = Lambda \ LambdaMu;
    
    % ---------------------------------------------------------------------

    
    % ---------------------------------------------------------------------
    % Update b0 (mode, closed-form)

    b0(k)= 0;
    for s=1:S
        [m,b,W,n] = get_posteriors(cluster,s);
        m1 = m(:,k) - m0(:,k);
        b0(k) = b0(k) + m1.' * (n(k)*W(:,:,k)) * m1 + N/b(k);
    end
    b0(k) = N*S/b0(k);
    
    % ---------------------------------------------------------------------

end


% =========================================================================
% NOT CONSTRAINED
if ~constrained
    
    % ---------------------------------------------------------------------
    %   Gauss-Wishart "precision" parameters
    % ---------------------------------------------------------------------

    for k=1:K
        
        % ---
        % Set up some constants
        sumLogDet = 0;
        sumPsi    = 0;
        Wn        = 0;
        for s=1:S
            [~,~,W,n] = get_posteriors(cluster,s);
            sumLogDet = sumLogDet + spm_matcomp('LogDet', W(:,:,k));
            sumPsi    = sumPsi    + spm_prob('DiGamma', n(k)/2, N);
            Wn        = Wn        + n(k)*W(:,:,k);
        end
        sumLogDet = sumLogDet/S;
        sumPsi    = sumPsi/S;
        Wn        = Wn/S;
            
        % -----------------------------------------------------------------
        % Update n0 (mode, Gauss-Newton [convex])
        E = inf;
        for gniter=1:1000
            
            % -------------------------------------------------------------
            % Update W0 (mode, closed-form)
            W0(:,:,k)   = Wn/n0(k);
            LogDetW0(k) = spm_matcomp('LogDet', W0(:,:,k));
            % -------------------------------------------------------------
            
            % ---
            % Objective function
            Eprev = E;
            E = 0.5*S*n0(k)*( LogDetW0(k) - sumLogDet - sumPsi ) ...
                + S*spm_prob('LogGamma', n0(k)/2, N);
            
            if E == Eprev
                break;
            end

            % ---
            % Gradient & Hessian
            g = 0.5*S*( LogDetW0(k) - sumLogDet - sumPsi ...
                         + spm_prob('DiGamma', n0(k)/2, N) );
            H = S/4*spm_prob('DiGamma', n0(k)/2, N, 1);

            % ---
            % Update
            n0(k) = max(n0(k) - H\g, N-1+2*eps);
            
        end
        % -----------------------------------------------------------------
    
    end
    
    % ---------------------------------------------------------------------
    %   Save results
    % ---------------------------------------------------------------------
    extras.b   = b0;
    extras.m   = m0;
    extras.n   = n0;
    extras.W   = W0;
    extras.ldW = LogDetW0;
    extras.lb  = 0;
    
    GaussPrior{1} = m0;
    GaussPrior{2} = b0;
    GaussPrior{3} = W0;
    GaussPrior{4} = n0;
    
% =========================================================================
% CONSTRAINED
else

    lb = -inf;
    for em=1:50
        % ---
        % Starting estimate
        if p0 == 0
            p0 = 0;
            V0 = 0;
            for k=1:K
                p0 = p0 + S*n0(k);
                for s=1:S
                    [~,~,W,n] = get_posteriors(cluster,s);
                    V0 = V0 + spm_matcomp('Inv', n(k)*W(:,:,k));
                end
            end
            p0 = p0/K;
            V0 = V0/K;
        end

        % -----------------------------------------------------------------
        %   Gauss-Wishart "precision" parameters
        % -----------------------------------------------------------------

        for k=1:K

            % ---
            % Set up some constants
            % > compute sum E[logdet W] and sum psi(nu/2)
            logDetW  = 0;
            psiN     = 0;
            Lambda   = 0;
            for s=1:S
                [~,~,W,n] = get_posteriors(cluster,s);
                logDetW = logDetW  + spm_matcomp('Logdet', W(:,:,k));
                psiN    = psiN     + spm_prob('DiGamma', n(k)/2, N);
                Lambda  = Lambda   + n(k)*W(:,:,k);
            end
            logDetW  = logDetW/S;
            psiN = psiN/S;


            % -------------------------------------------------------------
            % Update n0 (mode, Gauss-Newton [convex])
            E = inf;
            for gniter=1:100

                % ---------------------------------------------------------
                % Update {p,V} for W0 (posterior, closed form)
                p(k)       = p0 + S*n0(k);
                V(:,:,k)   = spm_matcomp('Inv', spm_matcomp('Inv', V0) + Lambda);
                % Useful values
                W0(:,:,k)   = spm_matcomp('Inv', spm_prob('W', 'E', V(:,:,k), p(k)));
                LogDetW0(k) = -spm_prob('W', 'Elogdet', V(:,:,k), p(k));
                % ---------------------------------------------------------

                % ---
                % Objective function                
                E1 = S*n0(k)/2 * (LogDetW0(k) - logDetW - psiN) ...
                     + S*spm_prob('LogGamma', n0(k)/2, N);
                E = [E E1];
                
                subgain = spm_misc('get_gain',E);                
                if subgain < 1e-6
                    % Finished
                    break
                end

                % ---
                % Gradient & Hessian
                g = S/2*(LogDetW0(k) - logDetW - psiN + spm_prob('DiGamma', n0(k)/2, N));
                H = S/4 * spm_prob('DiGamma', n0(k)/2, N, 1);

                % ---
                % Update
                n0(k) = max(n0(k) - H\g, N-1+2*eps);
            end
            % ------------------------------------------------------------

        end


        % -----------------------------------------------------------------
        %   Inverse-Wishart parameters
        % -----------------------------------------------------------------

        % ---
        % Set up some constants
        % > compute sum Logdet(psi) and sum psi(m/2)
        sumlogV = 0;
        sumPsi  = 0;
        pV      = 0;
        for k=1:K
            sumlogV = sumlogV + spm_matcomp('LogDet', V(:,:,k));
            sumPsi  = sumPsi  + spm_prob('DiGamma', p(k)/2, N);
            pV      = pV      + p(k)*V(:,:,k);
        end
        sumlogV = sumlogV/K;
        sumPsi  = sumPsi/K;
        pV      = pV/K;


        % -----------------------------------------------------------------
        % Update p0 (mode, Gauss-Newton [convex])
        E = inf;
        for gniter=1:1000

            % -------------------------------------------------------------
            % Update V0 (closed-form)
            V0 = pV/p0;
            LogDetV0 = spm_matcomp('LogDet', V0);
            % -------------------------------------------------------------

            % ---
            % Objective function
            Eprev = E;
            E = p0*K/2*( N*LogDetV0 - sumlogV - sumPsi ) + K*spm_prob('LogGamma', p0/2, N);
            if E == Eprev
                break;
            end

            % ---
            % Gradient & Hessian
            g = K/2*( LogDetV0 - sumlogV - sumPsi + spm_prob('DiGamma', p0/2, N) );
            H = K/4*spm_prob('DiGamma', p0/2, N, 1);

            % ---
            % Update
            p0 = max(p0 - H\g, N-1+2*eps);

        end
        % -----------------------------------------------------------------
    
        
        % ---
        % Objective function        
        nlb  = 0;
        for k=1:K
            nlb  = nlb - spm_prob('Wishart', 'kl', V(:,:,k), p(k), V0, p0);
        end
        
        lb   = [lb nlb];      
        gain = spm_misc('get_gain',lb);        
        
        if gain < 1e-3
            % Finished
            break
        end
        
    end % < "EM" loop

    % ---------------------------------------------------------------------
    %   Save results
    % ---------------------------------------------------------------------
    extras.b   = b0;
    extras.m   = m0;
    extras.n   = n0;
    extras.W   = W0;
    extras.ldW = LogDetW0;
    extras.V   = V;
    extras.p   = p;
    extras.V0  = V0;
    extras.p0  = p0;
    extras.lb  = 0;
    for k=1:K
        extras.lb  = extras.lb - spm_prob('Wishart', 'kl', V(:,:,k), p(k), V0, p0);
    end  
    
    GaussPrior{1} = m0;
    GaussPrior{2} = b0;
    GaussPrior{3} = W0;
    GaussPrior{4} = n0;
end

if verbose
    % Visualise results
    spm_gmm_lib('plot','gaussprior',GaussPrior,lkp,figname);
end
% =========================================================================

%--------------------------------------------------------------------------
% Lower bound 
% -------------------------------------------------------------------------


% =========================================================================
function [lb,const] = marginalsum(SS0, SS1, SS2, mean, prec, L, SS2b)
% [lb,const] = spm_gmm_lib('marginalsum', SS0, SS1, SS2, MU, A, L, SS2b)
% [lb,const] = spm_gmm_lib('marginalsum', SS0, SS1, SS2, {MU,b}, {V,n}, L, SS2b)
% 
% SS0       - {1xK}   Zero-th order moment (per config)
% SS1       - {PxK}   First   order moment (per config)
% SS1       - {PxPxK} Second  order moment (per config)
% MU        - PxK     Means
% A/V       - PxPxK   Precision/Scale matrices
% b         - 1xK     Mean degrees of freedom
% n         - 1xK     Precision degrees of freedom
% L         - Mx1     List of existing codes
% SS2b      - PxPxK   Binning uncertainty
%
% lb        -         Sum of (expected) marginal likelihoods
% const     - MxK     Constant terms
%
% Compute the expected log-likelihood of each observation belonging to each
% cluster: lb = sum_{i,k} E[z_ik] E[ln p(g(i) | MU_k,A_k)]


MU = [];
A  = [];
V  = [];
n  = [];

%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(mean)
    MU = mean;
else
    if numel(mean) >= 1
        MU = mean{1};
        if numel(mean) >= 2
            b = mean{2};
        end
    end
end
if ~iscell(prec)
    A = prec;
else
    if numel(prec) >= 1
        A = prec{1};
        if numel(prec) >= 2
            V = A;
            n = prec{2};
        end
    end
end
if nargin <= 6
    L = [];
end
if nargin <= 7
    SS2b = 0;
end

% -------------------------------------------------------------------------
% Dimensions
P  = size(MU,1);
K  = size(MU,2);
if isempty(L), L = 2^P - 1;    end % None missing

% -------------------------------------------------------------------------
% Constant term
const = spm_gmm_lib('const', mean, prec, L);

lb  = 0;

% -------------------------------------------------------------------------
% For each combination of missing voxels
for i=1:numel(L)
    
    % ---------------------------------------------------------------------
    % Get mask of missing values and modalities (with this particular code)
    c        = L(i);
    observed = code2bin(c, P);
    Po       = sum(observed);
    if Po == 0, continue; end
    missing  = ~observed;
    Pm       = P-Po;
    
    % ---------------------------------------------------------------------
    % Initialise with constant term
    lb = lb + sum(const(i,:) .* SS0{i});
    
    % ---------------------------------------------------------------------
    % Non constant terms
    for k=1:K
        
        % /!\ Sub-covariance is different from the inverse sub-precision
        % ML case:
        %   inv(S(o,o)) = A(o,o) - A(o,m)*A(m,m)\A(m,o)
        % Bayesian case:
        %   inv(S(o,o)) ~ W(V(o,o) - V(o,m)*V(m,m)\V(m,o), n - Pm)
        if sum(n) > 0
%             Ao = V(observed,observed,k) - V(observed,missing,k)*(V(missing,missing,k)\V(missing,observed,k));
            Ao = V(observed,observed,k) - V(observed,missing,k)*(spm_matcomp('inv',V(missing,missing,k))*V(missing,observed,k));
            Ao = (n(k)-Pm) * Ao;
        else
%             Ao = A(observed,observed,k) - A(observed,missing,k)*(A(missing,missing,k)\A(missing,observed,k));
            Ao = A(observed,observed,k) - A(observed,missing,k)*(spm_matcomp('inv',A(missing,missing,k))*A(missing,observed,k));
        end
        
        % 1) obs x mean
        lb = lb + SS1{i}(:,k).' * Ao * MU(observed,k);
        
        % 1) obs x obs
        lb = lb - 0.5 * trace(Ao * SS2{i}(:,:,k));
        
        % 3) Binning uncertainty
        lb = lb - 0.5 * trace(Ao * SS2b);
        
    end
        
end

% =========================================================================
function varargout = kl(varargin)
% Useful KL-divergences for Gaussian Mixture modelling
%
% [klMU,klA] = spm_gmm_lib('kl', 'GaussWishart', {MU,b}, {V,n}, {MU0,b0}, {V0,n0})
% > KL divergence between two Gauss-Wishart distributions
%
% klP = spm_gmm_lib('kl', 'Dirichlet', a, a0)
% > KL divergence between two Dirichlet distributions
%
% klZ = spm_gmm_lib('kl', 'Categorical', Z, W, logPI)
% > KL divergence between two Categorical distributions

if nargin == 0
    help spm_gmm_lib>kl
    error('Not enough argument. Type ''help spm_gmm_lib>kl'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case {'gausswishart','gw'}
        [varargout{1:nargout}] = kl_gausswishart(varargin{:});
    case {'dirichlet','d'}
        [varargout{1:nargout}] = kl_dirichlet(varargin{:});
    case {'categorical','cat','c'}
        [varargout{1:nargout}] = kl_categorical(varargin{:});              
    otherwise
        help spm_gmm_lib>kl
        error('Unknown function %s. Type ''help spm_gmm_lib>kl'' for help.', id)
end

% =========================================================================
function klZ = kl_categorical(Z, W, logPI)

% Initialise
klZ = zeros(1, 'like', Z);

% E[ln p(Z|PI)] (prior ~ responsibilities)
klZ = klZ + sum(sum(bsxfun(@times,Z,logPI), 2) .* W, 1);

% -E[ln q(Z)] (posterior ~ responsibilities))
klZ = klZ - sum(sum(Z .* log(max(Z,eps)), 2) .* W, 1);

% =========================================================================
function klP = kl_dirichlet(a, a0)

klP = zeros(1, 'like', a);
K   = numel(a);
if sum(a0) > 0
    % prior
    klP = gammaln(sum(a0)) - sum(gammaln(a0));
    klP = klP + sum((a0-1) .* (psi(a) - K*psi(sum(a))));
    % posterior
    klP = klP - gammaln(sum(a)) - sum(gammaln(a));
    klP = klP - sum((a-1) .* (psi(a) - K*psi(sum(a))));
end

% =========================================================================
function [klMU,klA] = kl_gausswishart(mean,prec,mean0,prec0)

MU  = [];
b   = [];
V   = []; % It can actually be A (when n == 0)
n   = [];
MU0 = [];
b0  = [];
V0  = [];
n0  = [];

%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(mean)
    MU = mean;
else
    if numel(mean) >= 1
        MU = mean{1};
        if numel(mean) >= 2
            b = mean{2};
        end
    end
end
if ~iscell(prec)
    V = prec;
else
    if numel(prec) >= 1
        V = prec{1};
        if numel(prec) >= 2
            n = prec{2};
        end
    end
end
if nargin >= 3
    if ~iscell(mean0)
        MU0 = mean0;
    else
        if numel(mean0) >= 1
            MU0 = mean0{1};
            if numel(mean0) >= 2
                b0 = mean0{2};
            end
        end
    end
end
if nargin >=4
    if ~iscell(prec0)
        V0 = prec0;
    else
        if numel(prec0) >= 1
            V0 = prec0{1};
            if numel(prec0) >= 2
                n0 = prec0{2};
            end
        end
    end
end

%--------------------------------------------------------------------------
% Read input arguments
P = size(MU,1);
K = size(MU,2);
LogDetA = zeros(1,K, 'like', V);
if sum(n) > 0
    A = bsxfun(@times, V, reshape(n, [1 1 K]));
    for k=1:K
        LogDetA(k) = spm_prob('W','ELogDet',V(:,:,k),n(k));
    end
else
    A = V;
    for k=1:K
        LogDetA(k) = spm_matcomp('LogDet',A(:,:,k));
    end
end

% Lower bound
klMU = zeros(1, 'like', MU);
klA  = zeros(1, 'like', A);
for k=1:K
    % + prior
    if sum(b0) > 0
        % prior
        klMU = klMU - P*log(2*pi) ...
                    + P*log(b0(k)) ...
                    + LogDetA(k) ...
                    - b0(k)*(MU(:,k)-MU0(:,k)).'*A(:,:,k)*(MU(:,k)-MU0(:,k)) ...
                    - P*b0(k)/b(k);
        % posterior
        klMU = klMU + P*log(2*pi) ...
                    - P*log(b(k)) ...
                    - LogDetA(k) ...
                    + P;
    end
    if sum(n0) > 0
        klA = klA - spm_prob('W', 'kl', V(:,:,k), n(k), V0(:,:,k), n0(k));
    end
end
klMU = 0.5 * klMU;

%--------------------------------------------------------------------------
% "Missing code" image 
% -------------------------------------------------------------------------

% =========================================================================
function code = obs2code(X)
% FORMAT code = spm_gmm_lib('obs2code', X)
%
% Compute a "missing code" image for the input observation matrix.

code = double2int(sum(bsxfun(@times, ~isnan(X), 2.^(0:size(X,2)-1)), 2));

% =========================================================================
function bin = code2bin(code, length)
% FORMAT bin = spm_gmm_lib('code2bin', code, length)

bin = dec2bin(code,length) == '1';
bin = bin(end:-1:1);

% =========================================================================
function L = double2int(L)
% FORMAT L = spm_gmm_lib('double2int', L)
%
% Find the best suited integer type to convert L, based on min and max
% values

minval = min(L(:));
maxval = max(L(:));
type   = range2int(maxval,minval);
func   = str2func(type);
L      = func(L);

% =========================================================================
function type = range2int(maxval,minval)
% FORMAT type = spm_gmm_lib('range2int', maxval, minval)
%
% Find the best suited integer type to store integer values in the range
% [minval,maxval]

if nargin < 2
    minval = 0;
end

type     = 'int';
unsigned = minval >= 0;
if unsigned
    type   = ['u' type];
    minval = 0;
else
    minval = numel(dec2base(-minval,2));
end
maxval = numel(dec2base(maxval,2));
nbits  = max(minval,maxval);
if unsigned
    nbits = nbits + 1;
end
if nbits <= 8
    type = [type '8'];
elseif nbits <= 16
    type = [type '16'];
elseif nbits <= 32
    type = [type '32'];
elseif nbits <= 64
    type = [type '54'];
else
    type = 'double';
end

% =========================================================================
function varargout = gmmplot(varargin)
% Custom visualisation tools for Gaussian Mixture modelling
%
% spm_gmm_lib('plot', 'lb', lb, (wintitle))
% > Plot lower bound
%
% spm_gmm_lib('plot', 'gmm', {X,W}, {MU,A}, PI, (wintitle))
% > Plot mixture fit
%
% spm_gmm_lib('plot', 'cat', dm, Z, Template, (wintitle))
% > Plot (categorical) responsibilities and template (if available)
%
% spm_gmm_lib('plot', 'gaussprior', GaussPrior, (wintitle))
% > Plot VB-GMM hyper-parameters
%
% c = spm_gmm_lib('plot', 'cat2rgb', f, pal)
% > Generate an RGB volume from a categorical (e.g. responsibilities) volume.
%
% c = spm_gmm_lib('plot', 'show_cat_img', img, title_nam)
% > Show categorical images
%

if nargin == 0
    help spm_gmm_lib>plot
    error('Not enough argument. Type ''help spm_gmm_lib>plot'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case {'lowerbound','lb'}
        [varargout{1:nargout}] = plot_lowerbound(varargin{:});
    case {'gmm'}
        [varargout{1:nargout}] = plot_gmm(varargin{:});     
    case {'cat'}
        [varargout{1:nargout}] = plot_cat(varargin{:});           
    case {'gaussprior'}
        [varargout{1:nargout}] = plot_GaussPrior(varargin{:});              
    case {'cat2rgb'}
        [varargout{1:nargout}] = cat2rgb(varargin{:});         
    case {'showcatimg'}
        [varargout{1:nargout}] = show_cat_img(varargin{:});              
    otherwise
        help spm_gmm_lib>plot
        error('Unknown function %s. Type ''help spm_gmm_lib>plot'' for help.', id)
end

% =========================================================================
function plot_lowerbound(lb, figname)

% -------------------------------------------------------------------------
% Get figure (create if it does not exist)
if nargin < 2
    figname = '(SPM) Plot GMM Lower Bound';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

% -------------------------------------------------------------------------
% Choose type
if isfield(lb, 'B')
    nrow = 2;
    ncol = 4;
else
    nrow = 2;
    ncol = 3;
end

% -------------------------------------------------------------------------
% Plots
subplot(nrow, ncol, sub2ind([ncol nrow], 1, 1));
plot(lb.sum)
title('Lower Bound')
subplot(nrow, ncol, sub2ind([ncol nrow], 2, 1));
if isfield(lb, 'B')
    plot(sum(lb.X,1) + sum(lb.XB,1));
else
    plot(sum(lb.X,1))
end
box on
title('Observations (E)')
subplot(nrow, ncol, sub2ind([ncol nrow], 3, 1));
plot(sum(lb.Z,1))
box on
title('Responsibilities (KL)')
subplot(nrow, ncol, sub2ind([ncol nrow], 1, 2));
plot(sum(lb.P,1))
box on
title('Proportions (KL)')
subplot(nrow, ncol, sub2ind([ncol nrow], 2, 2));
plot(sum(lb.MU,1))
box on
title('Means (KL)')
subplot(nrow, ncol, sub2ind([ncol nrow], 3, 2));
plot(sum(lb.A,1))
box on
title('Precisions (KL)')
if isfield(lb, 'B')
    subplot(nrow, ncol, sub2ind([ncol nrow], 4, 2));
    plot(sum(lb.B,1))
    box on
    title('Bias Prior')
end
drawnow

% =========================================================================
function plot_cat(Z,Template,ticklabels,figname)
if nargin < 3, ticklabels = {}; end

% -------------------------------------------------------------------------
% Get figure (create if it does not exist)
if nargin<4
    figname = '(SPM) Plot GMM Categorical';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

show_cat_img({Z,Template},{'Z','Template'},ticklabels);
% =========================================================================

% =========================================================================
function plot_gmm(obs, cluster, PI, part, figname)

X  = [];
W  = [];
MU = [];
A  = [];

if ~iscell(obs)
    X = obs;
else
    if numel(obs) >= 1
        X = obs{1};
        if numel(obs) >= 2
            W = obs{2};
        end
    end
end
if ~iscell(cluster)
    MU = cluster;
else
    if numel(cluster) >= 1
        MU = cluster{1};
        if numel(cluster) >= 2
            A = cluster{2};
        end
    end
end

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
if nargin < 5
    figname = '(SPM) Plot GMM';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

% ---------------------------------------------------------------------
% Sizes / colors
P = size(X, 2);
K = size(MU, 2);

if nargin < 4
    lkp = 1:K;
    mg  = ones(1,K);
else
    lkp = part.lkp;
    mg  = part.mg;  
end

PI = mg.*PI(:,lkp);

% Set colors
colors = hsv(max(lkp));

% ---------------------------------------------------------------------
% For each input dimension
for p=1:P
    % -----------------------------------------------------------------
    % Plot histogram and marginal density
    subplot(2, P, p)
    hold on
    % ---------
    % Histogram
    if any(W > 1)
        % Input is already an histogram
        X1 = X(isfinite(X(:,p)),p);
        weights = W(isfinite(X(:,p)));
        [idx,edges] = discretize(X1, 64);
        centres = (edges(2:end) + edges(1:end-1))/2;
        idx2 = zeros(numel(idx), max(idx));
        idx2(sub2ind(size(idx2), 1:numel(idx), idx')) = 1;
        weights = sum(bsxfun(@times, weights, idx2), 1);
        clear idx1 idx2
        minx  = min(X1);
        maxx  = max(X1);
        weights = weights ./ (sum(weights)*(maxx-minx)/numel(weights));
        bar(centres, weights, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
        ymax = max(weights);
    else
        % Input is a list of observations
        [H, edges] = histcounts(X(:,p), 64, 'Normalization', 'pdf');
        centres = (edges(1:end-1) + edges(2:end))/2;
        bar(centres, H, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
        ymax = max(H);
    end
    xlims = [inf -inf];
    % -----------
    % GMM Density
    for k=1:K
        x = linspace(MU(p,k)-3*A(p,p,k)^(-0.5),MU(p,k)+3*A(p,p,k)^(-0.5),100);
        y = PI(k)*normpdf(x, MU(p,k), A(p,p,k)^(-0.5));
        plot(x, y, 'Color', colors(lkp(k),:), 'LineWidth', 1)
        xlims = [min([xlims(1) x]) max([xlims(2) x])];
    end
    xlabel(sprintf('x%d',p))
    ylabel('density')
    xlim(xlims);
    ylim([0 1.1*ymax]);
    box on
    hold off

    % -----------------------------------------------------------------
    % Plot joint density (X1 vs Xj)
    if p > 1
        subplot(2, P, P+p)
        hold on
        for k=1:K
            Mu1     = MU([1 p],k);
            Sigma2  = spm_matcomp('Inv', A([1 p],[1 p],k));
            Sigma   = sqrt(Sigma2);
            [x1,x2] = meshgrid(linspace(Mu1(1)-3*Sigma(1,1),Mu1(1)+3*Sigma(1,1),100)', ...
                               linspace(Mu1(2)-3*Sigma(2,2),Mu1(2)+3*Sigma(2,2),100)');
            y = mvnpdf([x1(:) x2(:)],Mu1',Sigma2);
            contour(x2, x1, reshape(y, [100 100])', 1, 'color', colors(lkp(k),:), 'LineWidth', 1);
        end
        xlabel(sprintf('x%d',p))
        ylabel('x1')
        xlim(xlims); % use same scale as histogram plot for comparison
        box on
        hold off
    end
end

% ---------------------------------------------------------------------
% Plot proportions in the remaining square
subplot(2, P, P+1)
hold on
for k=1:K
    bar(k, PI(k), 'FaceColor', colors(lkp(k),:));
end
xlabel('class')
ylabel('proportion')
box on
hold off
drawnow
% =========================================================================

% =========================================================================
function plot_GaussPrior(GaussPrior,lkp,figname)

if nargin < 2, lkp = []; end

figname0 = '(SPM) GaussPrior';
if nargin==3
    figname0 = [figname0 ' ' figname];
end

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname0);
if isempty(f)
    f = figure('Name', figname0, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

m0 = GaussPrior{1};
W0 = GaussPrior{3};
n0 = GaussPrior{4};
    
P      = size(m0,1);
K      = size(m0,2);

if isempty(lkp)
    lkp = 1:K;
end

colors = hsv(max(lkp));

MU = m0;
A  = bsxfun(@times,W0,reshape(n0,[1 1 K]));

% ---------------------------------------------------------------------
% For each input dimension
for p=1:P
    % -----------------------------------------------------------------
    % Plot histogram and marginal density
    if P>1
        subplot(2, P, p)
    end
    hold on

    xlims = [inf -inf];
    % -----------
    % GMM Density
    for k=1:K
        x = linspace(MU(p,k)-3*A(p,p,k)^(-0.5),MU(p,k)+3*A(p,p,k)^(-0.5),100);
        y = 1/K*normpdf(x, MU(p,k), A(p,p,k)^(-0.5));
        plot(x, y, 'Color', colors(lkp(k),:), 'LineWidth', 1)
        xlims = [min([xlims(1) x]) max([xlims(2) x])];
    end
    xlabel(sprintf('x%d',p))
    ylabel('density')
    xlim(xlims);
    box on
    hold off

    % -----------------------------------------------------------------
    % Plot joint density (X1 vs Xj)
    if p > 1
        subplot(2, P, P+p)
        hold on
        for k=1:K
            Mu1     = MU([1 p],k);
            Sigma2  = spm_matcomp('Inv', A([1 p],[1 p],k));
            Sigma   = sqrt(Sigma2);
            [x1,x2] = meshgrid(linspace(Mu1(1)-3*Sigma(1,1),Mu1(1)+3*Sigma(1,1),100)', ...
                               linspace(Mu1(2)-3*Sigma(2,2),Mu1(2)+3*Sigma(2,2),100)');
            y = mvnpdf([x1(:) x2(:)],Mu1',Sigma2);
            contour(x2, x1, reshape(y, [100 100])', 1, 'color', colors(lkp(k),:), 'LineWidth', 1);
        end
        xlabel(sprintf('x%d',p))
        ylabel('x1')
        xlim(xlims); % use same scale as histogram plot for comparison
        box on
        hold off
    end
end

drawnow
%==========================================================================

%==========================================================================
function c = cat2rgb(f, pal)
% FORMAT c = cat2rgb(f, pal)
% f   - categorical (4D) image.
% pal - palette (Mx3 array or handle to palette function) [hsv]
%
% Generate an RGB volume from a categorical (e.g. responsibilities) volume.

if nargin < 2
    pal = @hsv;
end

if size(f,3)>1
    z = floor(size(f,3)/2) + 1;
    f = f(:,:,z,:);
end

tri = false;
if numel(size(f)) == 4 && size(f, 3) == 1
    tri = true;
    dm  = [size(f) 1 1];
    f   = reshape(f, [dm(1:2) dm(4)]);
end
if isa(pal, 'function_handle')
    pal = pal(size(f,3));
end

dm = [size(f) 1 1];
c   = zeros([dm(1:2) 3]); % output RGB image
s   = zeros(dm(1:2));     % normalising term

for k=1:dm(3)
    s = s + f(:,:,k);
    color = reshape(pal(k,:), [1 1 3]);
    c = c + bsxfun(@times, f(:,:,k), color);
end
if dm(3) == 1
    c = c / max(1, max(s(:)));
else
    c = bsxfun(@rdivide, c, s);
end

if tri
    c = reshape(c, [size(c, 1) size(c, 2) 1 size(c, 3)]);
end
%==========================================================================

%==========================================================================
function show_cat_img(img,title_nam,ticklabels)
if isnumeric(img)
    img = {img};
end 
N       = numel(img);
if nargin < 2, title_nam = cell(1,N); end

dm0    = size(img{1});
K      = dm0(4);
colors = hsv(K);
            
if nargin < 3 || isempty(ticklabels), ticklabels = 1:K; end

if numel(ticklabels) ~= K
    ticklabels = 1:K
end

if dm0(3)==1
    % 2d 
    %----------------------------------------------------------------------
    
    for n=1:N            
        
        subplot(1,N,n)
        
        slice = img{n}(:,:,1,:);
        slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
        slice = squeeze(slice(:,:,:,:));
        slice = permute(slice,[2 1 3]);
        imagesc(slice); axis off xy;      
        
        if ~isempty(title_nam{n})
            title(title_nam{n})
        end
    end
    
    colormap(colors);
    cb = colorbar;
    set(gca, 'clim', [0.5 K+0.5]);
    set(cb, 'ticks', 1:K, 'ticklabels', ticklabels); 
else
    % 3d
    %----------------------------------------------------------------------
    
    for n=1:N  
        subplot(N,3,3*(n - 1) + 1)
        
        slice = img{n}(:,:,floor(dm0(3)/2) + 1,:);
        slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
        slice = squeeze(slice(:,:,:,:));
        slice = permute(slice,[2 1 3]);
        imagesc(slice); axis off xy;  

        if ~isempty(title_nam{n})
            title(title_nam{n})
        end
        
        subplot(N,3,3*(n - 1) + 2)
        
        slice = permute(img{n}(:,floor(dm0(2)/2) + 1,:,:),[3 1 2 4]);
        slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
        slice = squeeze(slice(:,:,:,:));    
        imagesc(slice); axis off xy;  

        if ~isempty(title_nam{n})
            title(title_nam{n})
        end
        
        subplot(N,3,3*(n - 1) + 3)
        
        slice = permute(img{n}(floor(dm0(1)/2) + 1,:,:,:),[2 3 1 4]);
        slice = spm_gmm_lib('plot', 'cat2rgb', slice, colors);
        slice = squeeze(slice(:,:,:,:));
        slice = permute(slice,[2 1 3]);
        imagesc(slice); axis off xy;  
        
        if ~isempty(title_nam{n})
            title(title_nam{n})
        end
        
        colormap(colors);
        cb = colorbar;
        set(gca, 'clim', [0.5 K+0.5]);
        set(cb, 'ticks', 1:K, 'ticklabels', ticklabels);             
    end    
end

drawnow
%==========================================================================

% =========================================================================
function varargout = gmm_extras(varargin)
% Additional GMM functions
%
% gmm = spm_gmm_lib('extras', 'more_gmms', gmm, part)
% > A crude heuristic to replace a single Gaussian by a bunch of Gaussians.
% gmm = spm_gmm_lib('extras', 'collapse_gmms', gmm, part)
% > A crude heuristic to replace a bunch of Gaussian with a single Gaussian.
%

if nargin == 0
    help spm_gmm_lib>extras
    error('Not enough argument. Type ''help spm_gmm_lib>plot'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case {'more_gmms'}
        [varargout{1:nargout}] = more_gmms(varargin{:});      
    case {'collapse_gmms'}
        [varargout{1:nargout}] = collapse_gmms(varargin{:});         
    otherwise
        help spm_gmm_lib>extras
        error('Unknown function %s. Type ''help spm_gmm_lib>extras'' for help.', id)
end
% =========================================================================

% =========================================================================
function [gmm,mg] = more_gmms(gmm,lkp)
% FORMAT gmm = spm_gmm_lib('extras', 'more_gmms', gmm, lkp)
%
% gmm  - Cell with the following format {m,b,W,n}, where there are K
%        Gaussians.
% lkp - [1,K_p] vector partitioning a K GMM into a K_p GMM (K_P>=K). E.g.
%        [1 1 1 2 3 4 5 6 6] means that the first Gaussian will be divided
%        into 3 and the last into 2. The rest will remain the same.
%
% gmm  - Cell with the following format {m,b,W,n}, where there are K_p 
%        Gaussians.
%
% A crude heuristic to replace a single VB Gaussian by a bunch of VB Gaussians.
% If there is only one Gaussian, then it should be the same as the
% original distribution.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

MU0 = gmm{1};
b0  = gmm{2};
W0  = gmm{3};
n0  = gmm{4};

Kb = numel(n0);
K  = numel(lkp);
C  = size(MU0,1);

A0 = bsxfun(@times,W0,reshape(n0,[1 1 Kb])); % E[Lambda]

m  = zeros(C,K);
b  = zeros(1,K);
W  = zeros(C,C,K);
n  = zeros(1,K); % A = nW, W = 1/n*inv(Cov)
mg = ones(1,K);

for k=1:Kb    
    kk = sum(lkp==k);
%     w  = 1./(1 + exp(-(kk - 1)*0.25)) - 0.5;
    w  = 1./(1 + exp(-(kk - 1)*1)) - 0.5;
    mn = MU0(:,k);
    vr = inv(A0(:,:,k));
    
    mn = sqrtm(vr)*randn(C,kk)*w + repmat(mn,[1 kk]);
    vr = vr*(1 - w);
    pr = inv(vr);
    W1 = (1/n0(k))*pr;
    
    m(:,lkp==k)   = mn;
    b(lkp==k)     = b0(k)*kk;
    W(:,:,lkp==k) = repmat(W1,[1 1 kk]);
    n(lkp==k)     = n0(k)*kk;
    
    mg(lkp==k) = 1/kk;
end

gmm{1} = m;
gmm{2} = b;
gmm{3} = W;
gmm{4} = n;
%==========================================================================

% =========================================================================
function [gmm,mg] = collapse_gmms(gmm,lkp)
% FORMAT gmm = spm_gmm_lib('extras', 'collapse', gmm, lkp)
%
% A crude heuristic to replace a bunch of Gaussian with a single Gaussian.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

MU0 = gmm{1};
b0  = gmm{2};
W0  = gmm{3};
n0  = gmm{4};

Kb = max(lkp);
K  = numel(lkp);
C  = size(MU0,1);

A0 = bsxfun(@times,W0,reshape(n0,[1 1 K]));
vr = zeros(size(A0));
for k=1:K
   vr(:,:,k) = inv(A0(:,:,k)); 
end

m  = zeros(C,Kb);
b  = zeros(1,Kb);
W  = zeros(C,C,Kb);
n  = zeros(1,Kb);
mg = ones(1,Kb);

for k=1:Kb                   
    kk = sum(lkp==k);
    
    mn  = mean(MU0(:,lkp == k),2);
    vr1 = mean(vr(:,:,lkp == k),3);
    
    pr = inv(vr1);    
    W1 = (1/n0(k))*pr;
    
    m(:,k)   = mn;    
    W(:,:,k) = W1;
    
    b(k) = mean(b0(lkp == k))/kk;
    n(k) = mean(n0(lkp == k))/kk;
end

gmm{1} = m;
gmm{2} = b;
gmm{3} = W;
gmm{4} = n;
%==========================================================================

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

%==========================================================================
function [m,b,W,n] = get_posteriors(cluster,s)
m = cluster{s}{1}{1};
b = cluster{s}{1}{2};
W = cluster{s}{2}{1};
n = cluster{s}{2}{2};
%==========================================================================