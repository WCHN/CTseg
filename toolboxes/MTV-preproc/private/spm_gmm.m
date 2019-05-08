function varargout = spm_gmm(X, varargin)
%__________________________________________________________________________
%
% Fit a [Bayesian] Gaussian mixture model to observed [weighted] data.
%
% FORMAT [Z,MU,A,PI,...] = spm_gmm(X,...)
% 
% MANDATORY
% ---------
% X - NxP matrix of observed values
% 
% OPTIONAL
% --------
% K  - Number of cluster [0=guess from options, 2 if cannot guess]
% W  - [N]x1 Vector of weights associated with each observation [1]
%
% KEYWORD
% -------
% PropPrior  - 1x[K] vector of Dirichlet priors [0=ML]
%                or NxK matrix of fixed observation-wise proportions.
% GaussPrior - {MU(PxK),b(1x[K]),V([PxP]x[K]),n(1x[K])} Gauss-Wishart prior
%              [{}=ML]
% Prune      - Threshold on proportions to prune uninformative clusters
%              [0=no pruning]
% Missing    - Infer missing data [true]
% Start      - Starting method: METHOD or {METHOD, PRECISION} with
%   METHOD    = ['kmeans'],'linspace','prior','sample','uniform',
%               or provided: MU(PxK) or {MU(PxK),b}
%   PRECISION = [kmeans2gmm or diag(a) with a = (range/(2K))^(-2)],
%               or provided: A([PxP]x[K]) or {V([PxP]x[K]),n}
% KMeans     - Cell of KMeans options [{}].
% IterMax    - Maximum number of EM iterations [1000]
% Tolerance  - Convergence criterion (~ lower bound gain) [1e-4]
% BinWidth   - 1x[P] Bin width (histogram mode: add bits of variance) [0]
% InputDim   - Input space dimension [0=try to guess]
% Verbose    - Verbosity level: [0]= quiet
%                                1 = write (lower bound)
%                                2 = plot (lower bound)
%                                3 = plot more (gmm fit)
%
% OUTPUT
% ------
% Z    - NxK   cluster reponsibility
% MU   - PxK   means                                    [E[MU] if Bayesian]
% A    - PxPxK precision matrices                        [E[A] if Bayesian]
% PI   - 1xK   proportions                              [E[PI] if Bayesian]
% b    - 1xK   posterior mean degrees of freedom              [if Bayesian]
% V    - PxPxK posterior scale matrix                         [if Bayesian]
% n    - 1xK   posterior precision degrees of freedom         [if Bayesian]
% a    - 1xK   posterior Dirichlet                            [if Bayesian]
% X    - NxP   obs and inferred values                           [if infer]
%__________________________________________________________________________
%
% Use a learned mixture to segment an image.
%
% FORMAT [Z, X] = spm_gmm('apply',X,MU,A,PI,...)         > Classical
% FORMAT [Z, X] = spm_gmm('apply',X,{MU,b},{V,n},a,...)  > Bayesian
%
% KEYWORD
% -------
% Missing    - Infer missing data: ['infer']/'remove'
% BinWidth   - 1x[P] Bin width (histogram mode: add bits of variance) [0]
%__________________________________________________________________________
%
% help spm_gmm>Options
% help spm_gmm>TellMeMore
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% TODO
% - Recheck the different Start options

% Convention:
% N: Number of observations
% P: Dimension of observation space
% K: Number of clusters


% -------------------------------------------------------------------------
% Special case: Apply model
% > Here, we use a learned GMM to segment an image
if ischar(X) && strcmpi(X, 'apply')
    [varargout{1:nargout}] = gmm_apply(varargin{:});
    return
end

% -------------------------------------------------------------------------
% Parse inputs
p = inputParser;
p.FunctionName = 'spm_gmm';
p.addRequired('X',                          @isnumeric);
p.addOptional('K',           0,             @(X) isscalar(X) && isnumeric(X));
p.addOptional('W',           1,             @isnumeric);
p.addParameter('PropPrior',  0,             @isnumeric);
p.addParameter('GaussPrior', {},            @iscell);
p.addParameter('Prune',      0,             @(X) (isscalar(X) && isnumeric(X)) || islogical(X));
p.addParameter('Missing',    true,          @islogical);
p.addParameter('Start',      'kmeans',      @(X) ischar(X) || isnumeric(X));
p.addParameter('KMeans',     {},            @iscell);
p.addParameter('IterMax',    1000,          @(X) isscalar(X) && isnumeric(X));
p.addParameter('Tolerance',  1e-4,          @(X) isscalar(X) && isnumeric(X));
p.addParameter('BinWidth',   0,             @isnumeric);
p.addParameter('InputDim',   0,             @(X) isscalar(X) && isnumeric(X));
p.addParameter('Verbose',    0,             @(X) isscalar(X) && (isnumeric(X) || islogical(X)));
p.parse(X, varargin{:});
W          = p.Results.W;
K          = p.Results.K;
P          = p.Results.InputDim;
E          = p.Results.BinWidth;
Start      = p.Results.Start;
KMeans     = p.Results.KMeans;
PropPrior  = p.Results.PropPrior;
GaussPrior = p.Results.GaussPrior;
Prune      = p.Results.Prune;
Missing    = p.Results.Missing;

% -------------------------------------------------------------------------
% A bit of formatting
if ~iscell(Start)
    Start = {Start};
end
if ~iscell(GaussPrior)
    GaussPrior = {GaussPrior};
end
if islogical(Prune)
    Prune = 1e-7 * Prune;
end

% -------------------------------------------------------------------------
% Guess dimension/clusters from provided initial values
[K,P,Start]      = dimFromStart(K,P,Start);
[K,P,GaussPrior] = dimFromGaussPrior(K,P,GaussPrior);
[K,PropPrior]    = dimFromPropPrior(K,PropPrior);
[P,dimX]         = dimFromObservations(P, X);
% ---
% default value
if K == 0, K = 2; end
        

% -------------------------------------------------------------------------
% Proportion / Dirichlet
if size(PropPrior, 1) > 1
    % Class probability map (fixed proportions)
    PI              = reshape(PropPrior, [], K);
    logPI           = log(max(PI,eps));
    a0              = zeros(1,K, 'like', PI);
else
    % Dirichlet prior
    a0              = PropPrior(:)';
    if numel(a0) < K
        a0          = padarray(a0, [0 K - numel(a0)], 'replicate', 'post');
    end
    PI              = [];
    logPI           = [];
end

% -------------------------------------------------------------------------
% Reshape X (observations)
if dimX(1) == 1 && numel(dimX) == 2
    % row-vector case
    latX = dimX;
else
    latX = dimX(1:end-1);
end
X  = reshape(X, [], P);
N  = size(X, 1); % Number of observations
N0 = N;          % Original number of observations


% -------------------------------------------------------------------------
% Reshape W (weights)
W = W(:);

% -------------------------------------------------------------------------
% Prepare missing data stuff (code image, mask, ...)
if ~any(any(isnan(X)))
    Missing = false;
end
if Missing
    % Deal with missing data
    code      = spm_gmm_lib('obs2code', X);     % Code image
    code_list = unique(code);                   % List of codes
    missmsk   = [];
else
    % Compute mask of removed rows
    missmsk   = any(isnan(X),2);
    code      = spm_gmm_lib('double2int', (2^P-1) * ones(sum(~missmsk),1));
    code_list = 2^P-1;
end
% Discard rows with missing values
if ~isempty(missmsk)
    X         = X(~missmsk,:);
    if size(W,1) > 1
        W     = W(~missmsk);
    end
    if size(PI,1) > 1
        PI0   = PI(missmsk,:);
        PI    = PI(~missmsk,:);
    end
    N         = sum(~missmsk);
end
missmsk = find(missmsk); % saves a bit of memory

% -------------------------------------------------------------------------
% "Bin" variance
% > When we work with histograms, a bit of variance is lost due to the
% binning. Here, we assume some kind of uniform distribution inside the bin
% and consequently add the corresponding variance to the 2nd order moment.
if numel(E) < P
    E = padarray(E, [0 P - numel(E)], 'replicate', 'post');
end
E = (E.^2)/12;

% -------------------------------------------------------------------------
% Initialise Gauss-Wishart prior
[MU0,b0,V0,n0] = initialise_prior(GaussPrior, K, P);
pr = struct('MU', MU0, 'b', b0, 'V', V0, 'n', n0);
    
% -------------------------------------------------------------------------
% Initialise mixture
[~, MU, b, A, V, n, PI00,logPI00] = start(Start, X, W, K, a0, pr, KMeans);
if isempty(PI)
    PI    = PI00;
    logPI = logPI00;
end
clear PI00 logPI00

% -------------------------------------------------------------------------
% Default prior mean/precision if needed
if isempty(MU0) && sum(b0) > 0
    MU0 = MU;
end
if isempty(V0) && sum(n0) > 0
    if sum(n) > 0
        V0 = V;
    else
        V0 = bsxfun(@rdivide, A, reshape(n0, [1 1 K]));
    end
end

% -------------------------------------------------------------------------
% Initialise posterior
if sum(b) == 0
    b = b0;
end
if sum(n) == 0
    n = n0;
    if sum(n) > 0, V = bsxfun(@rdivide,A,reshape(n,[1 1 K])); end
end
if sum(b) > 0, mean = {MU, b};
else,          mean = MU;        end
if sum(n) > 0, prec = {V, n};
else,          prec = {A};       end

% -------------------------------------------------------------------------
% MAIN LOOP
[Z,cluster,prop] = spm_gmm_loop({X,W}, {mean,prec}, ...
    {'LogProp', logPI, 'Prop', PI}, ...
    'GaussPrior',     {MU0, b0, V0, n0}, ...
    'PropPrior',      a0, ...
    'Missing',        Missing, ...
    'MissingCode',    {code,code_list}, ...
    'IterMax',        p.Results.IterMax, ...
    'Tolerance',      p.Results.Tolerance, ...
    'SubIterMax',     p.Results.IterMax, ...
    'SubTolerance',   p.Results.Tolerance, ...
    'BinUncertainty', E, ...
    'Verbose',        p.Results.Verbose);

MU = cluster.MU;
b  = cluster.b;
A  = cluster.A;
V  = cluster.V;
n  = cluster.n;
PI = prop.Prop;
a  = prop.Dir;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Infer missing values
if nargout >= 9 && Missing
    X = spm_gmm_lib('InferMissing', X, Z, {MU,A}, {code,code_list});
end

% -------------------------------------------------------------------------
% Prune clusters
if Prune > 0
    [K,PI,Z,MU,b,A,V,n] = prune(Prune, PI, Z, {MU,b}, {A,V,n});
end

% -------------------------------------------------------------------------
% Replace discarded missing values
if sum(missmsk) > 0
    present = ones(N0, 1, 'logical');
    present(missmsk) = false;
    clear missing
    
    Z = expand(Z, present, N0, K, 0);
    if size(PI,1) > 1 && nargout >= 4
        PI = expand(PI, present, N0, K, PI0);
    end
    if nargout >= 9
        X = expand(X, present, N0, P, NaN);
    end
end
    
% -------------------------------------------------------------------------
% Reshape everything (use input lattice)
Z = reshape(Z, [latX K]);
if size(PI,1) > 1 && nargout >= 4
    PI = reshape(PI, [latX K]);
end
if nargout >= 9
    X = reshape(X, [latX P]);
end

% -------------------------------------------------------------------------
% Push results in output object
if nargout >= 1, varargout{1} = Z;
if nargout >= 2, varargout{2} = MU;
if nargout >= 3, varargout{3} = A;
if nargout >= 4, varargout{4} = PI;
if nargout >= 5, varargout{5} = b;
if nargout >= 6, varargout{6} = V;
if nargout >= 7, varargout{7} = n;
if nargout >= 8, varargout{8} = a;
if nargout >= 9, varargout{9} = X;
end;end;end;end;end;end;end;end;end


% =========================================================================
function [Z,X] = gmm_apply(X, mean, prec, prop, varargin)
% FORMAT [Z,X] = gmm_apply(X, MU, A, PI, ...)
% FORMAT [Z,X] = gmm_apply(X, {MU,b}, {V,n}, a, ...)

% -------------------------------------------------------------------------
% Parse inputs
p = inputParser;
p.FunctionName = 'gmm_apply';
p.addRequired('X',                     @isnumeric);
p.addRequired('mean',                  @(X) isnumeric(X) || iscell(X));
p.addRequired('precision',             @(X) isnumeric(X) || iscell(X));
p.addRequired('prop',                  @(X) isnumeric(X) || iscell(X));
p.addParameter('Missing',        true, @islogical);
p.addParameter('BinWidth', 0,    @isnumeric);
p.parse(X, mean, prec, prop, varargin{:});
E          = p.Results.BinWidth;
Missing    = p.Results.Missing;

MU = [];
b  = [];
A  = [];
V  = [];
N  = [];
PI = [];
a  = [];
logPI = [];

% -------------------------------------------------------------------------
% Read arguments
if ~iscell(mean)
    MU = mean;
elseif ~isempty(mean)
    MU = mean{1};
    if numel(mean) > 1
        b = mean{2};
    end
end
if ~iscell(prec)
    A = prec;
elseif ~isempty(prec)
    A = prec{1};
    if numel(prec) > 1
        n = prec{2};
        if sum(b) > 0
            V = A;
            A = bsxfun(@times, V, reshape(n,1,1,[]));
            prec = {V n};
        end
    end
end
if ~iscell(prop)
    PI = prop;
elseif ~isempty(prop)
    PI = prop{1};
end

% -------------------------------------------------------------------------
% Dimensions
P = size(MU,1);
K = size(MU,2);

% -------------------------------------------------------------------------
% Proportions/Dirichlet
if numel(PI) <= K
    PI = PI(:)';
    PI = padarray(PI, [0 K - numel(PI)], 'replicate', 'post');
    
    if abs(sum(PI)-1) > eps('single')
        % Dirichlet prior
        a     = PI;
        logPI = psi(a) - psi(sum(a));
    else
        logPI = log(max(PI,eps));
    end
else
    logPI = log(max(PI,eps));
end
clear PI

% -------------------------------------------------------------------------
% Reshape X (observations)
dimX = size(X);
if P == 1
    latX = dimX;
    if latX(2)==1
        latX = latX(1);
    end
else
    latX = dimX(1:end-1);
end
X = reshape(X, [], P);
N  = size(X, 1); % Number of observations
N0 = N;          % Original number of observations

% -------------------------------------------------------------------------
% Prepare missing data stuff (code image, mask, ...)
if ~any(any(isnan(X)))
    Missing = false;
end
if Missing
    % Deal with missing data
    code      = spm_gmm_lib('obs2code', X);     % Code image
    code_list = unique(code);                   % List of codes
    missmsk   = [];
else
    % Compute mask of removed rows
    missmsk   = any(isnan(X),2);
    code      = spm_gmm_lib('double2int', (2^P-1) * ones(sum(~missmsk),1));
    code_list = 2^P-1;
end
% Discard rows with missing values
if ~isempty(missmsk)
    X         = X(~missmsk,:);
    if size(logPI,1) > 1
        logPI = logPI(~missmsk,:);
    end
    N         = sum(~missmsk);
end
missmsk = find(missmsk); % saves a bit of memory

% -------------------------------------------------------------------------
% "Bin" variance
% > When we work with histograms, a bit of variance is lost due to the
% binning. Here, we assume some kind of uniform distribution inside the bin
% and consequently add the corresponding variance to the 2nd order moment.
if numel(E) < P
    E = padarray(E, [0 P - numel(E)], 'replicate', 'post');
end
E = (E.^2)/12;

% -------------------------------------------------------------------------
% Compute marginal log-likelihood
if Missing
    const = spm_gmm_lib('Const', mean, prec, code_list);
else
    const  = spm_gmm_lib('Const', mean, prec);
end
logpX = spm_gmm_lib('Marginal', X, [{MU} prec], const, {code,code_list}, E);

% -------------------------------------------------------------------------
% Compute responsibilities
Z = spm_gmm_lib('Responsibility', logpX, logPI);
clear logpX logPI
 
% -------------------------------------------------------------------------
% Infer missing values
if nargout >= 2 && Missing
    X = spm_gmm_lib('InferMissing', X, Z, {MU,A}, {code,code_list});
end
 
% -------------------------------------------------------------------------
% Replace discarded missing values
if sum(missmsk) > 0
    present = ones(N0, 1, 'logical');
    present(missmsk) = false;
    clear missing
    
    Z = expand(Z, present, N0, K, 0);
    if nargout >= 2
        X = expand(X, present, N0, P, NaN);
    end
end
    
% -------------------------------------------------------------------------
% Reshape everything (use input lattice)
Z = reshape(Z, [latX K]);
if nargout >= 2
    X = reshape(X, [latX P]);
end

% =========================================================================
function TellMeMore
% _________________________________________________________________________
%
% Gaussian mixture model (GMM)
% ----------------------------
% The  Gaussian  Mixture   relies  on  a  generative  model  of  the  data.
% Each  observation is assumed  to stem  from one of  K clusters,  and each
% cluster possesses a Gaussian density.
%
% Classical GMM
% -------------
% With  the  most  basic  (and  well  known)  model, we  look  for  maximum
% likelihood  values  for  the  model  parameters,  which are  the mean and
% precision matrix of each cluster: {Mu, A}_k = argmax p({x}_n | {Mu, A}_k)
% To compute  this probabilitity,  we need to integrate  over all  possible
% cluster  responsibility  (to which cluster belongs  a given observation),
% which are unknown:
%   p({x}_n | {Mu, A}_k) = int p({x}_n | {z}_n, {Mu, A}_k) p({z}_n) d{z}_n
% Since  this integral is intractable,  we use the Expectation-Maximisation
% algorithm:  we alternate between computing  the posterior distribution of
% responsibilities (given known means and precision matrices)  and updating 
% the mean and precisions (given the known posterior).
%
% Bayesian GMM
% ------------
% When we have some idea  about how these mean and precision look like,  we
% can take it  into account in the form  of Bayesian beliefs,  i.e.,  prior
% probability distributions  over the parameters we want to estimate.  What
% we  are now  looking  for  are the  posterior  distributions  (given some
% observed data) of these parameters:
%   p({Mu, A}_k | {x}_n) = p({x}_n | {Mu, A}_k) p({Mu, A}_k) / p({x}_n)
% To  make everything  tractable,  these prior  beliefs are  chosen  to  be
% conjugate priors. It is not very important to know what it means,  except
% that it makes computing the posterior probabilities easier.  In our case,
% we can use a  Gaussian prior distribution for the means,  a Wishart prior
% distribution  for   the  precision  matrices   (we  talk  of  Gauss-Prior
% distirbution  when they  are combined)  and a  Dirichlet distribution for
% clusters' proportion.
% Despite all that,  we still cannot compute  these posteriors.  We make an
% additional assumption,  which is that there is  some sort of independence
% between  the  parameters to estimate;  we say that  the posterior  can be
% factorised:
%   q({Mu, A, Pi}_k, {z}_n) = q({Mu, A}_k) q({Pi}_k) q({z}_n)
% We  can then  use a technique  called variational  Bayesian inference  to
% estimate, in turn, these distributions.
%
% Histogram GMM
% -------------
% To speed up computation,  we sometimes prefer  to use an histogram  (that
% is,  binned observations)  rather than  the full set  of observations.  A 
% first way of doing this,  is by assuming that each bin centre corresponds 
% to several identical observations (the bin count). In other words, we now
% have weighted observations.
% However, using weighted observations makes estimating the precisions less
% robust,  as we artificially  reduce the  variance by  assigning different
% values to the same bin centre.  To lower this effect,  we can assume that
% each observation  is  actually  a  distribution, i.e.,  there is a bit of
% uncertainty about the "true" value.  When computing the model, we want to
% integrate over all possible values.
% Here,  we assume that these distributions  are uniform in each bin.  This
% makes the implementation  easy and efficient:  the expected value of each
% observation is the bin centre, and their variance is (w^2)/12, where w is
% the  bin width.  This consists of adding a bit of  jitter when  computing
% second order statistics to update the precision matrices.
%
% Missing values
% --------------
% In the classical GMM case, in the presence of partial observations  (in a
% given voxel, some modalities might be missing), it is easy to compute the
% conditional  likelihood  of a  given  voxel,  but computing  ML mean  and
% precision estimates is more tricky.
% With one Gaussian,  an EM scheme can be designed  to obtain ML parameters
% by  alternating   between  inferring  missing  values and   updating  the
% parameters.   This scheme can be  extended to the mixture case,  in which
% case the EM algorithm relies on a joint posterior distribution over class
% repsonsibilities and missing values.
% _________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% =========================================================================
function Options
% _________________________________________________________________________
%
% FORMAT [Z,MU,A,PI,...] = spm_gmm(X,...)
% 
% MANDATORY
% ---------
% X is a set of multidimensional observations. It usually takes the form of
%   a  matrix,  in which  case  its first  dimension (N)  is the number  of
%   observations and its second dimension (P) is their dimension.  However, 
%   it is possible  to provide a multidimensional array,  in which case the
%   dimension P  must be provided  in some way.  This might be  through the 
%   size of user-provided  starting estimates or priors,  or directly using 
%   the 'InputDim' option. If the 'Missing' option is activated, X might
%   contain missing values (NaN), which will be inferred by the model.  If
%   the 'Missing'  option is  not activated,  all rows  containing  missing
%   values will be excuded.
% 
% OPTIONAL
% --------
% K is the  number of  clusters in the  mixture.  If not provided,  we will
%   first try to guess it from user-provided options (starting estimates or 
%   priors). If this is not possible, the default value is 2.
%
% W is a  vector of  weights  associated  with each  observation.  Usually,
%   weights are used when the input dataset is an histogram. Suitable X and
%   W arrays  can be obtained  with the  spm_imbasics('hist') function.  In
%   this  case,  it is  advised to  also use the  'BinWidth' option,  which
%   prevents  variances to collapse  due to the pooling of values.  If W is
%   not provided, all observations have weight 1.
%
% KEYWORD
% -------
% PropPrior  may take two forms (+ the default):
%            0) It can be empty (the default),  in which case the algorithm
%            searches for maximum-likelihood proportion.
%            1) It  can  be  a  vector of  concentration  parameters  for a 
%            Dirichlet prior.  Dirichlet distributions are conjugate priors 
%            for proportions,  e.g.  Categorical parameters.  In this case, 
%            it  must  consist  of   K  strictly  positive  values.   Their 
%            normalised value is the prior expected proportion, while their 
%            sum  corresponds  to the precision  of the  distribution  (the 
%            larger the precision, the stonger the prior).
%            2) It can be a  matrix  (or multidimensional array)  of  fixed
%            observation-wise    proportions.    We   often   talk    about
%            non-stationary  class  proportions.  In  this  case,  it  must
%            contain N (the number of observations) times K  (the number of
%            clusters) elements. Elements must sum to one along the cluster
%            dimension.
%
% GaussPrior {Mu(PxK),b(1x[K]),V([PxP]x[K]),n(1x[K])}
%            The Gauss-Wishart  distribution is a  conjugate prior  for the 
%            parameters  (mean  and  precision matrix)  of  a  multivariate
%            Gaussian distribution.  Its  parameters are  Mu,  the expected
%            mean,  b,  the prior degrees of freedom for the mean,  V,  the 
%            scale  matrix  and n,  the prior degrees  of freedom  for  the
%            precision.  The expected precision matrix is n*V.  This option
%            allows to defines these parameters for each cluster: it should
%            be a cell  of arrays with  the dimensions  stated above.  Note
%            that some  parameters might be left empty,  in which case they
%            will be  automatically determined.  Typically,  one could only
%            provide the  degrees of freedom  (b and n),  in which case the
%            starting  estimates will be  used as  prior expected means and
%            precisions.  Also, dimensions that are written in brackets are 
%            automatically expanded if needed.
%            By  default,   this  option   is  empty,   and  the  algorithm
%            searches for maximum-likelihood parameters.
%
% Prune      contains a  threshold  for  the  final  estimated proportions. 
%            Classes  that are under this threshold  will be considered  as
%            non-informative and pruned out. By default, the threshold is 0
%            and no pruning is performed.
%
% Missing     'infer': (default) missing values are inferred by making  use
%                      of  some properties of  the  Gaussian  distribution. 
%                      However, the fit is  slower due to  this  inferrence 
%                      scheme.
%            'remove': all rows  with missing values are exclude,  and  the
%                      fit   is   performed    without   inferrence.    For 
%                      computational efficiency, when no values are missing 
%                      from the input set  (i.e., there are no NaNs),  this 
%                      option is activated by default.
%
% Start      Method used to select starting estimates:
%              'kmeans': A K-means  clustering is  used to  obtain a  first
%                        classification of the observations.  Centroids are
%                        used  to  initialise   the  means,   intra-cluster
%                        co-variance   is  used   to  initialise  precision 
%                        matrices  (unless initial  precision matrices  are
%                        user-provided)   and  cluster  size   is  used  to 
%                        initialise proportions. Options can be provided to
%                        the K-means algorithm using the 'KMeans' option.
%            'linspace': Centroids  are  chosen so  that they are  linearly
%                        spaced along the input range of values.  Precision
%                        matrices are set as explained below.
%               'prior': The prior  expected means  and  precision matrices
%                        are used as initial estimates.
%              'sample': Random samples are selected from the input dataset
%                        and used as initial means.  Precision matrices are 
%                        set as explained below.
%             'uniform': Random values are uniformly sampled from the input
%                        range  of  values   and  used  as  initial  means. 
%                        Precision matrices are set as explained below.
%               MU(PxK): Initial means are user-provided.
%            The  method  can  be  provided along  with  initial  precision
%            matrices,  in  a  cell:  {METHOD, A([PxP]x[K])}.  By  default, 
%            the  initial  precision matrix  is  common to all classes  and 
%            chosen so that the input range is well covered:
%            -> A = diag(a) with a(p) = (range(p)/(2K))^(-2)
%
% KMeans     is a cell of options for the K-means algorithm.
%
% IterMax    is the maximum number of EM iterations. Default is 1000.
%
% Tolerance  is the  convergence  threshold  below  which the algorithm  is 
%            stopped. The convergence criterion is a normalised lower bound
%            gain (lcur-lprev)(lmax-lmin). The default tolerance is 1e-4.
%
% BinWidth   is a vector of bin widths.  It is useful when  the input is an
%            histogram  (e.g.  weighted  observations)  to  regularise  the
%            variance estimation.  If only  one width  is  provided,  it is
%            automatically expanded to all dimensions.
%
% InputDim   Number of dimensions in the input space.  If not provided,  we
%            will try  to infer  it from the  input  arrays, e.g.  starting
%            estimates, last dimension of the input array, etc.
%
% Verbose    0: Quiet mode. Nothing is written or plotted.
%            1: Verbose  mode.  The  lower  bound  is  written  after  each
%               iteration. This slows down significantly the algorithm.
%            2: Graphical mode. The evolution of the lower bound and its
%               different components is plotted.
%            3: Graphical mode +.  A representation  of the  fit  (marginal
%               distribution,  joint  2D  distributions,  proportions)   is 
%               plotted.  This slows things down  dramatically  and  should 
%               only be used for education or debugging purpose.
% _________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% =========================================================================
function [K,P,Start] = dimFromStart(K, P, Start)
% FORMAT [K,P,Start] = dimFromStart(K, P, Start)
% K - Number of clusters (previous guess)
% K - Number of dimensions (previous guess)
% Start - Cell of "starting estimates" options
%
% Guess input dimension and number of clusters from starting estimates.

% ---
% Guess from starting mean
if numel(Start) > 1 && ~ischar(Start{2})
    if size(Start{2}, 1) > 1
        if P == 0
            P = size(Start{2}, 1);
        elseif P ~= size(Start{2}, 1)
            warning(['Input space dimension does not agree with starting ' ...
                     'precision: %d vs. %d'], P, size(Start{2}, 1));
            Start{2} = Start{2}(1:P,1:P,:);
        end
    end
    if size(Start{2}, 3) > 1
        if K == 0
            K = size(Start{2}, 3);
        elseif K ~= size(Start{2}, 3)
            warning(['Number of clusters does not agree with starting ' ...
                     'precision: %d vs. %d'], K, size(Start{2}, 3));
            Start{2} = Start{2}(:,:,1:K);
        end
    end
end
if ~ischar(Start{1})
    if size(Start{1}, 1) > 0
        if P == 0
            P = size(Start{1}, 1);
        elseif P ~= size(Start{1}, 1)
            warning(['Input space dimension does not agree with starting ' ...
                     'mean: %d vs. %d'], P, size(Start{1}, 1));
            Start{1} = Start{1}(1:P,:);
        end
    end
    if size(Start{1}, 2) > 0
        if K == 0
            K = size(Start{1}, 2);
        elseif K ~= size(Start{1}, 2)
            warning(['Number of clusters does not agree with starting ' ...
                     'mean: %d vs. %d'], K, size(Start{1}, 2));
            Start{1} = Start{1}(:,1:K);
        end
    end
end

% ---
% Guess from starting precision
if numel(Start) > 1 && ~ischar(Start{2})
    if size(Start{2}, 1) > 1
        if P == 0
            P = size(Start{2}, 1);
        elseif P ~= size(Start{2}, 1)
            warning(['Input space dimension does not agree with starting ' ...
                     'precision: %d vs. %d'], P, size(Start{2}, 1));
            Start{2} = Start{2}(1:P,1:P,:);
        end
    end
    if size(Start{2}, 3) > 1
        if K == 0
            K = size(Start{2}, 3);
        elseif K ~= size(Start{2}, 3)
            warning(['Number of clusters does not agree with starting ' ...
                     'precision: %d vs. %d'], K, size(Start{2}, 3));
            Start{2} = Start{2}(:,:,1:K);
        end
    end
end

% =========================================================================
function [K,P,GaussPrior] = dimFromGaussPrior(K, P, GaussPrior)
% FORMAT [K,P,GaussPrior] = dimFromGaussPrior(K, P, GaussPrior)
% K          - Number of clusters (previous guess)
% K          - Number of dimensions (previous guess)
% GaussPrior - Cell of "prior" options
%
% Guess input dimension and number of clusters from Gauss-Wishart prior

% ---
% Guess from prior mean
if numel(GaussPrior) > 1
    if size(GaussPrior{1}, 1) > 0
        if P == 0
            P = size(GaussPrior{1}, 1);
        elseif P ~= size(GaussPrior{1}, 1)
            warning(['Input space dimension does not agree with prior ' ...
                     'mean: %d vs. %d'], P, size(GaussPrior{1}, 1));
            GaussPrior{1} = GaussPrior{1}(1:P,:);
        end
    end
    if size(GaussPrior{1}, 2) > 0
        if K == 0
            K = size(GaussPrior{1}, 2);
        elseif K ~= size(GaussPrior{1}, 2)
            warning(['Number of clusters does not agree with prior ' ...
                     'mean: %d vs. %d'], K, size(GaussPrior{1}, 2));
            GaussPrior{1} = GaussPrior{1}(:,1:K);
        end
    end
end
% ---
% Guess from prior precision
if numel(GaussPrior) > 3
    if size(GaussPrior{2}, 1) > 1
        if P == 0
            P = size(GaussPrior{2}, 1);
        elseif P ~= size(GaussPrior{2}, 1)
            warning(['Input space dimension does not agree with prior ' ...
                     'precision: %d vs. %d'], P, size(GaussPrior{2}, 1));
            GaussPrior{2} = GaussPrior{2}(1:P,1:P,:);
        end
    end
    if size(GaussPrior{2}, 3) > 1
        if K == 0
            K = size(GaussPrior{2}, 3);
        elseif K ~= size(GaussPrior{2}, 3)
            warning(['Number of clusters does not agree with prior ' ...
                     'precision: %d vs. %d'], K, size(GaussPrior{2}, 3));
            GaussPrior{2} = GaussPrior{2}(:,:,1:K);
        end
    end
end


% =========================================================================
function [K,PropPrior] = dimFromPropPrior(K, PropPrior)
% Guess number of clusters from Dirichlet prior

% ---
% Guess from prior proportion
if size(PropPrior, 2) > 1
    if K == 0
        K = size(PropPrior, 2);
    elseif K ~= size(PropPrior, 2)
        warning(['Number of clusters does not agree with prior ' ...
                 'proportion: %d vs. %d'], K, size(PropPrior, 2));
        PropPrior = PropPrior(:,1:K);
    end
end

% =========================================================================
function [P,dimX] = dimFromObservations(P, X)
dimX = size(X);
if P == 0
    if numel(dimX) == 2
        if dimX(1) == 1
            % row-vector case
            P = dimX(1);
        else
            % matrix case
            P = dimX(2);
        end
    else
        % N-array case
        dim = size(X);
        P = dim(end);
    end
end

% =========================================================================
function varargout = prune(threshold, PI, Z, mean, prec)
% FORMAT prune(threshold, PI, Z, {MU,b}, {A,V,n})
%
% Remove classes with proportion <= threshold

kept = sum(PI,1)/sum(PI(:)) >= threshold;
K    = sum(kept);

Z    = Z(:,kept);
PI   = PI(:,kept);

if ~iscell(mean)
    mean = {mean(:,kept)};
else
    if numel(mean) >= 1
        mean{1} = mean{1}(:,kept);
        if numel(mean) >= 2 && sum(mean{2}) > 0
            mean{2} = mean{2}(kept);
        end
    end
end

if ~iscell(prec)
    prec = {prec(:,kept)};
else
    if numel(prec) >= 1
        if ~isempty(prec{1})
            prec{1} = prec{1}(:,:,kept);
        end
        if numel(prec) >= 2
            if ~isempty(prec{2})
                prec{2} = prec{2}(:,:,kept);
            end
            if numel(prec) >= 3 && sum(prec{3}) > 0
                prec{3} = prec{3}(kept);
            end
        end
    end
end

varargout = [{K} {PI} {Z} mean prec];

% =========================================================================
function X = expand(X, msk, N, P, val)
X1 = X;
val = cast(val, 'like', X1);
switch val
    case 0
        X = zeros(N, P, 'like', X1);
    case 1
        X = ones(N, P, 'like', X1);
    case Inf
        X = Inf(N, P, 'like', X1);
    otherwise
        if numel(val) == 1
            X = val * ones(N, P, 'like', X1);
        elseif numel(val) == P
            X = repmat(val, [N 1]);
        end
end
X(msk,:) = X1; clear X1

% =========================================================================
function [Z,MU,b,A,V,n,PI,logPI] = start(method, X, W, K, a0, pr, kmeans)
% FORMAT C = start(method, X, W, K, a0, pr, kmeans)
%
% method - Method to use to select starting centroids
%               'kmeans', 'sample', 'uniform' or provided matrix
% X      - Vector of NxP observed values
% W      - Vector of Nx1 weights
% K      - Number of clusters
% a0     - 1xK dirichlet priors (or empty)
% pr     - Structure of Gauss-Wishart prior parameters
% kmeans - Options for spm_kmeans
%
% MU - PxK   means
% A  - PxPxK precision matrices
% PI - 1xK   proportions
%
% Compute starting estimates

if ~iscell(method)
    method = {method};
end

b = [];
n = [];
V = [];

if ischar(method{1})
    switch method{1}

        case 'kmeans'
        % Use kmeans to produce a first clustering of the data
            P = size(X,2);
            N = size(X,1);
            % Get labelling and centroids from K-means
            [L,MU] = spm_kmeans(X, K, W, kmeans{:});
            MU = MU';
            % Convert to "responsibility"
            Z = zeros(N,K);
            for k=1:K
                Z(:,k) = (L == k) + eps;
            end
            Z(all(isnan(X),2),:) = 1/K;
            clear L
            Z = bsxfun(@rdivide, Z, sum(Z,2));
            % Compute proportion from labelling
            PI = zeros(1,K);
            for k=1:K
                PI(k) = sum(Z(:,k)) + a0(k);
            end
            if sum(a0) > 0
                logPI = psi(PI) - psi(sum(PI));
                PI    = PI ./ sum(PI);
            else
                PI    = max(PI,eps);
                PI    = PI ./ sum(PI);
                logPI = log(PI);
            end
            if numel(method) == 1
                % Compute precision from intra-class sample variance
                % 1) SuffStat2 = sum_i { Wi * (Xi-MU)*(Xi-MU)' }
                X2 = bsxfun(@minus,X,reshape(MU, [1 P K]));
                A  = zeros(N,P,P,K);
                for c1=1:P
                    A(:,c1,c1,:) = reshape(X2(:,c1,:).^2, [N 1 1 K]);
                    for c2=c1+1:P
                        A(:,c1,c2,:) = reshape(X2(:,c1,:) .* X2(:,c2,:), [N 1 1 K]);
                        A(:,c2,c1,:) = A(:,c1,c2,:);
                    end
                end
                clear X2
                A = sum(bsxfun(@times, A, reshape(bsxfun(@times, Z, W), [N 1 1 K])), 1, 'omitnan');
                A = reshape(A, [P P K]);
                % 2) SuffStat0 = sum_i { Wi * Present*Present' }
                X   = ~isnan(X);
                SUM = zeros(N,P,P);
                for c1=1:P
                    SUM(:,c1,c1) = X(:,c1,:).^2;
                    for c2=c1+1:P
                        SUM(:,c1,c2) = X(:,c1) .* X(:,c2);
                        SUM(:,c2,c1) = SUM(:,c1,c2);
                    end
                end
                SUM = sum(bsxfun(@times, SUM, reshape(bsxfun(@times, Z, W), [N 1 1 K])), 1, 'omitnan');
                SUM = reshape(SUM, [P P K]);
                % 3) Normalise and invert
                A = A ./ SUM; clear SUM X
                for k=1:K
                    A(:,:,k) = A(:,:,k) + A(:,:,k)'; % Ensure symmetric
                    A(:,:,k) = spm_matcomp('Inv', A(:,:,k));
                end
            end

        case 'prior'
        % Use prior expected value
            if isempty(pr.MU) || isempty(pr.V) || isempty(pr.n)
                error('To initialise from prior, a full prior must be provided')
            end
            MU = pr.MU;
            b  = pr.b;
            V  = pr.V;
            n  = pr.n;
            A  = bsxfun(@times, V, reshape(n, [1 1 K]));

        case 'sample'
        % Sample uniform
        % All centroids are selected at random from the observed values.
        % They are all unique (to avoid starting with several identical
        % centroids)
            X   = X(~(any(W==0,2)|any(isnan(X),2)),:); % Remove all rows w/ NaNs or w/o obs
            i   = randperm(size(X,1));
            i   = i(1:K);
            MU = X(i,:)';

        case 'uniform'
        % Range uniform
        % All centroids are selected at random from the continuous range of 
        % observed values.
        % They are all unique (to avoid starting with several identical
        % centroids)
            minval = min(X, [], 1, 'omitnan');
            maxval = max(X, [], 1, 'omitnan');
            MU = rand(K,size(X,2));
            MU = bsxfun(@times, MU, maxval - minval) + minval;
            MU = MU';

        case 'linspace'
        % Range uniform
        % All centroids are selected at random from the continuous range of 
        % observed values.
        % They are all unique (to avoid starting with several identical
        % centroids)
            minval = min(X, [], 1, 'omitnan');
            maxval = max(X, [], 1, 'omitnan');
            MU = rand(size(X,2),K);
            for p=1:size(X,2)
                ticks = linspace(minval(p), maxval(p), 2*K+1);
                MU(p,:) = ticks(2:2:end-1);
            end
    end
else
   % Provided
    if ~iscell(method{1})
        MU = method{1};
    else
        if numel(method{1}) >= 1
            MU = method{1}{1};
            if numel(method{1}) >= 2
                b = method{1}{2};
            end
        end
    end         
end

% Precision matrix
% - If nothing is provided, we split the range in 2K pieces and set the 
%   standard deviation (sigma) to that value. This way, half of each 
%   Gaussian covers (1/K)th of the range of values.
% - If scalars,are provided, we use them as precisions in diagonal
%   precision matrices.
% - If matrices are provided, we just copy them.
if numel(method) >= 2 || ~any(strcmpi(method{1}, {'kmeans','prior'}))
    if numel(method) >= 2 && (~iscell(method{2}) || numel(method{2}) >= 1)
        if ~iscell(method{2})
            A = method{2};
        else
            A = method{2}{1};
            if numel(method{2}) >= 2
                n = method{2}{2};
                if sum(n) > 0
                    V = A;
                    A = bsxfun(@rdivide, V, n);
                end
            end
        end
    else
        minval = min(X, [], 1, 'omitnan');
        maxval = max(X, [], 1, 'omitnan');
        A = diag(((maxval - minval)./(2*K)).^(-2));
    end
    if size(A,3) < K
        A = padarray(A, [0 0 K-size(A,3)], 'replicate', 'post');
    end
    if size(A,1) == 1
        A = bsxfun(@times, A, eye(size(X,2)));
    end
end
% Proportion
% - If a Dirichlet prior is provided, we use its expected value
% - Else, we use 1/K
if ~any(strcmpi(method{1}, {'kmeans'}))
    if sum(a0) > 0
        PI    = a0;
        logPI = psi(PI) - psi(sum(PI));
        PI    = PI ./ sum(PI);
    else
        PI    = 1/K * ones(1,K);
        PI    = PI ./ sum(PI);
        logPI = log(PI);
    end
end
% Responsibility
% - We use the expected proportion.
if ~any(strcmpi(method{1}, {'kmeans'}))
    Z = repmat(PI(:)', [size(X,1) 1]);
    allnan = all(isnan(X), 2);
    Z(allnan,:) = repmat([0.9 repmat(0.1/(K-1), [1 K-1])], [sum(allnan) 1]);
end

% =========================================================================
function [MU,b,V,n] = initialise_prior(pr0, K, P)
% FORMAT pr = initialise_prior(pr0)
%
% pr0 - Cell of paramters {MU, b, V, n}
% K   - Number of classes
% P   - Number of channels
%
% MU  - Mean
% b   - Mean degrees of freedom
% V   - Scale matrix (E[A] = nV)
% n   - Precision degrees of freedom
%
% Each input parameter can be missing or empty, in which case:
%   isempty(MU) ->  isempty(b)  -> no prior over MU (ML optimisation)
%               -> ~isempty(b)  -> 'Start' method
%   isempty(b)  -> ~isempty(MU) -> eps (most uninformative prior)
%   isempty(V)  ->  isempty(n)  -> no prior over V (ML optimisation)
%               -> ~isempty(n)  -> eye(P)
%   isempty(n)  -> ~isempty(V)  -> P (most uninformative prior)

MU = [];
b  = [];
V  = [];
n  = [];
if isempty(pr0)
    return
end

% -------------------------------------------------------------------------
% Mean
if numel(pr0) > 0
    MU = pr0{1};
end
% -------------------------------------------------------------------------
% Mean df
if numel(pr0) > 1
    b = pr0{2};
end
% -------------------------------------------------------------------------
% Mean df [default]
if ~isempty(MU) && isempty(b)
    b = eps;
end
if ~isempty(b)
    b = b(:)';
    if numel(b) < K
        b = padarray(b, [0 K - numel(b)], 'replicate', 'post');
    end
end
% -------------------------------------------------------------------------
% Scale
if numel(pr0) > 2
    V = pr0{3};
end
% -------------------------------------------------------------------------
% Scale df
if numel(pr0) > 3
    n = pr0{4};
end
% -------------------------------------------------------------------------
% Scale df [default]
if ~isempty(V) && isempty(n)
    n = P;
end
if ~isempty(n)
    n = n(:)';
    if numel(n) < K
        n = padarray(n, [0 K - numel(n)], 'replicate', 'post');
    end
end
% -------------------------------------------------------------------------
% Scale [default]
if ~isempty(V)
    if size(V,3) < K
        V = padarray(V, [0 0 K - size(V,3) ], 'replicate', 'post');
    end
    if size(V,1) == 1
        V = bsxfun(@times, V, eye(P));
    end
end



