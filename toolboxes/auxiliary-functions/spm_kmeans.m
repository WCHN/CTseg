function [L,C,SUMD,D] = spm_kmeans(X, K, varargin)
% _________________________________________________________________________
%
% Learn a K-means clustering from observed [weighted] data.
%
% FORMAT [L,C,SUMD,D] = spm_kmeans(X,K,...)
% 
% MANDATORY
% ---------
% X - NxP matrix of observed values
% K - Number of cluster
% 
% OPTIONAL
% --------
% W - Nx1 Vector of weights associated with each observation [1]
%
% KEYWORD
% -------
% Distance   - Distance between points: ['sqeuclidian'], 'cityblock'
% Start      - Starting method: ['plus'], 'sample', 'uniform'
%                   or a KxPxR array with R the number of replicates
% Replicates - Number of replicates with different starts [1]
% Missing    - Keep rows with missing data [true]
% Order      - Centroid ordering method: ['magnitude'], 'total', 'random',
%              'intensity'
% Verbose    - Verbosity level [0]
%
% OUTPUT
% ------
% L    - Nx1 labelling
% C    - KxP centroids
% SUMD - Kx1 sum of distances inside each cluster
% D    - NxK distances to each centroid
% _________________________________________________________________________
%
% Use learned clusters to segment an image.
%
% FORMAT [L,D] = spm_kmeans(X,C,...)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Parse inputs
p = inputParser;
p.FunctionName = 'spm_kmeans';
p.addRequired('X', @isnumeric);
p.addRequired('K', @isnumeric);
p.addOptional('W',           1,             @isnumeric);
p.addParameter('Distance',   'sqeuclidian', @ischar);
p.addParameter('Start',      'plus',        @(X) ischar(X) || isnumeric(X));
p.addParameter('Replicates', 1,             @isnumeric);
p.addParameter('Order',      'magnitude',   @ischar);
p.addParameter('Missing',    true,          @isscalar);
p.addParameter('Verbose',    0,             @isscalar);
p.parse(X, K, varargin{:});
W = p.Results.W;
X = single(X);

Replicates = p.Results.Replicates;

% -------------------------------------------------------------------------
% Special case: Apply model
% > Here, we use learned centroids to segment an image
if numel(K) > 1
    [L,C] = kmeans_apply(X, K, p.Results.Distance, p.Results.Missing);
    return
end

% -------------------------------------------------------------------------
% Guess cluster/replicates from provided centroids
if ~ischar(p.Results.Start)
    Replicates = size(p.Results.Start, 3);
    K          = size(p.Results.Start, 1);
end

% -------------------------------------------------------------------------
% Vector case
row_vector = size(X,1) == 1 && numel(size(X)) == 2;
if row_vector
    X = X';
end
dim = size(X);
X   = reshape(X, [], dim(end));

% -------------------------------------------------------------------------
% Prepare weights
if numel(W) == 1
    W = W * ones([size(X,1) 1],'single');
end

if p.Results.Missing % Deal with missing data
    WW = repmat(W, [1 size(X,2)]);
    WW(isnan(X)) = 0;
else % Discard rows with missing values
    N0      = size(X,1);
    missing = any(isnan(X),2);
    W       = W(~missing);
    WW      = W;
    X       = X(~missing,:);
end

% -------------------------------------------------------------------------
% Convert distance function to handle
dist_fun = str2func(p.Results.Distance);

E00 = inf;
r0  = 0;
% -------------------------------------------------------------------------
% For each replicate
for r=1:Replicates
    
    % ---------------------------------------------------------------------
    % Initial centroids
    C1 = start(p.Results.Start, X, WW, K, dist_fun, r);
    
    % ---------------------------------------------------------------------
    % Iterate
    E0 = inf;
    for i=1:1000
        % -----------------------------------------------------------------
        % Distance to previous centroids
        D1 = dist_fun(X,C1);
        % -----------------------------------------------------------------
        % Clustering to nearest centroid
        [MinD, L1] = min(D1, [], 2);
        L1         = single2int(L1);
        % -----------------------------------------------------------------
        % Check convergence
        E = sum(MinD.*W);
        if (E0-E)/E0 < 1e-7
            break
        end
        E0 = E;
        % -----------------------------------------------------------------
        % Update centroids
        % > C = sum(W*X)/sum(W)
        for k=1:K
            C1(k,:) = bsxfun(@rdivide,sum(bsxfun(@times,X(L1==k,:),WW(L1==k,:)),1,'omitnan'),sum(WW(L1==k,:),1));
        end
    end
    if p.Results.Verbose
        fprintf('r = %2d | imax = %3d | E = %3g\n', r, i, E);
    end
    
    % ---------------------------------------------------------------------
    % Keep best replicate
    if E < E00
        E00  = E;
        C    = C1;    clear C1
        L    = L1;    clear L1
        D    = D1;    clear D1
        SUMD = zeros(K,1);
        for k=1:K
            SUMD(k) = sum(MinD(L==k) .* W(L==k));
        end
        r0 = r;
    end
    
end
if p.Results.Verbose && Replicates > 1
    fprintf('Best | r = %2d | E = %3g\n', r0, E00);
end

if ~isempty(p.Results.Order)
    % -------------------------------------------------------------------------
    % Order centroids
    [L,C,SUMD,D] = order(p.Results.Order,L,C,SUMD,D,W);
end

% -------------------------------------------------------------------------
% Replace discarded missing values
if ~p.Results.Missing
    L1            = L;
    L             = zeros([N0 1], 'single');
    L(~missing)   = L1;
    clear L1;
    D1            = D;
    D             = NaN([N0 size(D1, 2)], 'single');
    D(~missing,:) = D1;
    clear D1;
end

% -------------------------------------------------------------------------
% Reshape output
L = reshape(L, [dim(1:end-1) 1]);
D = reshape(D, [dim(1:end-1) size(C,1)]);
if row_vector
    L = L';
end

% =========================================================================
function [L,D] = kmeans_apply(X, C, dist, missing)
% FORMAT [L,D] = kmeans_apply(X, C)
% Classify observations based on known centroids.

dist = str2func(dist);

% Reshape input
dim = size(X);
if dim(end) == size(C,2)
    dim = dim(1:end-1);
end
X = reshape(X, [], size(C,2));

% Compute distance
D = dist(X,C);
% Get closest centroid
[~, L] = min(D, [], 2);
L      = single2int(L);
% Replace missing data
if ~missing
    msk      = any(isnan(X),2);
    L(msk)   = 0;
    D(msk,:) = NaN; 
end
% Reshape output
L = reshape(L, [dim 1]);
D = reshape(D, [dim size(C,1)]);

% =========================================================================
function [L,C,SUMD,D] = order(method,L,C,SUMD,D,W)
% FORMAT [L,C,SUMD,D] = order(method,L,C,SUMD,D,W)
% 
% Order centroids (and centroid related data) w.r.t. a given measure.

switch method
    case 'total'
        measure = zeros([size(C,1) 1],'single');
        for k=1:size(C,1)
            measure(k) = sum(W(L == k));
        end
    case 'magnitude'
        measure = sqrt(sum(C.^2, 2));
    case 'intensity'
        measure = C(:,1);        
    otherwise
        return
end

[~,I] = sort(measure);
C     = C(I,:);
SUMD  = SUMD(I);
D     = D(:,I);
oldL  = L;
L     = zeros(size(oldL), 'single');
for k=1:size(C,1)
    L(oldL == I(k)) = k;
end

% =========================================================================
function [C,K] = start(method, X, WW, K, dist, r)
% FORMAT C = start(method, X, K, dist)
% method - Method to use to select starting centroids
%               'plus', 'sample', 'uniform' or provided matrix
% X      - Vector of NxP observed values
% K      - Number of clusters
% dist   - Distance function handle
% r      - Replicate index
%
% Compute starting centroids

if nargin < 5
    r = 1;
end

if isnumeric(method)
    % Provided
    K = size(method,1);
    if size(method, 3) == 1
        C = method;
    else
        C = method(:,:,r);
    end

    return
end

switch method
    
    case 'plus'
    % K-means ++
    % The first centroids is selected at random from the observed values
    % Subsequent centroids are selected at random with probability
    % proportional to theirs distance to the closest centroid.
    % This allows to start with centroids that are reasonably far from one
    % another.
        X = X(~any(WW==0,2),:); % Remove all rows w/ NaNs or w/o obs
        C = zeros([K size(X,2)],'single');
        i = randi(size(X,1));
        C(1,:) = X(i,:);
        for k=2:K
            P = dist(X, C(1:k-1,:));
            P = min(P, [], 2);
            i = randvalue(P);
            C(k,:) = X(i,:);
        end
        
    case 'sample'
    % Sample uniform
    % All centroids are selected at random from the observed values.
    % They are all unique (to avoid starting with several identical
    % centroids)
        X = X(~any(WW==0,2),:); % Remove all rows w/ NaNs or w/o obs
        i = randperm(size(X,1));
        i = i(1:K);
        C = X(i,:);
        
    case 'uniform'
    % Range uniform
    % All centroids are selected at random from the continuous range of 
    % observed values.
    % They are all unique (to avoid starting with several identical
    % centroids)
        minval = min(X, [], 1, 'omitnan');
        maxval = max(X, [], 1, 'omitnan');
        C = rand([K size(X,2)],'single');
        C = bsxfun(@times, C, maxval - minval) + minval;
        
    otherwise
        error('Undefined method!')
end

% =========================================================================
function x = randvalue(P,X)
% FORMAT x = randvalue(P,[X])
% P - A vector of probabilities. If it does not sum to 1, it will be
%     normalised.
% X - Values to sample. By default: indices of P.
%
% Return a random value according to known probabilities.

if nargin < 2
    X = 1:numel(P);
end

p      = cumsum([0 P(:).'/sum(P(:))]);
p(end) = 1e3*eps + p(end);
i      = 0;
while i == 0
    [~, ~, i] = histcounts(rand,p);
end
x      = X(i);

% =========================================================================
function D = sqeuclidian(X,C)
% FORMAT D = sqeuclidian(X,C)
% X - NxP observations
% C - KxP centroids
% D - NxK distance to centroids
%
% Compute L2 distance between each observation and each centroid.

X = reshape(X, [size(X,1) 1 size(X,2)]);
C = reshape(C, [1 size(C,1) size(C,2)]);
D = bsxfun(@minus, X, C);
D(isnan(D)) = 0;
D = sqrt(sum(D.^2, 3));

% =========================================================================
function D = cityblock(X,C)
% FORMAT D = cityblock(X,C)
% X - NxP observations
% C - KxP centroids
% D - NxK distance to centroids
%
% Compute L1 distance between each observation and each centroid.

X = reshape(X, [size(X,1) 1 size(X,2)]);
C = reshape(C, [1 size(C,1) size(C,2)]);
D = bsxfun(@minus, X, C);
D(isnan(D)) = 0;
D = sum(abs(D), 3);

% =========================================================================
function L = single2int(L)
% FORMAT L = single2int(L)
%
% Find the best suited integer type to convert L, based on min and max
% values

minval   = min(L(:));
maxval   = max(L(:));
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
    type = 'single';
end
func = str2func(type);
L = func(L);