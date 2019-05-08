function varargout = spm_bias(X, varargin)
%__________________________________________________________________________
%
% Estimate a spatial bias field by fitting a GMM to the data.
%
% FORMAT [B,X,...] = spm_bias(X,...)
% 
% MANDATORY
% ---------
% X - [LAT]xP Multichannel image or volume
%
% OPTIONAL
% ---------
% K - Number of GMM classes [10]
%
% KEYWORD
% -------
% VoxelSize  - 1x[D] Vector of voxel sizes [1] (mm/vox)
% FWHM       - 1x[D] Vector of FWHM to choose the number of bases [60] (mm)
% NbBases    - 1x[D] Number of bases along each dimension [0=use FWHM]
% RegParam   - 1x3   Regularisation parameters [0 0 10] (abs/mem/ben)
% IterMax    - Maximum number of EM iterations [1000]
% Tolerance  - Convergence criterion (lower bound gain) [1e-4]
% BinWidth   - 1x[P] Bin width (histogram mode: add bits of variance) [0]
% InputDim   - Input space dimension [0=try to guess]
% Verbose    - Verbosity level: [0]= quiet
%                                1 = write (lower bound)
%                                2 = plot (lower bound)
%                                3 = plot more (gmm fit)
% GMM        - Cell of GMM options:
%              GaussPrior     - {MU0,b0,V0,n0} [{}=ML]
%              PropPrior      - a0 [0=ML]
%              IterMax        - Max number of EM iterations [1024]
%              Tolerance      - Gain tolerance to stop the EM algorithm [1e-4]
%              SubIterMax     - Max number of sub-EM iterations [1024]
%              SubTolerance   - Sub-EM gain tolerance [1e-4]
% Bias       - Cell of Bias options:
%              IterMax        - Max number of GN iterations [1024]
%              Tolerance      - Gain tolerance to stop the GN algorithm [1e-4]
%              LineSearch     - Max number of line search iterations [6]
%              JointOptim     - Optimisation strategy [true=Joint]/false=Iterative
%
% OUTPUT
% ------
% B    - [LAT]xP   bias field
% X    - [LAT]xP   corrected input volume
% ...  - Same output as spm_gmm, in the same order > help spm_gmm
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Parse inputs
p = inputParser;
p.FunctionName = 'spm_bias';
p.addRequired('X',                      @isnumeric);
p.addOptional('K',           10,        @(X) isscalar(X) && isnumeric(X));
p.addParameter('VoxelSize',  1,         @isnumeric);
p.addParameter('FWHM',       60,        @isnumeric);
p.addParameter('NbBases',    0,         @isnumeric);
p.addParameter('RegParam',   1E7,       @isnumeric);
p.addParameter('GMM',        {},        @iscell);
p.addParameter('Bias',       {},        @iscell);
p.addParameter('IterMax',    1000,      @(X) isscalar(X) && isnumeric(X));
p.addParameter('Tolerance',  1e-4,      @(X) isscalar(X) && isnumeric(X));
p.addParameter('BinWidth',   0,         @isnumeric);
p.addParameter('InputDim',   0,         @(X) isscalar(X) && isnumeric(X));
p.addParameter('Verbose',    0,         @(X) isnumeric(X) || islogical(X));
p.parse(X, varargin{:});
K          = p.Results.K;
vs         = p.Results.VoxelSize;
fwhm       = p.Results.FWHM;
nbbases    = p.Results.NbBases;
RegParam   = p.Results.RegParam;
GMMopt     = p.Results.GMM;
Biasopt    = p.Results.Bias;
P          = p.Results.InputDim;
BinWidth   = p.Results.BinWidth;
IterMax    = p.Results.IterMax;
Verbose    = p.Results.Verbose;
Tolerance  = p.Results.Tolerance;

% -------------------------------------------------------------------------
% Prepare input
[P,latX] = dimFromObservations(P, X);
dim      = numel(latX);
C        = spm_gmm_lib('obs2code', reshape(X, [], P));     % Code image
L        = unique(C);                                      % List of codes

% -------------------------------------------------------------------------
% Prepare bases
if sum(nbbases) == 0
    nbbases = spm_bias_lib('fwhm2nbcomp', latX, vs, fwhm);
end
[Bx,By,Bz] = spm_bias_lib('dcbasis', latX, nbbases);
switch dim
    case 2, bases = {Bx,By};
    case 3, bases = {Bx,By,Bz};
end
ICO   = spm_bias_lib('regulariser', 'bending', latX, nbbases, vs);
coeff = zeros([nbbases P], 'like', X);
field = spm_bias_lib('reconstruct', bases, coeff, 'mult');

% -------------------------------------------------------------------------
% GMM prior
MU0 = zeros(P,K);
BX  = reshape(field.*X, [], P);
minval = min(BX, [], 'omitnan');
maxval = max(BX, [], 'omitnan');
clear BX
for p=1:P
    tmp = linspace(minval(p), maxval(p), K+1);
    MU0(p,:) = (tmp(1:end-1) + tmp(2:end))/2;
end
A0 = diag(((maxval - minval)./(2*K)).^(-2));
A0 = repmat(A0, [1 1 K]);
n0 = 10 * ones(1,K);
b0 = 10 * ones(1,K);
V0 = bsxfun(@rdivide, A0, reshape(n0, [1 1 K]));
a0 = 0 * ones(1,K);

mean    = {MU0, b0};
prec    = {V0, n0};
cluster = {mean prec};
a       = a0;
PI      = repmat(1/K, [1 K]);
Z       = [];

% -------------------------------------------------------------------------
% EM loop
lb  = struct('sum', NaN, 'X', [], 'XB', [], 'B', [], ...
                         'Z', [], 'P', [], 'MU', [], 'A', []);
ThisRegParam = RegParam(1);
RegParam     = RegParam(2:end);
for em=1:IterMax
    
    
    % ---------------------------------------------------------------------
    % Save previous lower bound value
    obj0 = lb.sum(end);
    
    % ---------------------------------------------------------------------
    % Optimise GMM
    [Z,cluster,prop,lb] = spm_gmm_loop(reshape(field, [], P).*reshape(X, [], P), ...
        cluster, {'Prop', PI, 'Dir', a}, ...
        'PropPrior',      a0, ...
        'LowerBound',     lb, ...
        'Missing',        true, ...
        'MissingCode',    {C,L}, ...
        'BinUncertainty', (bsxfun(@times, reshape(field, [], P), reshape(BinWidth, 1, [])).^2)/12, ...
        'Verbose',        Verbose, ...
        GMMopt{:});
    MU = cluster.MU;
    b  = cluster.b;
    A  = cluster.A;
    V  = cluster.V;
    n  = cluster.n;
    PI = prop.Prop;
    a  = prop.Dir;
    if sum(b) > 0, mean = {MU b};
    else,          mean = {MU};
    end
    if sum(n) > 0, prec = {V n};
    else,          prec = {A};
    end
    cluster = {mean prec};
    
    
    % ---------------------------------------------------------------------
    % Plot Bias
    if Verbose(1) >= 3
        spm_bias_lib('Plot', 'Bias', X, field, Z, latX);
    end
    
    % ---------------------------------------------------------------------
    % Optimise Bias Field
    [field,coeff,lb,ok] = spm_bias_loop(X, Z, cluster, bases, ...
        'Coefficients',   coeff, ...
        'BiasField',      field, ...
        'LowerBound',     lb, ...
        'RegPrecision',   ICO, ...
        'RegParam',       ThisRegParam, ...
        'MissingCode',    {C,L}, ...
        'BinWidth',       BinWidth, ...
        'Verbose',        Verbose, ...
        Biasopt{:});
    if ~ok
        break;
    end
    
    % ---------------------------------------------------------------------
    % Check convergence
    obj  = lb.sum(end);
    den  = max(lb.sum(:), [], 'omitnan')-min(lb.sum(:), [], 'omitnan');
    gain = check_convergence(obj, obj0, den, em, Verbose(1));
    if gain < Tolerance
        if ~isempty(RegParam)
            PrevRegParam = ThisRegParam;
            ThisRegParam = RegParam(1);
            RegParam     = RegParam(2:end);
            lb.B(:,end+1)  = lb.B(:,end) * ThisRegParam / PrevRegParam;
        else
            break
        end
    end
end

% -------------------------------------------------------------------------
% Prepare output
Z             = reshape(Z, [latX K]);
X             = reshape(X, [latX P]);
field         = reshape(field, [latX P]);
varargout{1}  = field;
varargout{2}  = field.*X;
varargout{3}  = Z;
varargout{4}  = MU;
varargout{5}  = A;
varargout{6}  = PI;
varargout{7}  = b;
varargout{8}  = V;
varargout{9}  = n;
varargout{10} = a;

% =========================================================================
function gain = check_convergence(obj, obj0, den, it, verbose)
% FORMAT gain = check_convergence(obj, obj0, den, it, verbose)

gain = (obj - obj0)/den;
if verbose >= 1
    switch sign(gain)
        case  1,   incr = '(+)';
        case -1,   incr = '(-)';
        case  0,   incr = '(=)';
        otherwise, incr = '';
    end
    fprintf('%-5s | %4d | lb = %-10.6g | gain = %-10.4g | %3s\n', 'b+g', it, obj, gain, incr);
end

% =========================================================================
function [P,latX] = dimFromObservations(P, X)
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
        P = dimX(end);
    end
end
if P == 1
    latX = dimX;
else
    latX = dimX(1:end-1);
end
