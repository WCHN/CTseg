function [B,coeff,lb,ok] = spm_bias_loop(X, Z, cluster, bases, varargin)
%__________________________________________________________________________
%
% Fit a multiplicative bias field to observed data by optimising a Gaussian
% mixture model.
%
% This function is the core of the fitting process. However, it needs all
% inputs to be well formatted and initialised and is, thus, not the usual
% entry point. To fit a bias field without having to bother with these 
% issues, use spm_bias instead.
%
% FORMAT [bias,coeff,lb] = spm_bias_loop(obs,resp,cluster,bases,...)
%
% MANDATORY
% ---------
%
% obs      - NxP observations
% resp     - NxK responsibilities
% cluster <- {MU,A}, {{MU,b},A}, {MU,{V,n}}, or {{MU,b},{V,n}}
%   MU - PxK   means
%   b  - 1xK   mean d.f. [0=ML]
%   A  - PxPxK precision matrices
%   V  - PxPxK scale matrices
%   n  - 1xK   precision d.f. [0=ML]
% bases   <- {Bx, By, ...}, where Bd are NdxJd basis functions
%       
% KEYWORD
% -------
% Coefficients   - Pre-computed bias coefficients [0]
% BiasField      - Pre-computed bias field [1]
% LowerBound     - Pre-computed lower bound structure with fields:
%                   sum, last, X, B
% RegPrecision   - Precision matrix for the bias coefficients [ML]
% RegParam       - Regularisation parameter for the bias coefficients [1]
% MissingCode    - C, {C}, or {C,L} [recompute]
%   C - Nx1 Image of missing code
%   L - List of unique codes
% IterMax        - Max number of GN iterations [1024]
% Tolerance      - Gain tolerance to stop the GN algorithm [1e-4]
% LineSearch     - Max number of line search iterations [6]
% BinWidth       - 1xP Bin width [0]
% JointOptim     - Optimisation strategy [true=Joint]/false=Iterative
% Verbose        - Verbosity level: [0]= quiet
%                                    1 = write (lower bound)
%                                    2 = plot (lower bound)
%                                    3 = plot more (gmm fit)
% 
% OUTPUT
% ------
% 
% bias  - NxP         Bias field
% coeff - prod(Jd)xP  Bias coefficients
% lb    - Structure with fields: sum, last, X, B
% ok    - True if at least one improved value was found
% 
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

lb0 = struct('sum', NaN, 'X', [], 'B', [], 'XB', []);

% -------------------------------------------------------------------------
% Parse inputs
q = inputParser;
q.FunctionName = 'spm_bias_loop';
q.addParameter('LowerBound',     lb0,   @isstruct);
q.addParameter('Coefficients',   [],    @isnumeric);
q.addParameter('BiasField',      [],    @isnumeric);
q.addParameter('RegPrecision',   [],    @isnumeric);
q.addParameter('RegParam',       1,     @(X) isscalar(X) && isnumeric(X));
q.addParameter('MissingCode',    {},    @(X) isnumeric(X) || iscell(X));
q.addParameter('IterMax',        1024,  @(X) isscalar(X) && isnumeric(X));
q.addParameter('Tolerance',      1e-4,  @(X) isscalar(X) && isnumeric(X));
q.addParameter('LineSearch',     6,     @(X) isscalar(X) && isnumeric(X));
q.addParameter('BinWidth',       0,     @isnumeric);
q.addParameter('JointOptim',     [],    @(X) isscalar(X) && islogical(X));
q.addParameter('Verbose',        0,     @(X) isnumeric(X) || islogical(X));
q.parse(varargin{:});
lb    = q.Results.LowerBound;
B     = q.Results.BiasField;
coeff = q.Results.Coefficients;
ICO   = q.Results.RegPrecision;
prms  = q.Results.RegParam;
BW    = q.Results.BinWidth;
C     = q.Results.MissingCode;
IterMax    = q.Results.IterMax;
Tolerance  = q.Results.Tolerance;
LineSearch = q.Results.LineSearch;
JointOptim = q.Results.JointOptim;
Verbose    = q.Results.Verbose;

% -------------------------------------------------------------------------
% Unfold inputs
b   = 0;  % Mean degrees of freedom (posterior)
n   = 0;  % Precision degrees of freedom (posterior)
V   = []; % Scale matrix (posterior)
L   = []; % List of unique codes

if iscell(C)
    codes = C;
    if numel(codes) >= 1
        C = codes{1};
    else
        C = [];
    end
    if numel(codes) >= 2
        L = codes{2};
    else
        L = unique(C);
    end
    clear codes
end
if ~iscell(cluster) || numel(cluster) < 2
    error('[spm_bias_loop] At least one mean and one precision matrix are needed.');
else
    if ~iscell(cluster{1})
        MU = cluster{1};
    else
        MU = cluster{1}{1};
        if numel(cluster{1}) >= 2
            b = cluster{1}{2};
        end
    end
    if ~iscell(cluster{2})
        A = cluster{2};
    else
        A = cluster{2}{1};
        if numel(cluster{2}) >= 2
            n = cluster{2}{2};
            if sum(n) > 0
                V = A;
                A = bsxfun(@times, V, reshape(n, 1, 1, []));
            end
        end
    end
end
if sum(b) > 0
    mean    = {MU,b};
else
    mean = {MU};
end
if sum(n) > 0
    prec    = {V,n};
else
    prec = {A};
end
cluster = {MU, prec{:}};

% -------------------------------------------------------------------------
% Get lattice
lattice = zeros(1,numel(bases));
nbcomp  = zeros(1,numel(bases));
for i=1:numel(bases)
    lattice(i) = size(bases{i},1);
    nbcomp(i)  = size(bases{i},2);
end
P = size(MU,1);
K = size(MU,2);
J = prod(nbcomp);


if isempty(C) && any(any(isnan(reshape(X, [], P))))
    C = spm_gmm_lib('obs2code', reshape(X, [], P));
    L = unique(C);
end

% -------------------------------------------------------------------------
% Reshape everything in matrix form
X = reshape(X, [], P);
Z = reshape(Z, [], K);


% -------------------------------------------------------------------------
% Initialise coefficients if needed
if isempty(coeff)
    coeff = zeros([nbcomp P], 'like', B);
end

% -------------------------------------------------------------------------
% Initialise lower bound + reconstruct bias field
if isnan(lb.sum(end))
    B   = spm_bias_lib('reconstruct', bases, coeff, 'add');
    B   = reshape(B, [], P);
    B(~isfinite(X)) = 0;
    LXB = sum(B, 1)';
    B   = exp(B);
    LX  = spm_bias_lib('objective', X, Z, B, mean, prec, {C, L}, BW);
    LB  = zeros(P,1);
    for p=1:P
        c1 = reshape(coeff(:,:,:,p), [], 1);
        LB(p) = -0.5 * get_param(prms, 1) * (c1' * ICO * c1);
    end
    llb = sum([sum(LX,'omitnan') sum(LXB,'omitnan') sum(LB,'omitnan')], 'omitnan');
    lb.X(:,end+1)  = LX;
    lb.XB(:,end+1) = LXB;
    lb.B(:,end+1)  = LB;
    lb.sum(end+1)  = llb;
else
    B   = spm_bias_lib('reconstruct', bases, coeff, 'mult');
    B   = reshape(B, [], P);
    if ~isempty(lb.X),  LX  = lb.X(:,end);  else, LX  = 0; end
    if ~isempty(lb.XB), LXB = lb.XB(:,end); else, LXB = zeros(P,1); end
    if ~isempty(lb.B),  LB  = lb.B(:,end);  else, LB  = zeros(P,1); end
    llb = sum([sum(LX,'omitnan') sum(LXB,'omitnan') sum(LB,'omitnan')], 'omitnan');
end

% -------------------------------------------------------------------------
% Select optimisation strategy
if isempty(JointOptim)
    JointOptim = numel(coeff) <= 3000;
end
if JointOptim
    R = 1;  % < Iterations needed (Joint: 1, Iter: P)
    Q = P;  % < Size of Grad/Hess (Joint: P, Iter: 1)      
else
    R = P;
    Q = 1;
end
gain = zeros(1,P); % < Allocate array to store gain

% -------------------------------------------------------------------------
% Gauss-Newton loop
oneok = false; % < At least one success?
for gnit=1:IterMax

    % ---------------------------------------------------------------------
    % Choose the right regularisation
    prm = get_param(prms, gnit);
    
    for r=1:R % < Joint: one loop, Iter: P loops
    
        if JointOptim
            list_p = 1:P;
            whichd = [];
        else
            list_p = r;
            whichd = r;
        end
        
        % -----------------------------------------------------------------
        % Compute grad/hess
        [g,H] = spm_bias_lib('derivatives', whichd, B.*X, bases, Z, cluster, {C,L}, (bsxfun(@times, B, BW).^2)/12);
        g = reshape(g, J, Q);
        H = reshape(H, J, Q, J, Q);
        for q=1:Q
            p          = list_p(q);
            c1         = reshape(coeff(:,:,:,p), [], 1);
            g(:,q)     = g(:,q)     + prm * ICO * c1;
            H(:,q,:,q) = H(:,q,:,q) + prm * reshape(ICO, [J 1 J]);
        end

        % -----------------------------------------------------------------
        % Gauss-Newton
        g = g(:);
        H = reshape(H, J*Q, J*Q);
        H = (H+H')/2;
        % H = spm_matcomp('LoadDiag', H);
        dc = H\g;
        dc = reshape(dc, [], Q);

        % -----------------------------------------------------------------
        % Line search
        armijo = 1;
        coeff0 = reshape(coeff, [], P);
        llb0   = llb;
        ok     = false;
        for ls=1:LineSearch
            % -------------------------------------------------------------
            % Update bias field
            coeff           = reshape(coeff,  [], P);
            coeff(:,list_p) = coeff0(:,list_p) - armijo * dc;
            coeff           = reshape(coeff,  [nbcomp P]);
            for q=1:Q
                p      = list_p(q);
                c1     = coeff(:,:,:,p);
                B(:,p) = reshape(spm_bias_lib('reconstruct', bases, c1, 'add'), [], 1);
                B(~isfinite(X(:,p)),p) = 0;
            end

            % -------------------------------------------------------------
            % Update objective function
            LXB(list_p,1) = sum(B(:,list_p),1)'; % < Update LogDet(B)
            B(:,list_p) = exp(B(:,list_p));    % < Exponentiate
            LX  = spm_bias_lib('objective', X, Z, B, mean, prec, {C, L}, BW);
            for q=1:Q
                p     = list_p(q);
                c1    = reshape(coeff(:,:,:,p), [], 1);
                LB(p) = -0.5 * prm * (c1' * ICO * c1);
            end
            llb = sum([sum(LX,'omitnan') ...
                       sum(LXB,'omitnan') ...
                       sum(LB,'omitnan')], 'omitnan');

            % -------------------------------------------------------------
            % Check improvement
            if llb > llb0
                ok = true;
                break
            else
                armijo = armijo/2;
            end
        end
        % -----------------------------------------------------------------
        % If no improvement, use previous value
        if ~ok
            ls    = 0;
            llb   = llb0;
            coeff = reshape(coeff0,  [nbcomp P]);
            for q=1:Q
                p      = list_p(q);
                c1     = coeff(:,:,:,p);
                B(:,p) = reshape(spm_bias_lib('reconstruct', bases, c1, 'mult'), [], 1);
                B(~isfinite(X(:,p)),p) = 1;
            end
            lb.X(:,end+1)  = lb.X(:,end);
            lb.XB(:,end+1) = lb.XB(:,end);
            lb.B(:,end+1)  = lb.B(:,end);
        else
            lb.X(:,end+1)  = LX;
            lb.XB(:,end+1) = LXB;
            lb.B(:,end+1)  = LB;
        end
        oneok = oneok || ok;

        
        % -----------------------------------------------------------------
        % Plot Bias
        if Verbose(1) >= 3
            spm_bias_lib('Plot', 'Bias', X, B, Z, lattice);
        end
        
        % -----------------------------------------------------------------
        % Check convergence
        [lb,gain(r)] = check_convergence(lb, gnit, ls, list_p, Verbose(1));
    end
    if sum(gain) < Tolerance
        break;
    end
    
end
ok = oneok; % At least one success

% =========================================================================
function prm = get_param(prms, i)

if numel(prms) > i
    prm = prms(i);
else
    prm = prms(end);
end

% =========================================================================
function [lb,gain] = check_convergence(lb, em, ls, whichp, verbose)
% FORMAT [lb,gain] = check_convergence(lb, em, verbose)
% lb      - Lower bound structure with fields X, Z, P, MU, A, sum, last
% em      - EM iteration
% verbose - Verbosity level (>= 0)
%
% Compute lower bound (by summing its parts) and its gain
% + print info

fields = fieldnames(lb);
lb.sum(end+1) = 0;
for i=1:numel(fields)
    field = fields{i};
    if ~any(strcmpi(field, {'sum' 'last'})) && ~isempty(lb.(field)) && ~isnan(lb.(field)(end))
        lb.sum(end) = lb.sum(end) + sum(lb.(field)(:,end));
    end
end
gain = (lb.sum(end) - lb.sum(end-1))/(max(lb.sum(:), [], 'omitnan')-min(lb.sum(:), [], 'omitnan'));
if verbose >= 1
    if verbose >= 2
        spm_gmm_lib('plot', 'LB', lb)
    end
    switch sign(gain)
        case  1,    incr = '(+)';
        case -1,    incr = '(-)';
        case  0,    incr = '(=)';
        otherwise,  incr = '';
    end
    if ls == 0
        lsres = ':(';
    else
        lsres = sprintf(':D (%d)', ls);
    end
    if numel(whichp) == 1
        name = sprintf('bias%d', whichp);
    else
        name = 'bias';
    end
    fprintf('%-5s | %4d | lb = %-10.6g | gain = %-10.4g | %3s | %-7s\n', name, em, lb.sum(end), gain, incr, lsres);
end
gain = abs(gain);