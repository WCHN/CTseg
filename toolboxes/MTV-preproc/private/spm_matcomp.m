function varargout = spm_matcomp(varargin)
%__________________________________________________________________________
% Collection of tools for matrices and vector/tensor fields.
%
% FORMAT out = spm_matcomp('name', input)
%
% FORMAT      A = spm_matcomp('LoadDiag', A)
% FORMAT     ld = spm_matcomp('LogDet', A)
% FORMAT      s = spm_matcomp('logsumexp',b,dm)
% FORMAT     iA = spm_matcomp('Inv', A)
% FORMAT      c = spm_matcomp('Pointwise3', a, (b), (op))
% FORMAT      c = spm_matcomp('Pointwise', a, (b), (op), (sym))
% FORMAT [Q,ll] = spm_matcomp('softmax',Q,scl)
% FORMAT    ind = spm_matcomp('SymIndices', k)
% FORMAT    ind = spm_matcomp('SymIndices', k)
% FORMAT      M = spm_matcomp('threshold_eig',M)
%
% FORMAT help spm_matcomp>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging
    if nargin == 0
        help spm_matcomp
        error('Not enough argument. Type ''help spm_matcomp'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case 'pointwise3'
            [varargout{1:nargout}] = pointwise3(varargin{:});
        case 'pointwise'
            [varargout{1:nargout}] = pointwise(varargin{:});
        case 'logdet'
            [varargout{1:nargout}] = logdet(varargin{:});
        case 'symindices'
            [varargout{1:nargout}] = symIndices(varargin{:});
        case 'loaddiag'
            [varargout{1:nargout}] = loaddiag(varargin{:});
        case {'inv', 'inverse'}
            [varargout{1:nargout}] = inv(varargin{:});
        case 'threshold_eig'
            [varargout{1:nargout}] = threshold_eig(varargin{:});            
        case 'logsumexp'
            [varargout{1:nargout}] = logsumexp(varargin{:});                   
        case 'softmax'
            [varargout{1:nargout}] = softmax(varargin{:});                    
        otherwise
            help spm_matcomp
            error('Unknown function %s. Type ''help spm_matcomp'' for help.', id)
    end
end

%% === logdet =============================================================

function ld = logdet(A)
% FORMAT ld = spm_matcomp('LogDet', A)
% A  - A postive-definite square matrix
% ld - Logarithm of determinant of A
%
% Log-determinant of a positive-definite matrix.
% Cholesky factorisation is used to compute a more stable log-determinant.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

    % Cholseki decomposition of A (A = C' * C, with C upper-triangular)
    [C, p] = chol(numeric(A));
    
    if p > 0
       % A should usually be positive definite, but check anyway.
       warning(['Attempting to compute log determinant of matrix ' ...
                'that is not positive definite (p=%d).'], p);
    end

    % Because C is triangular, |C| = prod(diag(C))
    % Hence: log|C| = sum(log(diag(C)))
    % And:   log|A| = log|C'*C| = log(|C|^2) = 2 * sum(log(diag(C)))
    ld = 2 * sum(log(diag(C)));

end

%% === loaddiag ===========================================================

function A = loaddiag(A)
% FORMAT A = spm_matcomp('LoadDiag', A)
% A  - A square matrix
%
% Load A's diagonal until it is well conditioned for inversion.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

    factor = 1e-7;
    while rcond(A) < 1e-5
        A = A + factor * max([diag(A); eps]) * eye(size(A));
        factor = 10 * factor;
    end

end

%% === inv ================================================================

function A = inv(A)
% FORMAT iA = spm_matcomp('Inv', A)
% A  - A positive-definite square matrix
% iA - Its inverse
%
% Stable inverse of a positive-definite matrix.
% Eigendecomposition is used to compute a more stable inverse.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

    [V,D] = eig(A);
    if any(diag(D) <= 0)
        warning('[spm_matcomp::inv] Matrix has negative eigenvalues')
        D(D <= 0) = eps; % Threshold negative eigenvalues
    end
    D     = loaddiag(D);
    A     = real(V * (D \ V'));
%     A     = 0.5*(A + A.'); % Ensure symetric inverse

end

%% === symindices =========================================================

function [ind, n] = symIndices(k, side)
% FORMAT [ind, n] = spm_matcomp('SymIndices', k, ('k'))
% FORMAT [ind, k] = spm_matcomp('SymIndices', n, 'n')
% k - Length of the linearized and sparse symmetric matrix
% n - Side size of the corresponding square matrix
%
% Returns a converter between row/col and linear indices for sparse
% symmetric matrices.
%
% 1) If needed, finds n the size of the corresponding square matrix
% 2) Returns a n*n matrix which links subscripts (ex: (1,2) or (2,3)) to a
% corresponding linear index when the data is stored "sparsely" to save
% redunduncies.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

    if nargin == 2
        side = strcmpi(side, 'n');
    else
        side = false;
    end

    % 1) Compute n
    if ~side
        n = round(max(roots([1 1 -2*k])));
    else
        n = k;
        k = n*(n+1)/2;
    end
    
    % 2) Build the correspondance matrix
    ind = diag(1:n);
    l = n;
    for i1=1:n
        for i2=(i1+1):n
            l = l + 1;
            ind(i1, i2) = l;
            ind(i2, i1) = l;
        end
    end
    
    if side
        n = k;
    end
end

%% === pointwise3 =========================================================

function c = pointwise3(a, b, op)
% FORMAT c = spm_matcomp('Pointwise3', a, (b), (op))
% a  - [nx nz nz 3 (3)] ND-array of vectors/matrices
% b  - [nx nz nz 3 (3)] ND-array of vectors/matrices
% op - Operation(s) to apply to an element of a before multiplying it 
%      with an element of b:
%      't' (transpose), 'i' (invert), 'd' (determinant).
%      (they can be concatenated: 'itd' means det(inv(a).')
%
% Performs matrix operations at each point of a tensor field.
% If b is provided, at each point:
%    c = op(a) * b
% Else, at each point:
%    c = op(a)
%
% Note that
% - if one (or both) entry is a vector, the adequate transpose
%   operation is usually automatically detected.
% - a and b can be fields of sparse symmetric matrices. In this case, they
%   have the form of a 6 element vector, which is automatically detected.
% - Implicit expansion of the lattice is performed when one of the input
%   dimensions is singular.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

    % --- Default arguments
    if nargin < 3
        if nargin < 2
            op = '';
            binary = false;
        elseif ischar(b)
            op = b;
            binary = false;
        else
            op = '';
            binary = true;
        end
    else
        binary = true;
    end
    
    % --- Check dimensions
    ismat = struct;
    dim   = struct;
    sym   = struct;
    dim.a = size(a);
    sym.a = size(a, 4) > 3;
    ismat.a = sym.a || size(a, 5) > 1;
    if length(dim.a) < 4
        error('Inputs must be vector or tensor fields')
    end
    if binary
        dim.b = size(b);
        sym.b = size(b, 4) > 3;
        ismat.b = sym.b || size(b, 5) > 1;
        if length(dim.b) < 4
            error('Inputs must be vector or tensor fields')
        end
    end
    % The corresponding dimensions of A and B must be equal to each other, 
    % or equal to one.
    if binary
        dima = dim.a(1:3);
        dimb = dim.b(1:3);
        same = dima == dimb;
        diff = min(dima(~same), dimb(~same));
        if  any(diff ~= 1)
            error(['The corresponding lattice dimensions of A and B ' ...
                   'must be equal to each other, or equal to one.']);
        end
    end
    
    % --- Deal with cases
    switch lower(op)
        case ''
            if binary && ismat.a && ismat.b
                if sym.a && sym.b
                    c = pw3_sasb(a, b);
                elseif sym.a
                    c = pw3_samb(a, b);
                elseif sym.b
                    c = pw3_masb(a, b);
                else
                    c = pw3_mamb(a, b);
                end
            elseif binary && ismat.a && ~ismat.b
                if sym.a
                    c = pw3_savb(a, b);
                else
                    c = pw3_mavb(a, b);
                end
            elseif binary && ~ismat.a && ismat.b
                if sym.b
                    c = pw3_vasb(a, b);
                else
                    c = pw3_vamb(a, b);
                end
            elseif  binary && ~ismat.a && ~ismat.b
                c = pw3_vavb(a, b);
            else
                c = a;
            end
        case 't'
            if binary && ismat.a && ismat.b
                if sym.a && sym.b
                    c = pw3_sasb(a, b);
                elseif sym.b
                    c = pw3_masb(a, b);
                elseif sym.a
                    c = pw3_samb(a, b);
                else
                    c = pw3_tmamb(a, b);
                end
            elseif binary && ismat.a && ~ismat.b
                if syma
                    c = pw3_savb(a, b);
                else
                    c = pw3_tmavb(a, b);
                end
            elseif binary && ~ismat.a && ismat.b
                if sym.b
                    c = pw3_vasb(a, b);
                else
                    c = pw3_vamb(a, b);
                end
            elseif  binary && ~ismat.a && ~ismat.b
                c = pw3_vavb(a, b);
            elseif ismat.a
                if sym.a
                    c = a;
                else
                    c = pw3_tma(a);
                end
            else
                c = a;
            end
        case 'i'
            if binary && ismat.a && ismat.b
                if sym.a && sym.b
                    c = pw3_sasb(pw3_isa(a), b);
                elseif sym.a
                    c = pw3_samb(pw3_isa(a), b);
                elseif sym.b
                    c = pw3_masb(pw3_ima(a), b);
                else
                    c = pw3_mamb(pw3_ima(a), b);
                end
            elseif binary && ismat.a && ~ismat.b
                if sym.a
                    c = pw3_savb(pw3_isa(a), b);
                else
                    c = pw3_mavb(pw3_ima(a), b);
                end
            elseif ~binary && ismat.a
                if sym.a
                    c = pw3_isa(a);
                else
                    c = pw3_ima(a);
                end
            else
                error('This case does not make sense')
            end
        case 'd'
            if binary && ismat.a
                if sym.a
                    c = bsxfun(@times, pw3_dsa(a), b);
                else
                    c = bsxfun(@times, pw3_dma(a), b);
                end
            elseif ~binary && ismat.a
                if sym.a
                    c = pw3_dsa(a);
                else
                    c = pw3_dma(a);
                end
            else
                error('This case does not make sense')
            end
        case {'it', 'ti'}
            if binary && ismat.a && ismat.b
                if sym.a && sym.b
                    c = pw3_sasb(pw3_isa(a), b);
                elseif sym.a
                    c = pw3_samb(pw3_isa(a), b);
                elseif sym.b
                    c = pw3_tmasb(pw3_ima(a), b);
                else
                    c = pw3_tmamb(pw3_ima(a), b);
                end
            elseif  binary && ismat.a && ~ismat.b
                if sym.a
                    c = pw3_savb(pw3_isa(a), b);
                else
                    c = pw3_tmavb(pw3_ima(a), b);
                end
            elseif ~binary && ismat.a
                if sym.a
                    c = pw3_isa(a);
                else
                    c = pw3_tma(pw3_ima(a));
                end
            else
                error('This case does not make sense')
            end
        case {'itd', 'tid'}
            if binary && ismat.a
                if sym.a
                    c = bsxfun(@times, b, pw3_dsa(a));
                else
                    c = bsxfun(@ldivide, b, pw3_dma(a));
                end
            elseif ~binary && ismat.a
                if sym.a
                    c = pw3_dsa(a);
                else
                    c = 1 ./ pw3_dma(a);
                end
            else
                error('This case does not make sense')
            end
        otherwise
            error('The provided operation does not make sense')
    end
end

function c = pw3_mamb(a, b)
% a - 3x3 tensor field
% b - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise a * b
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim 3 3], 'like', a);
    
    c(:,:,:,1,1) = bsxfun(@times, a(:,:,:,1,1), b(:,:,:,1,1)) + bsxfun(@times, a(:,:,:,1,2), b(:,:,:,2,1)) + bsxfun(@times, a(:,:,:,1,3), b(:,:,:,3,1));
    c(:,:,:,1,2) = bsxfun(@times, a(:,:,:,1,1), b(:,:,:,1,2)) + bsxfun(@times, a(:,:,:,1,2), b(:,:,:,2,2)) + bsxfun(@times, a(:,:,:,1,3), b(:,:,:,3,2));
    c(:,:,:,1,3) = bsxfun(@times, a(:,:,:,1,1), b(:,:,:,1,3)) + bsxfun(@times, a(:,:,:,1,2), b(:,:,:,2,3)) + bsxfun(@times, a(:,:,:,1,3), b(:,:,:,3,3));
    
    c(:,:,:,2,1) = bsxfun(@times, a(:,:,:,2,1), b(:,:,:,1,1)) + bsxfun(@times, a(:,:,:,2,2), b(:,:,:,2,1)) + bsxfun(@times, a(:,:,:,2,3), b(:,:,:,3,1));
    c(:,:,:,2,2) = bsxfun(@times, a(:,:,:,2,1), b(:,:,:,1,2)) + bsxfun(@times, a(:,:,:,2,2), b(:,:,:,2,2)) + bsxfun(@times, a(:,:,:,2,3), b(:,:,:,3,2));
    c(:,:,:,2,3) = bsxfun(@times, a(:,:,:,2,1), b(:,:,:,1,3)) + bsxfun(@times, a(:,:,:,2,2), b(:,:,:,2,3)) + bsxfun(@times, a(:,:,:,2,3), b(:,:,:,3,3));
    
    c(:,:,:,3,1) = bsxfun(@times, a(:,:,:,3,1), b(:,:,:,1,1)) + bsxfun(@times, a(:,:,:,3,2), b(:,:,:,2,1)) + bsxfun(@times, a(:,:,:,3,3), b(:,:,:,3,1));
    c(:,:,:,3,2) = bsxfun(@times, a(:,:,:,3,1), b(:,:,:,1,2)) + bsxfun(@times, a(:,:,:,3,2), b(:,:,:,2,2)) + bsxfun(@times, a(:,:,:,3,3), b(:,:,:,3,2));
    c(:,:,:,3,3) = bsxfun(@times, a(:,:,:,3,1), b(:,:,:,1,3)) + bsxfun(@times, a(:,:,:,3,2), b(:,:,:,2,3)) + bsxfun(@times, a(:,:,:,3,3), b(:,:,:,3,3));
end

function c = pw3_masb(a, b)
% a - 3x3 tensor field
% b - 3x3 symmetric tensor field = 6d vector field
% c - 3x3 tensor field
% Returns pointwise a * b
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim 3 3], 'like', a);
    
    ind = symIndices(size(b, 4));
    
    c(:,:,:,1,1) = bsxfun(@times, a(:,:,:,1,1), b(:,:,:,ind(1,1))) + bsxfun(@times, a(:,:,:,1,2), b(:,:,:,ind(2,1))) + bsxfun(@times, a(:,:,:,1,3), b(:,:,:,ind(3,1)));
    c(:,:,:,1,2) = bsxfun(@times, a(:,:,:,1,1), b(:,:,:,ind(1,2))) + bsxfun(@times, a(:,:,:,1,2), b(:,:,:,ind(2,2))) + bsxfun(@times, a(:,:,:,1,3), b(:,:,:,ind(3,2)));
    c(:,:,:,1,3) = bsxfun(@times, a(:,:,:,1,1), b(:,:,:,ind(1,3))) + bsxfun(@times, a(:,:,:,1,2), b(:,:,:,ind(2,3))) + bsxfun(@times, a(:,:,:,1,3), b(:,:,:,ind(3,3)));
    
    c(:,:,:,2,1) = bsxfun(@times, a(:,:,:,2,1), b(:,:,:,ind(1,1))) + bsxfun(@times, a(:,:,:,2,2), b(:,:,:,ind(2,1))) + bsxfun(@times, a(:,:,:,2,3), b(:,:,:,ind(3,1)));
    c(:,:,:,2,2) = bsxfun(@times, a(:,:,:,2,1), b(:,:,:,ind(1,2))) + bsxfun(@times, a(:,:,:,2,2), b(:,:,:,ind(2,2))) + bsxfun(@times, a(:,:,:,2,3), b(:,:,:,ind(3,2)));
    c(:,:,:,2,3) = bsxfun(@times, a(:,:,:,2,1), b(:,:,:,ind(1,3))) + bsxfun(@times, a(:,:,:,2,2), b(:,:,:,ind(2,3))) + bsxfun(@times, a(:,:,:,2,3), b(:,:,:,ind(3,3)));
    
    c(:,:,:,3,1) = bsxfun(@times, a(:,:,:,3,1), b(:,:,:,ind(1,1))) + bsxfun(@times, a(:,:,:,3,2), b(:,:,:,ind(2,1))) + bsxfun(@times, a(:,:,:,3,3), b(:,:,:,ind(3,1)));
    c(:,:,:,3,2) = bsxfun(@times, a(:,:,:,3,1), b(:,:,:,ind(1,2))) + bsxfun(@times, a(:,:,:,3,2), b(:,:,:,ind(2,2))) + bsxfun(@times, a(:,:,:,3,3), b(:,:,:,ind(3,2)));
    c(:,:,:,3,3) = bsxfun(@times, a(:,:,:,3,1), b(:,:,:,ind(1,3))) + bsxfun(@times, a(:,:,:,3,2), b(:,:,:,ind(2,3))) + bsxfun(@times, a(:,:,:,3,3), b(:,:,:,ind(3,3)));
end

function c = pw3_samb(a, b)
% a - 3x3 symmetric tensor field = 6d vector field
% b - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise a * b
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim 3 3], 'like', a);
    
    ind = symIndices(size(a, 4));
    
    c(:,:,:,1,1) = bsxfun(@times, a(:,:,:,ind(1,1)), b(:,:,:,1,1)) + bsxfun(@times, a(:,:,:,ind(1,2)), b(:,:,:,2,1)) + bsxfun(@times, a(:,:,:,ind(1,3)), b(:,:,:,3,1));
    c(:,:,:,1,2) = bsxfun(@times, a(:,:,:,ind(1,1)), b(:,:,:,1,2)) + bsxfun(@times, a(:,:,:,ind(1,2)), b(:,:,:,2,2)) + bsxfun(@times, a(:,:,:,ind(1,3)), b(:,:,:,3,2));
    c(:,:,:,1,3) = bsxfun(@times, a(:,:,:,ind(1,1)), b(:,:,:,1,3)) + bsxfun(@times, a(:,:,:,ind(1,2)), b(:,:,:,2,3)) + bsxfun(@times, a(:,:,:,ind(1,3)), b(:,:,:,3,3));
    
    c(:,:,:,2,1) = bsxfun(@times, a(:,:,:,ind(2,1)), b(:,:,:,1,1)) + bsxfun(@times, a(:,:,:,ind(2,2)), b(:,:,:,2,1)) + bsxfun(@times, a(:,:,:,ind(2,3)), b(:,:,:,3,1));
    c(:,:,:,2,2) = bsxfun(@times, a(:,:,:,ind(2,1)), b(:,:,:,1,2)) + bsxfun(@times, a(:,:,:,ind(2,2)), b(:,:,:,2,2)) + bsxfun(@times, a(:,:,:,ind(2,3)), b(:,:,:,3,2));
    c(:,:,:,2,3) = bsxfun(@times, a(:,:,:,ind(2,1)), b(:,:,:,1,3)) + bsxfun(@times, a(:,:,:,ind(2,2)), b(:,:,:,2,3)) + bsxfun(@times, a(:,:,:,ind(2,3)), b(:,:,:,3,3));
    
    c(:,:,:,3,1) = bsxfun(@times, a(:,:,:,ind(3,1)), b(:,:,:,1,1)) + bsxfun(@times, a(:,:,:,ind(3,2)), b(:,:,:,2,1)) + bsxfun(@times, a(:,:,:,ind(3,3)), b(:,:,:,3,1));
    c(:,:,:,3,2) = bsxfun(@times, a(:,:,:,ind(3,1)), b(:,:,:,1,2)) + bsxfun(@times, a(:,:,:,ind(3,2)), b(:,:,:,2,2)) + bsxfun(@times, a(:,:,:,ind(3,3)), b(:,:,:,3,2));
    c(:,:,:,3,3) = bsxfun(@times, a(:,:,:,ind(3,1)), b(:,:,:,1,3)) + bsxfun(@times, a(:,:,:,ind(3,2)), b(:,:,:,2,3)) + bsxfun(@times, a(:,:,:,ind(3,3)), b(:,:,:,3,3));
end

function c = pw3_sasb(a, b)
% a - 3x3 symmetric tensor field = 6d vector field
% b - 3x3 symmetric tensor field = 6d vector field
% c - 3x3 tensor field
% Returns pointwise a * b
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim 3 3], 'like', a);
    
    ind = symIndices(size(a, 4));
    
    c(:,:,:,1,1) = bsxfun(@times, a(:,:,:,ind(1,1)), b(:,:,:,ind(1,1))) + bsxfun(@times, a(:,:,:,ind(1,2)), b(:,:,:,ind(2,1))) + bsxfun(@times, a(:,:,:,ind(1,3)), b(:,:,:,ind(3,1)));
    c(:,:,:,1,2) = bsxfun(@times, a(:,:,:,ind(1,1)), b(:,:,:,ind(1,2))) + bsxfun(@times, a(:,:,:,ind(1,2)), b(:,:,:,ind(2,2))) + bsxfun(@times, a(:,:,:,ind(1,3)), b(:,:,:,ind(3,2)));
    c(:,:,:,1,3) = bsxfun(@times, a(:,:,:,ind(1,1)), b(:,:,:,ind(1,3))) + bsxfun(@times, a(:,:,:,ind(1,2)), b(:,:,:,ind(2,3))) + bsxfun(@times, a(:,:,:,ind(1,3)), b(:,:,:,ind(3,3)));
    
    c(:,:,:,2,1) = bsxfun(@times, a(:,:,:,ind(2,1)), b(:,:,:,ind(1,1))) + bsxfun(@times, a(:,:,:,ind(2,2)), b(:,:,:,ind(2,1))) + bsxfun(@times, a(:,:,:,ind(2,3)), b(:,:,:,ind(3,1)));
    c(:,:,:,2,2) = bsxfun(@times, a(:,:,:,ind(2,1)), b(:,:,:,ind(1,2))) + bsxfun(@times, a(:,:,:,ind(2,2)), b(:,:,:,ind(2,2))) + bsxfun(@times, a(:,:,:,ind(2,3)), b(:,:,:,ind(3,2)));
    c(:,:,:,2,3) = bsxfun(@times, a(:,:,:,ind(2,1)), b(:,:,:,ind(1,3))) + bsxfun(@times, a(:,:,:,ind(2,2)), b(:,:,:,ind(2,3))) + bsxfun(@times, a(:,:,:,ind(2,3)), b(:,:,:,ind(3,3)));
    
    c(:,:,:,3,1) = bsxfun(@times, a(:,:,:,ind(3,1)), b(:,:,:,ind(1,1))) + bsxfun(@times, a(:,:,:,ind(3,2)), b(:,:,:,ind(2,1))) + bsxfun(@times, a(:,:,:,ind(3,3)), b(:,:,:,ind(3,1)));
    c(:,:,:,3,2) = bsxfun(@times, a(:,:,:,ind(3,1)), b(:,:,:,ind(1,2))) + bsxfun(@times, a(:,:,:,ind(3,2)), b(:,:,:,ind(2,2))) + bsxfun(@times, a(:,:,:,ind(3,3)), b(:,:,:,ind(3,2)));
    c(:,:,:,3,3) = bsxfun(@times, a(:,:,:,ind(3,1)), b(:,:,:,ind(1,3))) + bsxfun(@times, a(:,:,:,ind(3,2)), b(:,:,:,ind(2,3))) + bsxfun(@times, a(:,:,:,ind(3,3)), b(:,:,:,ind(3,3)));
end

function c = pw3_mavb(a, b)
% a - 3x3 tensor field
% b - 3d  vector field
% b - 3d  vector field
% Returns pointwise a * b
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim 3], 'like', a);
    
    c(:,:,:,1) = bsxfun(@times, a(:,:,:,1,1), b(:,:,:,1)) + bsxfun(@times, a(:,:,:,1,2), b(:,:,:,2)) + bsxfun(@times, a(:,:,:,1,3), b(:,:,:,3));
    c(:,:,:,2) = bsxfun(@times, a(:,:,:,2,1), b(:,:,:,1)) + bsxfun(@times, a(:,:,:,2,2), b(:,:,:,2)) + bsxfun(@times, a(:,:,:,2,3), b(:,:,:,3));
    c(:,:,:,3) = bsxfun(@times, a(:,:,:,3,1), b(:,:,:,1)) + bsxfun(@times, a(:,:,:,3,2), b(:,:,:,2)) + bsxfun(@times, a(:,:,:,3,3), b(:,:,:,3));
end

function c = pw3_savb(a, b)
% a - 3x3 symmetric tensor field = 6d vector field
% b - 3d  vector field
% b - 3d  vector field
% Returns pointwise a * b
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim 3], 'like', a);
    
    ind = symIndices(size(a, 4));
    
    c(:,:,:,1) = bsxfun(@times, a(:,:,:,ind(1,1)), b(:,:,:,1)) + bsxfun(@times, a(:,:,:,ind(1,2)), b(:,:,:,2)) + bsxfun(@times, a(:,:,:,ind(1,3)), b(:,:,:,3));
    c(:,:,:,2) = bsxfun(@times, a(:,:,:,ind(2,1)), b(:,:,:,1)) + bsxfun(@times, a(:,:,:,ind(2,2)), b(:,:,:,2)) + bsxfun(@times, a(:,:,:,ind(2,3)), b(:,:,:,3));
    c(:,:,:,3) = bsxfun(@times, a(:,:,:,ind(3,1)), b(:,:,:,1)) + bsxfun(@times, a(:,:,:,ind(3,2)), b(:,:,:,2)) + bsxfun(@times, a(:,:,:,ind(3,3)), b(:,:,:,3));
end

function c = pw3_vamb(a, b)
% a - 3d  vector field
% b - 3x3 tensor field
% b - 3d  vector field
% Returns pointwise a.' * b
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim 3], 'like', a);
    
    c(:,:,:,1) = bsxfun(@times, a(:,:,:,1), b(:,:,:,1,1)) + bsxfun(@times, a(:,:,:,2), b(:,:,:,2,1)) + bsxfun(@times, a(:,:,:,3), b(:,:,:,3,1));
    c(:,:,:,2) = bsxfun(@times, a(:,:,:,1), b(:,:,:,1,2)) + bsxfun(@times, a(:,:,:,2), b(:,:,:,2,2)) + bsxfun(@times, a(:,:,:,3), b(:,:,:,3,2));
    c(:,:,:,3) = bsxfun(@times, a(:,:,:,1), b(:,:,:,1,3)) + bsxfun(@times, a(:,:,:,2), b(:,:,:,2,3)) + bsxfun(@times, a(:,:,:,3), b(:,:,:,3,3));
end

function c = pw3_vasb(a, b)
% a - 3d  vector field
% b - 3x3 symmetric tensor field = 6d vector field
% b - 3d  vector field
% Returns pointwise a.' * b
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim 3], 'like', a);
    
    ind = symIndices(size(b, 4));
    
    c(:,:,:,1) = bsxfun(@times, a(:,:,:,1), b(:,:,:,ind(1,1))) + bsxfun(@times, a(:,:,:,2), b(:,:,:,ind(2,1))) + bsxfun(@times, a(:,:,:,3), b(:,:,:,ind(3,1)));
    c(:,:,:,2) = bsxfun(@times, a(:,:,:,1), b(:,:,:,ind(1,2))) + bsxfun(@times, a(:,:,:,2), b(:,:,:,ind(2,2))) + bsxfun(@times, a(:,:,:,3), b(:,:,:,ind(3,2)));
    c(:,:,:,3) = bsxfun(@times, a(:,:,:,1), b(:,:,:,ind(1,3))) + bsxfun(@times, a(:,:,:,2), b(:,:,:,ind(2,3))) + bsxfun(@times, a(:,:,:,3), b(:,:,:,ind(3,3)));
end

function c = pw3_vavb(a, b)
% a - 3d vector field
% b - 3d vector field
% b - 3d vector field
% Returns pointwise a.' * b
    c = bsxfun(@times, a(:,:,:,1), b(:,:,:,1)) + bsxfun(@times, a(:,:,:,2), b(:,:,:,2)) + bsxfun(@times, a(:,:,:,3), b(:,:,:,3));
end

function c = pw3_tmamb(a, b)
% a - 3x3 tensor field
% b - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise a.' * b
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim 3 3], 'like', a);
    
    c(:,:,:,1,1) = bsxfun(@times, a(:,:,:,1,1), b(:,:,:,1,1)) + bsxfun(@times, a(:,:,:,2,1), b(:,:,:,2,1)) + bsxfun(@times, a(:,:,:,3,1), b(:,:,:,3,1));
    c(:,:,:,1,2) = bsxfun(@times, a(:,:,:,1,1), b(:,:,:,1,2)) + bsxfun(@times, a(:,:,:,2,1), b(:,:,:,2,2)) + bsxfun(@times, a(:,:,:,3,1), b(:,:,:,3,2));
    c(:,:,:,1,3) = bsxfun(@times, a(:,:,:,1,1), b(:,:,:,1,3)) + bsxfun(@times, a(:,:,:,2,1), b(:,:,:,2,3)) + bsxfun(@times, a(:,:,:,3,1), b(:,:,:,3,3));
    
    c(:,:,:,2,1) = bsxfun(@times, a(:,:,:,1,2), b(:,:,:,1,1)) + bsxfun(@times, a(:,:,:,2,2), b(:,:,:,2,1)) + bsxfun(@times, a(:,:,:,3,2), b(:,:,:,3,1));
    c(:,:,:,2,2) = bsxfun(@times, a(:,:,:,1,2), b(:,:,:,1,2)) + bsxfun(@times, a(:,:,:,2,2), b(:,:,:,2,2)) + bsxfun(@times, a(:,:,:,3,2), b(:,:,:,3,2));
    c(:,:,:,2,3) = bsxfun(@times, a(:,:,:,1,2), b(:,:,:,1,3)) + bsxfun(@times, a(:,:,:,2,2), b(:,:,:,2,3)) + bsxfun(@times, a(:,:,:,3,2), b(:,:,:,3,3));
    
    c(:,:,:,3,1) = bsxfun(@times, a(:,:,:,1,3), b(:,:,:,1,1)) + bsxfun(@times, a(:,:,:,2,3), b(:,:,:,2,1)) + bsxfun(@times, a(:,:,:,3,3), b(:,:,:,3,1));
    c(:,:,:,3,2) = bsxfun(@times, a(:,:,:,1,3), b(:,:,:,1,2)) + bsxfun(@times, a(:,:,:,2,3), b(:,:,:,2,2)) + bsxfun(@times, a(:,:,:,3,3), b(:,:,:,3,2));
    c(:,:,:,3,3) = bsxfun(@times, a(:,:,:,1,3), b(:,:,:,1,3)) + bsxfun(@times, a(:,:,:,2,3), b(:,:,:,2,3)) + bsxfun(@times, a(:,:,:,3,3), b(:,:,:,3,3));
end

function c = pw3_tma(a)
% a - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise a.'
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) 3 3], 'like', a);
    
    c(:,:,:,1,1) = a(:,:,:,1,1);
    c(:,:,:,1,2) = a(:,:,:,2,1);
    c(:,:,:,1,3) = a(:,:,:,3,1);
    
    c(:,:,:,2,1) = a(:,:,:,1,2);
    c(:,:,:,2,2) = a(:,:,:,2,2);
    c(:,:,:,2,3) = a(:,:,:,3,2);
    
    c(:,:,:,3,1) = a(:,:,:,1,3);
    c(:,:,:,3,2) = a(:,:,:,2,3);
    c(:,:,:,3,3) = a(:,:,:,3,3);
end

function c = pw3_ima(a)
% a - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise inv(a)
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) 3 3], 'like', a);
    
    % - Cofactors
    
    c(:,:,:,1,1) = a(:,:,:,2,2) .* a(:,:,:,3,3) - a(:,:,:,2,3) .* a(:,:,:,3,2);
    c(:,:,:,1,2) = a(:,:,:,1,3) .* a(:,:,:,3,2) - a(:,:,:,1,2) .* a(:,:,:,3,3);
    c(:,:,:,1,3) = a(:,:,:,1,2) .* a(:,:,:,2,3) - a(:,:,:,1,3) .* a(:,:,:,2,2);
    
    c(:,:,:,2,1) = a(:,:,:,2,3) .* a(:,:,:,3,1) - a(:,:,:,2,1) .* a(:,:,:,3,3);
    c(:,:,:,2,2) = a(:,:,:,1,1) .* a(:,:,:,3,3) - a(:,:,:,1,3) .* a(:,:,:,3,1);
    c(:,:,:,2,3) = a(:,:,:,1,3) .* a(:,:,:,2,1) - a(:,:,:,1,1) .* a(:,:,:,2,3);
    
    c(:,:,:,3,1) = a(:,:,:,2,1) .* a(:,:,:,3,2) - a(:,:,:,2,2) .* a(:,:,:,3,1);
    c(:,:,:,3,2) = a(:,:,:,1,2) .* a(:,:,:,3,1) - a(:,:,:,1,1) .* a(:,:,:,3,2);
    c(:,:,:,3,3) = a(:,:,:,1,1) .* a(:,:,:,2,2) - a(:,:,:,1,2) .* a(:,:,:,2,1);
    
    % - Determinant
    
    c = bsxfun(@rdivide, c, a(:,:,:,1,1) .* c(:,:,:,1,1) + ...
                            a(:,:,:,1,2) .* c(:,:,:,2,1) + ...
                            a(:,:,:,1,3) .* c(:,:,:,3,1) );
end

function c = pw3_isa(a)
% a - 3x3 symmetric tensor field = 6d vector field
% c - 3x3 symmetric tensor field = 6d vector field
% Returns pointwise inv(a)
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) 6], 'like', a);
    
    ind = symIndices(size(a, 4));
    
    % - Cofactors
    
    c(:,:,:,ind(1,1)) = a(:,:,:,ind(2,2)) .* a(:,:,:,ind(3,3)) - a(:,:,:,ind(2,3)) .* a(:,:,:,ind(3,2));
    c(:,:,:,ind(1,2)) = a(:,:,:,ind(1,3)) .* a(:,:,:,ind(3,2)) - a(:,:,:,ind(1,2)) .* a(:,:,:,ind(3,3));
    c(:,:,:,ind(1,3)) = a(:,:,:,ind(1,2)) .* a(:,:,:,ind(2,3)) - a(:,:,:,ind(1,3)) .* a(:,:,:,ind(2,2));
    
    c(:,:,:,ind(2,2)) = a(:,:,:,ind(1,1)) .* a(:,:,:,ind(3,3)) - a(:,:,:,ind(1,3)) .* a(:,:,:,ind(3,1));
    
    c(:,:,:,ind(3,2)) = a(:,:,:,ind(1,2)) .* a(:,:,:,ind(3,1)) - a(:,:,:,ind(1,1)) .* a(:,:,:,ind(3,2));
    
    % - Determinant
    
    c = bsxfun(@rdivide, c, a(:,:,:,ind(1,1)) .* c(:,:,:,ind(1,1)) + ...
                            a(:,:,:,ind(1,2)) .* c(:,:,:,ind(2,1)) + ...
                            a(:,:,:,ind(1,3)) .* c(:,:,:,ind(3,1)) );
end

function c = pw3_dma(a)
% a - 3x3 tensor field
% c - scalar field
% Returns pointwise det(a)
    c = a(:,:,:,1,1) .* ( a(:,:,:,2,2) .* a(:,:,:,3,3) - a(:,:,:,2,3) .* a(:,:,:,3,2) ) + ...
        a(:,:,:,1,2) .* ( a(:,:,:,2,3) .* a(:,:,:,3,1) - a(:,:,:,2,1) .* a(:,:,:,3,3) ) + ...
        a(:,:,:,1,3) .* ( a(:,:,:,2,1) .* a(:,:,:,3,2) - a(:,:,:,2,2) .* a(:,:,:,3,1) );
end


function c = pw3_dsa(a)
% a - 3x3 symmetric tensor field = 6d vector field
% c - scalar field
% Returns pointwise det(a)
    ind = symIndices(size(a, 4));
    c = a(:,:,:,ind(1,1)) .* ( a(:,:,:,ind(2,2)) .* a(:,:,:,ind(3,3)) - a(:,:,:,ind(2,3)) .* a(:,:,:,ind(3,2)) ) + ...
        a(:,:,:,ind(1,2)) .* ( a(:,:,:,ind(2,3)) .* a(:,:,:,ind(3,1)) - a(:,:,:,ind(2,1)) .* a(:,:,:,ind(3,3)) ) + ...
        a(:,:,:,ind(1,3)) .* ( a(:,:,:,ind(2,1)) .* a(:,:,:,ind(3,2)) - a(:,:,:,ind(2,2)) .* a(:,:,:,ind(3,1)) );
end


%% === pointwise ==========================================================

function c = pointwise(a, varargin)
% FORMAT c = spm_matcomp('Pointwise', a, (b), (op), (sym))
% a  - [nx nz nz A (B)] ND-array of vectors/matrices
% b  - [nx nz nz C (D)] ND-array of vectors/matrices
% op - Operation(s) to apply to an element of a before multiplying it 
%      with an element of b:
%      't' (transpose).
% sym - If true, at least one of a and b is sparse symmetric.
%
% Performs matrix operations at each point of a tensor field.
% If b is provided, at each point:
%    c = op(a) * b
% Else, at each point:
%    c = op(a)
%
% Note that
% - if one (or both) entry is a vector, the adequate transpose
%   operation is usually automatically detected.
% - a and b can be fields of sparse symmetric matrices. In this case, they
%   have the form of a vector field, and sym must be set to true.
% - Implicit expansion of the lattice is performed when one of the input
%   dimensions is singular.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

    % --- Default arguments
    b = [];
    issym = false;
    op = '';
    while ~isempty(varargin)
        if ischar(varargin{1})
            op = varargin{1};
        elseif isscalar(varargin{1})
            issym = varargin{1};
        else
            b = varargin{1};
        end
        varargin = varargin(2:end);
    end
    binary = ~isempty(b);
    
    % --- Check dimensions
    ismat = struct;
    dim   = struct;
    sym   = struct;
    dim.a = [size(a) 1 1];
    if ~binary
        sym.a = issym;
    end
    if binary
        dim.b = [size(b) 1 1];
        if issym
            if size(a,5) == 1 && size(b,5) > 1
                sym.a = true;
                sym.b = false;
            elseif size(a,5) > 1 && size(b,5) == 1
                sym.b = true;
                sym.a = false;
            else
                sym.a = false;
                sym.b = false;
                if size(a,4) >=  size(b,4)
                    sym.a = true;
                end
                if size(b,4) >= size(a,4)
                    sym.b = true;
                end
            end
        else
            sym.a = false;
            sym.b = false;
        end
    end
    % The corresponding dimensions of A and B must be equal to each other, 
    % or equal to one.
    if binary
        dima = dim.a(1:3);
        dimb = dim.b(1:3);
        same = dima == dimb;
        diff = min(dima(~same), dimb(~same));
        if  any(diff ~= 1)
            error(['The corresponding lattice dimensions of A and B ' ...
                   'must be equal to each other, or equal to one.']);
        end
    end
    ismat.a = sym.a || size(a, 5) > 1;
    ismat.b = sym.b || size(b, 5) > 1;
    
    % --- Deal with cases
    switch lower(op)
        case ''
            if binary && ismat.a && ismat.b
                if sym.a && sym.b
                    c = pw_sasb(a, b);
                elseif sym.a
                    c = pw_samb(a, b);
                elseif sym.b
                    c = pw_masb(a, b);
                else
                    c = pw_mamb(a, b);
                end
            elseif binary && ismat.a && ~ismat.b
                if sym.a
                    c = pw_savb(a, b);
                else
                    c = pw_mavb(a, b);
                end
            elseif binary && ~ismat.a && ismat.b
                if sym.b
                    c = pw_vasb(a, b);
                else
                    c = pw_vamb(a, b);
                end
            elseif  binary && ~ismat.a && ~ismat.b
                c = pw_vavb(a, b);
            else
                c = a;
            end
        case 't'
            if binary && ismat.a && ismat.b
                if sym.a && sym.b
                    c = pw_sasb(a, b);
                elseif sym.b
                    c = pw_masb(a, b);
                elseif sym.a
                    c = pw_samb(a, b);
                else
                    c = pw_tmamb(a, b);
                end
            elseif binary && ismat.a && ~ismat.b
                if sym.a
                    c = pw_savb(a, b);
                else
                    c = pw_tmavb(a, b);
                end
            elseif binary && ~ismat.a && ismat.b
                if sym.b
                    c = pw_vasb(a, b);
                else
                    c = pw_vamb(a, b);
                end
            elseif  binary && ~ismat.a && ~ismat.b
                c = pw_vavb(a, b);
            elseif ismat.a
                if sym.a
                    c = a;
                else
                    c = pw_tma(a);
                end
            else
                c = a;
            end
        otherwise
            error('The provided operation is not supported')
    end
end

function c = pw_mamb(a, b)
% a - NxM tensor field
% b - MxK tensor field
% c - NxK tensor field
% Returns pointwise a * b

    if size(a,5) ~= size(b,4)
        error('Matrix dimensions not consistant: %d and %d', size(a,5), size(b,4))
    end
    
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim size(a,4) size(b,5)], 'like', a);
    
    for i=1:size(a,4)
        for j=1:size(b,5)
            for k=1:size(a,5)
                c(:,:,:,i,j) = c(:,:,:,i,j) + bsxfun(@times, a(:,:,:,i,k), b(:,:,:,k,j));
            end
        end
    end
end

function c = pw_masb(a, b)
% a - NxM tensor field
% b - MxM symmetric tensor field = 6d vector field
% c - NxM tensor field
% Returns pointwise a * b

    [ind, n] = symIndices(size(b, 4));
    if size(a,5) ~= n
        error('Matrix dimensions not consistant: %d and %d', size(a,5), n)
    end
    
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim size(a,4) n], 'like', a);
    
    for i=1:size(a,4)
        for j=1:n
            for k=1:n
                c(:,:,:,i,j) = c(:,:,:,i,j) + bsxfun(@times, a(:,:,:,i,k), b(:,:,:,ind(k,j)));
            end
        end
    end
end

function c = pw_samb(a, b)
% a - NxN symmetric tensor field = 6d vector field
% b - NxM tensor field
% c - NxM tensor field
% Returns pointwise a * b

    ind = symIndices(size(a, 4));
    if n ~= size(b,4)
        error('Matrix dimensions not consistant: %d and %d', n, size(b,4))
    end
    
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim n size(b,5)], 'like', a);
    
    for i=1:n
        for j=1:size(b,5)
            for k=1:n
                c(:,:,:,i,j) = bsxfun(@times, a(:,:,:,ind(i,k)), b(:,:,:,k,j));
            end
        end
    end
end

function c = pw_sasb(a, b)
% a - NxN symmetric tensor field = 6d vector field
% b - NxN symmetric tensor field = 6d vector field
% c - NxN tensor field
% Returns pointwise a * b

    
    [inda, na] = symIndices(size(a, 4));
    [indb, nb] = symIndices(size(b, 4));
    if na ~= nb
        error('Matrix dimensions not consistant: %d and %d', na, nb)
    end
    
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim na na], 'like', a);
    
    for i=1:na
        for j=1:na
            for k=1:na
                c(:,:,:,i,j) = c(:,:,:,i,j) + bsxfun(@times, a(:,:,:,inda(i,k)), b(:,:,:,indb(k,j)));
            end
        end
    end
end

function c = pw_mavb(a, b)
% a - NxM tensor field
% b - Md  vector field
% c - Nd  vector field
% Returns pointwise a * b

    transpose = false;
    n = size(a,4);
    if size(a,5) ~= size(b,4)
        if size(a,4) ~= size(b,4)
            error('Matrix dimensions not consistant: %d and %d', size(a,5), size(b,4))
        else
            n = size(a,5);
            transpose = true;
        end
    end
    
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim n], 'like', a);
    
    if ~transpose
        for i=1:size(a,4)
            for k=1:size(b,4)
                c(:,:,:,i) = c(:,:,:,i) + bsxfun(@times, a(:,:,:,i,k), b(:,:,:,k));
            end
        end
    else
        for i=1:size(a,5)
            for k=1:size(b,4)
                c(:,:,:,i) = c(:,:,:,i) + bsxfun(@times, a(:,:,:,k,i), b(:,:,:,k));
            end
        end
    end
end

function c = pw_savb(a, b)
% a - NxN symmetric tensor field = 6d vector field
% b - Nd  vector field
% c - Nd  vector field
% Returns pointwise a * b

    [ind, n] = symIndices(size(a, 4));
    if n ~= size(b,4)
        error('Matrix dimensions not consistant: %d and %d', n, size(b,4))
    end
    
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim n], 'like', a);
    
    for i=1:n
        for k=1:n
            c(:,:,:,i) = c(:,:,:,i) + bsxfun(@times, a(:,:,:,ind(i,k)), b(:,:,:,k));
        end
    end
end

function c = pw_vamb(a, b)
% a - Nd  vector field
% b - NxM tensor field
% b - Md  vector field
% Returns pointwise a.' * b

    transpose = false;
    n = size(b,5);
    if size(a,4) ~= size(b,4)
        if size(a,4) ~= size(b,5)
            error('Matrix dimensions not consistant: %d and %d', size(a,4), size(b,4))
        else
            n = size(b,4);
            transpose = true;
        end
    end
    
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim n], 'like', a);
    
    if ~transpose
        for i=1:size(b,5)
            for k=1:size(a,4)
                c(:,:,:,i) = c(:,:,:,i) + bsxfun(@times, a(:,:,:,k), b(:,:,:,k,i));
            end
        end
    else
        for i=1:size(b,4)
            for k=1:size(a,4)
                c(:,:,:,i) = c(:,:,:,i) + bsxfun(@times, a(:,:,:,k), b(:,:,:,i,k));
            end
        end
    end
end

function c = pw_vasb(a, b)
% a - Nd  vector field
% b - NxN symmetric tensor field = 6d vector field
% b - Nd  vector field
% Returns pointwise a.' * b

    [ind, n] = symIndices(size(b, 4));
    if size(a,4) ~= n
        error('Matrix dimensions not consistant: %d and %d', size(a,4), n)
    end
    
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim n], 'like', a);
    
    for i=1:n
        for k=1:n
            c(:,:,:,i) = c(:,:,:,i) + bsxfun(@times, a(:,:,:,k), b(:,:,:,ind(k,i)));
        end
    end
end

function c = pw_vavb(a, b)
% a - Nd vector field
% b - Nd vector field
% c -    scalar field
% Returns pointwise a.' * b
    
    if size(a,4) ~= size(b,4)
        error('Matrix dimensions not consistant: %d and %d', size(a,4), size(b,4))
    end
    
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros(dim, 'like', a);
    
    for k=1:size(a,4)
        c = c + bsxfun(@times, a(:,:,:,k), b(:,:,:,k));
    end
end

function c = pw_tmamb(a, b)
% a - MxN tensor field
% b - MxK tensor field
% c - NxK tensor field
% Returns pointwise a.' * b

    if size(a,4) ~= size(b,4)
        error('Matrix dimensions not consistant: %d and %d', size(a,4), size(b,4))
    end
    
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim size(a,5) size(b,5)], 'like', a);
    
    for i=1:size(a,5)
        for j=size(b,5)
            for k=1:size(a,4)
                c(:,:,:,i,j) = c(:,:,:,i,j) + bsxfun(@times, a(:,:,:,k,i), b(:,:,:,k,j));
            end
        end
    end
end

function c = pw_tmavb(a, b)
% a - MxN tensor field
% b - Md  vector field
% c - Nd  vector field
% Returns pointwise a.' * b
    
    transpose = false;
    n = size(a,5);
    if size(a,4) ~= size(b,4)
        if size(a,5) ~= size(b,4)
            error('Matrix dimensions not consistant: %d and %d', size(a,4), size(b,4))
        else
            n = size(a,4);
            transpose = true;
        end
    end
    
    dima = [size(a) 1 1];
    dimb = [size(b) 1 1];
    dim  = max(dima(1:3), dimb(1:3));
    c = zeros([dim n], 'like', numeric(a(1)));
    
    if ~transpose
        for i=1:size(a,5)
            for k=1:size(b,4)
                c(:,:,:,i) = c(:,:,:,i) + bsxfun(@times, a(:,:,:,k,i), b(:,:,:,k));
            end
        end
    else
        for i=1:size(a,4)
            for k=1:size(b,4)
                c(:,:,:,i) = c(:,:,:,i) + bsxfun(@times, a(:,:,:,i,k), b(:,:,:,k));
            end
        end
    end
end

function c = pw_tma(a)
% a - MxN tensor field
% c - NxM tensor field
% Returns pointwise a.' 

    dim = [size(a) 1 1];
    c = zeros([dim(1:3) size(a,5) size(a,4)], 'like', a);
    
    for i=1:size(a,5)
        for j=1:size(a,4)
            c(:,:,:,i,j) = a(:,:,:,j,i);
        end
    end
end

%==========================================================================
function M = threshold_eig(M)
% Stabilise matrix by thresholding eigenvalues
% FORMAT M = threshold_eig(M)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
[V,D1] = eig(M);
tol    = max(diag(D1))*eps('single');
D1     = diag(max(diag(D1),tol));
M      = real(V*D1*V');
end
%==========================================================================

%==========================================================================
function s = logsumexp(b,dm)
% Compute log-sum-exp trick: 
% https://en.wikipedia.org/wiki/LogSumExp#log-sum-exp_trick_for_log-domain_calculations
% FORMAT s = logsumexp(b,dm)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
B       = max(b,[],dm);
dms     = ones(1,ndims(b));
dms(dm) = size(b,dm);
b       = b - repmat(B,dms);
s       = B + log(nansum(exp(b),dm));
end
%==========================================================================

%==========================================================================
function Q = softmax(Q,prop)
% Compute soft-max
% FORMAT [Q,ll] = softmax(Q,scl)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<2, prop = 0; end

dm      = size(Q);
ndm     = numel(dm);
dm(:)   = 1;
dm(end) = numel(prop);

if prop ~= 0
    prop = reshape(prop,dm);
    Q   = bsxfun(@plus,Q,prop); % rescale log-probabilities
end

mxQ = max(Q,[],ndm);
Q   = exp(bsxfun(@minus,Q,mxQ));
sQ  = nansum(Q,ndm);

Q = bsxfun(@rdivide,Q,sQ);
end
%==========================================================================

% function a = softmax(a, scale)
%     if nargin < 2
%         scale = 1;
%     end
%     a = bsxfun(@plus, a, log(scale));     % rescale probabilities
%     a = bsxfun(@minus, a, max(a, [], 4)); % safe softmax -> avoid overflow
%     a = exp(a);                           % exponentiate
%     a = bsxfun(@rdivide, a, sum(a, 4));   % normalise
% end