% TODO: sliding boundary conditions

function varargout = spm_sparse(id, varargin)
%__________________________________________________________________________
% Collection of useful sparse matrices.
%
%--------------------------------------------------------------------------
% Smooth regularisation matrices.
% > L*v(:) is equivalent to spm_diffeo('vel2mom', v, ...)
%
% FORMAT L = spm_sparse('precision', 'diffeo', lat_dim, lat_vs, bnd)
% FORMAT K = spm_sparse('kernel',    'diffeo', lat_dim, lat_vs, bnd)
%
% FORMAT L = spm_sparse('precision', 'field',  lat_dim, lat_vs, bnd)
% FORMAT K = spm_sparse('kernel',    'field',  lat_dim, lat_vs, bnd)
%
%--------------------------------------------------------------------------
% Differential operators in matrix form.
% > Let v be the vectorized version of a scalar or vector field, 
%   and J be the jacobian operator, then J * v returns a vectorized  
%   version of the Jacobian matrix of v.
%
% FORMAT J = spm_sparse('jacobian',       lat_dim, lat_vs, vec_dim, bnd)
% FORMAT H = spm_sparse('hessian',        lat_dim, lat_vs, vec_dim, bnd)
% FORMAT L = spm_sparse('laplacian',      lat_dim, lat_vs, bnd)
% FORMAT D = spm_sparse('euclidian_dist', lat_dim, lat_vs)
% FORMAT A = spm_sparse('divergence',     lat_dim, lat_vs, bnd)
% FORMAT S = spm_sparse('symjac',         lat_dim, lat_vs, vec_dim, bnd)
% 
%--------------------------------------------------------------------------
% FORMAT help spm_sparse>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    switch id
        case 'precision'
            [varargout{1:nargout}] = precision(varargin{:});
        case 'kernel'
            [varargout{1:nargout}] = kernel(varargin{:});
        case 'jacobian'
            [varargout{1:nargout}] = jacobian(varargin{:});
        case 'hessian'
            [varargout{1:nargout}] = hessian(varargin{:});
        case 'euclidian_dist'
            [varargout{1:nargout}] = euclidian_dist(varargin{:});
        case 'laplacian'
            [varargout{1:nargout}] = laplacian(varargin{:});
        case 'divergence'
            [varargout{1:nargout}] = divergence(varargin{:});
        case 'symjac'
            [varargout{1:nargout}] = symjac(varargin{:});
        otherwise
            error('No function named %s\n', string(id));
    end
end

%% Precision matrices

function L = precision(mode, lat_dim, lat_vs, param, bnd)
% FORMAT L = spm_sparse('precision', mode, lat_dim, lat_vs, param, (bnd))
% mode     - 'diffeo'or 'field'
%            > If 'diffeo':  input must be a deformation/velocity field.
% lat_dim  - Dimensions of the field lattice ([nx ny ...])
% lat_vs   - Voxel size of the field lattice in mm/vox ([vx vy ...])
% param[1] - Weight of the absolute displacement penalty (mm^{-2})
% param[2] - Lambda parameter of the membrane energy.
% param[3] - Lambda parameter of the bending energy (mm^2).
% param[4] - [diffeo only] Mu parameter of the linear-elastic energy.
% param[5] - [diffeo only] Lambda parameter of the linear-elastic energy.
% bnd      - Boundary conditions:
%               * 0/'c'/'circulant'  (image wraps around) [default]
%               * 1/'n'/'neumann'    (mirror / boundary)
%               * 2/'d'/'dirichlet'  (zero outside FOV)
%
% Computes the inverse of the covariance matrix of the smooth prior 
% in the LDDMM framwork (i.e., the posdef, self-adjoint, differential 
% operator that defines a metric in the velocity space). It penalizes 
% absolute deformations and membrane, bending and linear-elastic energies. 
% See Ashburner, "A fast diffeomorphic image registration algorithm",
% NeuroImage (2007)
%
% Robustness is improved by averaging sum of square penalties obtained by
% different approximations of the Jacobian (right or left approximation),
% as done in spm_diffeo.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    im_dim = length(lat_dim);
    if nargin < 5
        bnd = 0;
        if nargin < 4
            param = [0.0001 0.001 0.2 0.05 0.2];
            if nargin < 3
                lat_vs = ones(1, im_dim);
            end
        end
    end
    param = [param 0 0 0];
    param = param(1:5);
    if strcmpi(mode, 'field')
        param(4) = 0;
        param(5) = 0;
    end
    
    lambda_abs      = param(1); % sqrt(lam) ~ 1/mm
    lambda_membrane = param(2); % sqrt(lam) ~ 1
    lambda_bending  = param(3); % sqrt(lam) ~ mm
    mu_elastic      = param(4); % sqrt(mu)  ~ 1
    lambda_elastic  = param(5); % sqrt(lam) ~ 1
    
    % Compute all different combination of Jacobian estimate
    % ([Jx+ Jy+], [Jx+ Jy-], [Jx- Jy+], [Jx- Jy-], ...)
    dirs = combdir(im_dim);
    
    if strcmpi(mode, 'field')
        im_dim = 1;
    end

    % Note that in the 'diffeo' case, I have to rescale the matrix by the
    % voxel size because displacement fields are expressed in voxels, not
    % millimetres.
    
    % Start with a zero matrix
    L = sparse(prod(lat_dim) * im_dim, prod(lat_dim) * im_dim);
    % Absolute displacement in each component (~= euclidian dist)
    if lambda_abs > 0
        E = speye(prod(lat_dim) * im_dim);
        if strcmpi(mode, 'diffeo')
            for j=1:im_dim
                sub = 1+(j-1)*prod(lat_dim):j*prod(lat_dim);
                E(sub,sub) = E(sub,sub) * lat_vs(j)^2;
            end
        end
        L = L + lambda_abs * E;
    end
    % Membrane energy
    if lambda_membrane > 0
        for i=1:numel(dirs)
            J = spm_sparse('jacobian', lat_dim, lat_vs, im_dim, dirs{i}, bnd);
            J = J' * J;
            if strcmpi(mode, 'diffeo')
                for j=1:im_dim
                    sub = 1+(j-1)*prod(lat_dim):j*prod(lat_dim);
                    J(sub,sub) = J(sub,sub) * lat_vs(j)^2;
                end
            end
            L = L + lambda_membrane / numel(dirs) * J;
        end
    end
    % Bending energy
    if lambda_bending > 0
        for i=1:numel(dirs)
            H = spm_sparse('hessian', lat_dim, lat_vs, im_dim, dirs{i}, bnd);
            H = H' * H;
            if strcmpi(mode, 'diffeo')
                for j=1:im_dim
                    sub = 1+(j-1)*prod(lat_dim):j*prod(lat_dim);
                    H(sub,sub) = H(sub,sub) * lat_vs(j)^2;
                end
            end
            L = L + lambda_bending / numel(dirs) * H;
        end
    end
    % Linear-elastic (symmetric part of the jacobian)
    if mu_elastic > 0
        for i=1:numel(dirs)
            S = spm_sparse('symjac', lat_dim, lat_vs, im_dim, dirs{i}, bnd);
            S = S' * S;
            if strcmpi(mode, 'diffeo')
                for j=1:im_dim
                    sub1 = 1+(j-1)*prod(lat_dim):j*prod(lat_dim);
                    for k=1:im_dim
                        sub2 = 1+(k-1)*prod(lat_dim):k*prod(lat_dim);
                        S(sub1,sub2) = S(sub1,sub2) * lat_vs(j) * lat_vs(k);
                    end
                end
            end
            L = L + 2 * mu_elastic / numel(dirs) * S;
        end
    end
    % Linear-elastic (divergence)
    if lambda_elastic > 0
        for i=1:numel(dirs)
            D = spm_sparse('divergence', lat_dim, lat_vs, dirs{i}, bnd);
            D = D' * D;
            if strcmpi(mode, 'diffeo')
                for j=1:im_dim
                    sub1 = 1+(j-1)*prod(lat_dim):j*prod(lat_dim);
                    for k=1:im_dim
                        sub2 = 1+(k-1)*prod(lat_dim):k*prod(lat_dim);
                        D(sub1,sub2) = D(sub1,sub2) * lat_vs(j) * lat_vs(k);
                    end
                end
            end
            L = L + lambda_elastic / numel(dirs) * D;
        end
    end
end

function CD = combdir(im_dim)
% Computes all possible combinations of directions for penalty operators
    sets = mat2cell(repmat([-1 1], im_dim, 1), ones(1,im_dim), 2);
    out = cell(im_dim, 1);
    [out{1:im_dim}] = ndgrid(sets{:});
    CD = cell(numel(out{1}), 1);
    for i=1:numel(out{1})
        combination = zeros(1, im_dim);
        for j=1:im_dim
            combination(j) = out{j}(i);
        end
        CD{i} = combination;
    end
end

function K = kernel(mode, L, lat_dim, lat_vs, bnd)
% FORMAT K = spm_sparse('kernel', mode, L, lat_dim)
% FORMAT K = spm_sparse('kernel', mode, param, lat_dim, lat_vs, bnd)
% mode      - 'diffeo' or 'field'
% L         - Regularization matrix (i.e. inverse of the prior covariance)
% param     - It is also possible to provide the penalty parameters instead
%             of the full matrix.
% lat_dim   - Dimension of the field lattice. Necessary because it is not
%             explicitely stored in L and is needed to extract the 
%             convolution operators from L.
% (lat_vs)  - Voxel size of the lattice (default: 1).
% (bnd)     - Boundary conditions:
%                * 0/'c'/'circulant'  (image wraps around) [default]
%                * 1/'n'/'neumann'    (mirror / boundary)
%                * 2/'d'/'dirichlet'  (zero outside FOV)
% K         - Direct convolution operators, equivalent to L.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    if nargin < 5
        bnd = 0;
        if nargin < 4
            lat_vs = ones(size(lat_dim));
        end
    end
    if size(L,1) == 1
        % Params were provided, not the matrix
        param = L;
        L = precision(mode, lat_dim, lat_vs, param, bnd);
    end
    im_dim = length(lat_dim);
    
    if strcmpi(mode, 'diffeo')
        K = zeros([lat_dim im_dim im_dim]);
        for i=1:im_dim
            K(:,:,:,:,i) = reshape(full(L(:,1+(i-1)*prod(lat_dim))), [lat_dim im_dim]);
        end
    else
        K = reshape(full(L(:,1)), lat_dim);
    end
end

%% Matricial operators

function J = jacobian(lat_dim, lat_vs, vec_dim, dir, bnd)
% FORMAT J = spm_sparse('jacobian', lat_dim, lat_vs, vec_dim, dir, bnd)
% lat_dim   - Dimensions of the lattice ([nx ny ...]).
% (lat_vs)  - Voxel size of the lattice (default: 1).
% (vec_dim) - Dimension of the vector space (default: 1).
% (dir)     - Approximation used in each direction ([dx dy ...]):
%             * -1: Left-hand approximation
%             *  1: Right-hand approximation
%             *  0: Mean approximation (default)
% (bnd)     - Boundary conditions:
%             * 0/'c'/'circulant'  (image wraps around) [default]
%             * 1/'n'/'neumann'    (mirror / boundary)
%             * 2/'d'/'dirichlet'  (zero outside FOV)
% J         - Jacobian matrix operator s.t. Jac(v) = J * v
%
% Return a matrix that allows estimating the Jacobian of a vector
% field at a given point. The partial derivative at a point is 
% estimated by:
% > dv/dx (p) ~= ( v(p + dx) - v(p) ) / dx , 
% where dx is one voxel.
% The returned matrix J can be used through 
% > Jac(v) = reshape(J * v(:), [lat_dim field_dim vec_dim]),
% where Jac(v)(x) is the Jacobian at point x.
% If vector_dim is provided, the field v is considered to be vectorized
% over both the lattice and the vectors (v = [v1; v2; v3]), and the 
% returned matrix J has dimensions 
% (prod(lattice_dim) * vector_dim) x (prod(lattice_dim)).
% Else, it is considered to be vectorized only over the lattice 
% (v = [v1 v2 v3]), and the returned matrix J has dimensions 
% (prod(lattice_dim)) x (prod(lattice_dim)).
%
% The sum of square derivatives at point v is: v'*(J'*J)*v
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    
    if nargin < 5
        bnd = 0;
        if nargin < 4
            dir = zeros(1, length(lat_dim));
            if nargin < 3
                vec_dim = 1;
                if nargin < 2
                    lat_vs = ones(1, length(lat_dim));
                end
            end
        end
    end
        
    % Multidimensional gradient functor
    CG = multi_gradient(lat_dim, dir, lat_vs, bnd);
    
    % Concatenate gradient functors
    % > CG * f allows to compute the gradient at f, where f is a scalar 
    %   field.
    CG = cell2mat(CG); 
    
    % Vector gradient functor
    % > J * v allows to compute the jacobian at v, where v is a vector 
    %   field.
    J = kron(speye(vec_dim), CG);
end

function H = hessian(lat_dim, lat_vs, vec_dim, dir, bnd)
% FORMAT H = spm_sparse('hessian', lat_dim, lat_vs, vec_dim, dir, bnd)
% lat_dim   - Dimensions of the lattice ([nx ny ...]).
% (lat_vs)  - Voxel size of the lattice (default: 1).
% (vec_dim) - Dimension of the vector space (default: 1).
% (dir)     - Approximation used in each direction ([dx dy ...]):
%             * -1: Left-hand approximation
%             *  1: Right-hand approximation
%             *  0: Mean approximation (default)
% (bnd)     - Boundary conditions:
%             * 0/'c'/'circulant'  (image wraps around) [default]
%             * 1/'n'/'neumann'    (mirror / boundary)
%             * 2/'d'/'dirichlet'  (zero outside FOV)
% H         - Hessian matrix operator s.t. Hess(v) = H * v
%
% Returns a matrix that allows estimating the Hessian of a vector 
% field at a given point. The partial derivative at a point is 
% estimated by d2v/dxdy (p) ~= ( dv/dx(p + dy) - dv/dx(p) ) / dy , 
% where dy is one voxel.
% The returned matrix J can be used through 
% > Hess(v) = reshape(H * v(:), [lat_dim n_comb vec_dim]), 
% where Hess(v)(i,:,:,j) is the Hessian at point v_i (stored sparse) 
% for component j, and n_comb = lat_dim * (lat_dim+1) / 2.
% If vector_dim is provided, the field v is considered to be vectorized
% over both the lattice and the vectors (v = [v1; v2; v3]), and the 
% returned matrix H has dimensions 
% (prod(lattice_dim) * vector_dim) x (prod(lattice_dim)).
% Else, it is considered to be vectorized only over the lattice 
% (v = [v1 v2 v3]), and the returned matrix J has dimensions 
% (prod(lattice_dim)) x (prod(lattice_dim)).
%
% The sum of square second derivatives at point v is: v'*(H'*H)*v
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    
    if nargin < 5
        bnd = 0;
        if nargin < 4
            dir = zeros(1, length(lat_dim));
            if nargin < 3
                vec_dim = 1;
                if nargin < 2
                    lat_vs = ones(1, length(lat_dim));
                end
            end
        end
    end
    im_dim = length(lat_dim);
    
    % Multidimensional gradient functor
    CG = multi_gradient(lat_dim, dir, lat_vs, bnd);
    
    % Multidimensional hessian functor
    % > CGG{index(i,j)} * f allows to compute the second derivative along 
    %   directions i and j, where f is a scalar field.
    % I don't use the symmetric approximation because the sqrt(2)
    % multiplication causes imprecision errors.
    CGG = cell(im_dim * im_dim, 1);
    index = 1;
    for i=1:im_dim
        for j=1:im_dim
            CGG{index} = CG{i} * CG{j};
            index = index + 1;
        end
    end
    clear CG
    
    % Concatenate hessian functors
    % > CGG * f allows to compute the hessian at f, where 
    %   f is a scalar field.
    CGG = -cell2mat(CGG); 
    
    % Vector hessian functor
    % > H * v allows to compute the hessian at v, where v is a vector 
    %   field.
    H = kron(speye(vec_dim), CGG);
end

function L = laplacian(lat_dim, lat_vs, dir, bnd)
% FORMAT L = spm_sparse('laplacian', lat_dim, lat_vs, dir, bnd)
% lat_dim   - Dimensions of the lattice ([nx ny ...]).
% (lat_vs)  - Voxel size of the lattice (default: 1).
% (dir)     - Approximation used in each direction ([dx dy ...]):
%             * -1: Left-hand approximation
%             *  1: Right-hand approximation
%             *  0: Mean approximation (default)
% (bnd)     - Boundary conditions:
%             * 0/'c'/'circulant'  (image wraps around) [default]
%             * 1/'n'/'neumann'    (mirror / boundary)
%             * 2/'d'/'dirichlet'  (zero outside FOV)
% L         - Laplacian matrix operator s.t. Lap(v) = L * f
%
% Returns a matrix that allows estimating the laplacian of a scalar field
% at a given point. The partial derivative at a point is 
% estimated by d2f/dx2 (p) ~= ( df/dx(p + dx) - df/dx(p) ) / dx , 
% where dx is one voxel.
% The returned matrix can be used through
% > Lapl(f) = reshape(L * f(:), [lat_dim]),
% where Lapl(f)(x) is the value of the Laplacian at point x.
%
% The sum of square laplacians at point f is: f'*(L'*L)*f
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    if nargin < 4
        bnd = 0;
        if nargin < 3
            dir = zeros(1, length(lat_dim));
            if nargin < 2
                lat_vs = ones(1, length(lat_dim));
            end
        end
    end
    im_dim = length(lat_dim);

    % Multidimensional gradient functor
    CG = multi_gradient(lat_dim, dir, lat_vs, bnd);
    
    % Multidimensional diagonal of hessian functor
    % > CGG{i} * f allows to compute the unmixed second derivative along 
    %   direction i,  where f is a scalar field.
    CGG = cell(1, im_dim);
    for i=1:im_dim
        CGG{i} = CG{i} * CG{i};
    end
    clear CG
    
    % Concatenate diagonal of hessian functors
    % > L * f allows to compute the laplacian at f, where 
    %   f is a scalar field.
    L = -cell2mat(CGG); 
end

function D = euclidian_dist(lat_dim, lat_vs)
% FORMAT D = spm_sparse('euclidian_dist', lat_dim, lat_vs)
% lat_dim   - Dimensions of the lattice ([nx ny ...]).
% (lat_vs)  - Voxel size of the lattice (default: 1).
% D         - Distance matrix operator s.t. Dist(v) = D * v
%
% Returns a matrix that allows *estimating* the euclidian norm of a  
% vector field at a given point. 
% The returned matrix D can be used through 
% > Dist(v) = reshape(D * v(:), [lat_dim vec_dim]), 
% where Dist(v)([x,j]) is the value (in mm) of component j at point x.
%
% The sum of square distance components at point v is: v'*(D'*D)*v,
% which is not the sum of square euclidian distances!
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    
    if nargin < 2
        lat_vs = ones(1, length(lat_dim));
    end
    vector_dim = length(lat_dim);
    
    D = cell(1, vector_dim);
    for i=1:vector_dim
        D{i} = speye(prod(lat_dim));
        if lat_vs(i) ~= 1
            D{i} = D{i} * lat_vs(i);
        end
    end
    D = blkdiag(D{:});
end

function A = divergence(lat_dim, lat_vs, dir, bnd)
% FORMAT A = spm_sparse('divergence', lat_dim, lat_vs, dir, bnd)
% lat_dim   - Dimensions of the lattice ([nx ny ...]).
% (lat_vs)  - Voxel size of the lattice (default: 1).
% (dir)     - Approximation used in each direction ([dx dy ...]):
%             * -1: Left-hand approximation
%             *  1: Right-hand approximation
%             *  0: Mean approximation (default)
% (bnd)     - Boundary conditions:
%             * 0/'c'/'circulant'  (image wraps around) [default]
%             * 1/'n'/'neumann'    (mirror / boundary)
%             * 2/'d'/'dirichlet'  (zero outside FOV)
% A         - Divergence matrix operator s.t. Div(v) = A * v
%
% Returns a matrix that allows estimating the divergence of a vector 
% field at a given point. 
% The returned matrix A can be used through 
% > Div(v) = reshape(A * v(:), lat_dim), 
% where Div(v)(x) is the divergence at point x.
%
% The sum of square divergences at point v is: v'*(A'*A)*v
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    if nargin < 4
        bnd = 0;
        if nargin < 3
            dir = zeros(1, length(lat_dim));
            if nargin < 2
                lat_vs = ones(1, length(lat_dim));
            end
        end
    end
    
    % Multidimensional gradient functor (column shape)
    CG = multi_gradient(lat_dim, dir, lat_vs, bnd);
    
    % Divergence functor (row shape)
    % > A * v allows to compute the divergence at v, where v is a vector 
    %   field.
    A = cell2mat(CG');
end

function S = symjac(lat_dim, lat_vs, vec_dim, dir, bnd)
% FORMAT S = spm_sparse('symjac', lat_dim, lat_vs, vec_dim, dir, bnd)
% lat_dim   - Dimensions of the lattice ([nx ny ...]).
% (lat_vs)  - Voxel size of the lattice (default: 1).
% (vec_dim) - Dimension of the vector space (default: 1).
% (dir)     - Approximation used in each direction ([dx dy ...]):
%             * -1: Left-hand approximation
%             *  1: Right-hand approximation
%             *  0: Mean approximation (default)
% (bnd)     - Boundary conditions:
%             * 0/'c'/'circulant'  (image wraps around) [default]
%             * 1/'n'/'neumann'    (mirror / boundary)
%             * 2/'d'/'dirichlet'  (zero outside FOV)
% S         - Jacobian matrix operator s.t. SymJac(v) = S * v
%
% Returns a matrix that allows estimating the symmetric part of the
% jacobian of a vector field at a given point. 
% The returned matrix S can be used through 
% > SymJac(v) = reshape(S * v, [lat_dim field_dim vec_dim]),
% where SymJac(v)(i) is the symmetric part of the Jacobian at point 
% v_i.
%
% The sum of square symmetric parts at point v is: v'*(S'*S)*v
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    
    if nargin < 5
        bnd = 0;
        if nargin < 4
            dir = zeros(1, length(lat_dim));
            if nargin < 3
                vec_dim = 1;
                if nargin < 2
                    lat_vs = ones(size(lat_dim));
                end
            end
        end
    end
    
    im_dim = length(lat_dim);
    J = jacobian(lat_dim, lat_vs, im_dim, dir, bnd);
    
    % Compute a matrix that links a the element of a Jacobian (x_i, f_j) 
    % with itself and its symmetric component (x_j, f_i).
    link_sym = spalloc(im_dim * vec_dim, im_dim * vec_dim, im_dim * vec_dim);
    for i=1:im_dim
        for j=1:vec_dim
            if i == j
                link_sym(sub2ind([im_dim vec_dim], i, j), sub2ind([im_dim vec_dim], i, j)) = 1;
            else
                link_sym(sub2ind([im_dim vec_dim], i, j), sub2ind([im_dim vec_dim], i, j)) = 0.5;
                link_sym(sub2ind([im_dim vec_dim], i, j), sub2ind([im_dim vec_dim], j, i)) = 0.5;
            end
        end
    end
    
    % Expand so that it works for any element of the vector field
    S = kron(link_sym, speye(prod(lat_dim)));
    
    % Compute the symmetric part
    S = S * J;
end

%% HELPERS

function G = simple_gradient(dim, dir, vs, bnd)
% Simple gradient functor
% > G * f allows to compute the gradient along the first direction at 
%   f, where f is a scalar field.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    if nargin < 4
        bnd = 0;
        if nargin < 3
            vs = 1;
            if nargin < 2
                dir = 0;
            end
        end
    end
    
    if dir > 0
        % Right-hand approximation
        G = spdiags(repmat([-1, 1], dim, 1), [0, 1], dim, dim);
        switch lower(bnd)
            case {0, 'c', 'circ', 'circulant'}
                G(end, 1) = 1;
            case {1, 'n', 'neumann'}
                G(end,end) = 0;
            case {2, 'd', 'dirichlet'}
                % nothing to do
        end
    elseif dir < 0
        % Left-hand approximation
        G = spdiags(repmat([-1, 1], dim, 1), [-1, 0], dim, dim);
        switch lower(bnd)
            case {0, 'c', 'circ', 'circulant'}
                G(1, end) = -1;
            case {1, 'n', 'neumann'}
                G(1,1) = 0;
            case {2, 'd', 'dirichlet'}
                % nothing to do
        end
    else
        % Mean
        Gp = spdiags(repmat([-1, 1], dim, 1), [0, 1], dim, dim);
        switch lower(bnd)
            case {0, 'c', 'circ', 'circulant'}
                Gp(end, 1) = 1;
            case {1, 'n', 'neumann'}
                Gp(end,end) = 0;
            case {2, 'd', 'dirichlet'}
                % nothing to do
        end
        Gm = spdiags(repmat([-1, 1], dim, 1), [-1, 0], dim, dim);
        switch lower(bnd)
            case {0, 'c', 'circ', 'circulant'}
                Gm(1, end) = -1;
            case {1, 'n', 'neumann'}
                Gm(1,1) = 0;
            case {2, 'd', 'dirichlet'}
                % nothing to do
        end
        G = 0.5 * (Gp + Gm);
    end
    if vs ~= 1
        G = G / vs;
    end
end

function CG = multi_gradient(dim, dir, vs, bnd)
% Multidimensional gradient functor
    % > CG{i} * f allows to compute the gradient along direction i at f,  
    %   where f is a scalar field.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    if nargin < 4
        bnd = zeros(1, length(lat_dim));
        if nargin < 3
            vs = ones(1, length(lat_dim));
            if nargin < 2
                dir = zeros(1, length(lat_dim));
            end
        end
    end
    im_dim = length(dim);
    if numel(bnd) == 1
        bnd = repmat(bnd, [1 im_dim]);
    end
    if numel(vs) == 1
        vs = repmat(vs, [1 im_dim]);
    end
    if numel(dir) == 1
        dir = repmat(dir, [1 im_dim]);
    end
    
    % Create gradient operators for each dimension
    G = cell(im_dim, 1);
    for i=1:im_dim
        G{i} = simple_gradient(dim(i), dir(i), vs(i), bnd(i));
    end
    
    % Extend matrix
    CG = cell(im_dim, 1); 
    for i=1:im_dim
        GorI = cell(im_dim, 1);
        for j=1:im_dim
            if i ~= j
                GorI{j} = speye(dim(j));
            else
                GorI{j} = G{i};
            end
        end
        CG{i} = speye(1);
        for j = 1:im_dim
            CG{i} = kron(GorI{j}, CG{i});
        end
    end
end