%==========================================================================
function F = blur_function(dm,mat,win)
% FORMAT F = blur_function(dm,mat,win)
%
% dm  - High-resolution image dimensions                    [64 64]
% mat - Ratio of voxel-to-world matrices mat_hr\mat_lr      [eye]
%       > Defines the shape of the convolution kernel
% win - A window function, or a cell of window functions    ['tukey']
%
% Returns the Fourier representation of the convolution kernel.

use_scale = false;

function f = rect(x)
    f = single(abs(x) < 0.5);
    f(x==0.5) = 0;
end
function f = gauss(x,sigma)
    if nargin < 2
        sigma = 1/6;
    end
    f = normpdf(single(x),0,sigma);
end
function f = tukey(x,alpha)
    if nargin < 2
        alpha = 0.1;
    end
    if alpha == 0
        f = rect(x);
        return
    end
    f = 0.5*(1 + cos(pi*(2*abs(single(x))+alpha-1)/alpha));
    f(abs(x) >= 0.5) = 0;
    f(abs(x) <= 0.5*(1-alpha)) = 1;
end

% -------------------------------------------------------------------------
% Default parameters
if nargin < 1, dm  = [64 64];           end
if nargin < 2, mat = eye(numel(dm));    end
if nargin < 3, win = @rect;             end

% if any(size(mat)~=numel(dm))
%     error('Incompatible dimensions.');
% end

ndm = numel(dm);

% -------------------------------------------------------------------------
% Replicate window if needed to match the number of dimensions
if ~iscell(win), win = {win}; end
win = padarray(win(:)', [0 max(0,ndm-numel(win))], 'replicate', 'post');

for i=1:numel(win)
    if ischar(win{i})
        win{i} = str2func(win{i});
    end
end

% -------------------------------------------------------------------------
% Compute rotation matrix and voxel size
M  = mat(1:ndm,1:ndm);
vs = sqrt(sum(M.^2));
vs = vs;
R  = (M/diag(vs));

% -------------------------------------------------------------------------
% Compute scaling factor that fillsthe FOV with the shape
% (to get additional precision)
if use_scale
    switch ndm
        case 2
            bbox = [-0.5 -0.5  0.5  0.5; ...
                    -0.5  0.5 -0.5  0.5];
        case 3
            bbox = [-0.5 -0.5 -0.5 -0.5  0.5  0.5  0.5  0.5; ...
                    -0.5 -0.5  0.5  0.5 -0.5 -0.5  0.5  0.5; ...
                    -0.5  0.5 -0.5  0.5 -0.5  0.5 -0.5  0.5];
    end
    bbox  = M*bbox;
    bbox  = minmax(bbox);
    scale = bbox(:,2) - bbox(:,1);
    scale = scale ./ (dm(:)-1);
    scale = max(scale);
else
    scale = 1;
end

% -------------------------------------------------------------------------
% Create grid
% > First, r contains the coordinates in each dimension. Zero is top-left.
r = cell(1,ndm);
for i=1:ndm 
    r{i} = single([0:ceil(dm(i)/2-1) -floor(dm(i)/2):-1]'); 
end
% > Then, X builds the coordinate system of the whole lattice from the 
%   separable coordinates.
X      = cell(1,ndm);
[X{:}] = ndgrid(r{:});
% clear r

% -------------------------------------------------------------------------
% Transform
% > Apply the rotation part of the slice selection matrix to the coordinate
%   system. Now, coordinates are aligned with those of the low-resolution
%   image, but on a high-resolution lattice.
R = (scale*eye(ndm))*diag(1./vs)*R';
Y = cell(size(X));
for i=1:ndm
    Y{i} = single(0);
    for j=1:ndm 
        Y{i} = Y{i} + R(i,j)*X{j}; 
    end
end
% clear X

% -------------------------------------------------------------------------
% Compute window function on rotated lattice
f = ones(dm, 'single');
for i=1:ndm
    f = f .* win{i}(Y{i});
end
clear Y

% -------------------------------------------------------------------------
% Fourier transform
if use_scale
    for i=1:ndm
        F = zeros(size(f), 'single');
        for z=1:dm(i)
            fz = select_slices(f, i, z);
            F  = F + scale .* bsxfun(@times, fz, exp(-2*1i*pi*X{i}*r{i}(z)*scale/dm(i)));
        end
        f = F;
    end
else
    F = prod(dm) * ifftn(f, 'symmetric');
end

end
%==========================================================================