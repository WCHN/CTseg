function dat = init_dat(method,Nii,mat,dm,window,gap,gapunit)
% Initialise projection matrices for super-resolution
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 5, window  = []; end
if nargin < 6, gap     = 0; end
if nargin < 7, gapunit = '%'; end

% Get rigid basis
is3d = dm(3) > 1;
B    = get_rigid_basis(is3d);

% Slice profile
window = get_window(window,Nii);

% Slice gap
gap = get_slice_gap(gap,Nii,gapunit);

C   = numel(Nii);
dat = struct('mat',[],'dm',[],'N',[],'A',[],'method',[]);
for c=1:C % Loop over channels
    
    matc = mat(:,:,min(c,size(mat,3)));
    
    if iscell(Nii(c))
        dat(c) = init_A(Nii{c},matc,dm,window{c},gap{c},B,method); 
    else
        dat(c) = init_A(Nii(c),matc,dm,window{c},gap{c},B,method);    
    end
end    
%==========================================================================

%==========================================================================
function dat = init_A(Nii,Mmu,dmmu,window,gap,B,method)
% Initialise projection matrices (stored in dat struct)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

N          = numel(Nii); % Number of LR images
Nq         = size(B,3);  % Number of affine parameters
dat.mat    = Mmu;
dat.dm     = dmmu;   
dat.N      = N;
dat.method = method;
vs         = sqrt(sum(Mmu(1:3,1:3).^2));
for n=1:N % Loop over LR images
    Mf  = Nii(n).mat;    
    dmf = Nii(n).dat.dim;
    dmf = [dmf 1];
    
    dat.A(n).mat = Mf;    
    dat.A(n).dm  = dmf;
    dat.A(n).win = window{n};
    dat.A(n).gap = gap{n};
    dat.A(n).q   = zeros([Nq 1]); 
    
    R          = spm_dexpm(dat.A(n).q,B); % Rigid matrix
    dat.A(n).R = R;
    M          = Mmu\R*Mf;
    
%     R          = (M(1:3,1:3)/diag(sqrt(sum(M(1:3,1:3).^2))))';
%     dat.A(n).S = blur_fun(dm,R,sqrt(sum(M(1:3,1:3).^2)));
%     dat.A(n).S = blur_function(dm,M);

    % Include slice-gap
    M = model_slice_gap(M,gap{n},vs);
    
    dat.A(n).J = single(reshape(M, [1 1 1 3 3]));
end
%==========================================================================

%==========================================================================
function gap = get_slice_gap(gap,Nii_x,gapunit)
% Construct slice-gap
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C = numel(Nii_x);

if isscalar(gap)
    % Find thick-slice direction automatically
    G = cell(1,C);
    for c=1:C
        N    = numel(Nii_x{c});
        G{c} = cell(1,N);
        for n=1:N
            G{c}{n} = zeros(1,3);
            
            thickslice_dir          = get_thickslice_dir(Nii_x{c}(n));
            G{c}{n}(thickslice_dir) = gap;            
        end
    end    
    gap = G;
else
    % x-, y-, z-directions given
    gap = padarray(gap, [0 max(0,C-numel(gap))], 'replicate', 'post');
    for c=1:C
        if ~iscell(gap{c})
            gap{c} = {gap{c}};
        end
        gap{c} = padarray(gap{c}, [0 max(0,numel(Nii_x{c})-numel(gap{c}))], 'replicate', 'post');
        for i=1:numel(gap{c})
            if isempty(gap{c}{i})
                gap{c}{i} = 0;
            end
            gap{c}{i} = padarray(gap{c}{i}, [0 max(0,3-numel(gap{c}{i}))], 'replicate', 'post');
        end
    end
end

if strcmp(gapunit,'%')
    % Convert from percentage to mm
    for c=1:C
        N    = numel(Nii_x{c});
        for n=1:N       
            mat = Nii_x{c}(n).mat;
            vx  = sqrt(sum(mat(1:3,1:3).^2));   
            
            gap{c}{n} = vx.*gap{c}{n};
        end
    end    
end
%==========================================================================

%==========================================================================
function window = get_window(window,Nii_x)
% Construct slice-profile
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C = numel(Nii_x);

if isempty(window)
    % Find thick-slice direction automatically and set defaults:
    % In-plane: Gaussian, Through-plane: Rectangle
    W = cell(1,C);
    for c=1:C
        N    = numel(Nii_x{c});
        W{c} = cell(1,N);
        for n=1:N
            W{c}{n} = ones(1,3); % All Gaussian
            
%             thickslice_dir          = get_thickslice_dir(Nii_x{c}(n));
%             W{c}{n}(thickslice_dir) = 2; % Rectangle
        end
    end    
    window = W;    
else
    % x-, y-, z-directions given
    if ~iscell(window)
        window = {window};
    end
    window = padarray(window, [0 max(0,C-numel(window))], 'replicate', 'post');
    for c=1:C
        if ~iscell(window{c})
            window{c} = {window{c}};
        end
        window{c} = padarray(window{c}, [0 max(0,numel(Nii_x{c})-numel(window{c}))], 'replicate', 'post');
        for i=1:numel(window{c})
            if isempty(window{c}{i})
                window{c}{i} = 2;
            end
            window{c}{i} = padarray(window{c}{i}, [0 max(0,3-numel(window{c}{i}))], 'replicate', 'post');
        end
    end
end
%==========================================================================

%==========================================================================
function f = blur_fun(dm,mat,vx)
if nargin<1, dm = [64 64]; end
if nargin<2, mat = eye(numel(dm)); end
if nargin<3, vx = ones(1,numel(dm)); end

if any(size(mat)~=numel(dm)) || numel(vx)~=numel(dm), error('Incompatible dimensions.'); end

% Grid in frequency space
r        = cell(1,numel(dm));
for i=1:numel(dm) 
    r{i} = single([0:ceil(dm(i)/2-1) -floor(dm(i)/2):-1]'*pi/dm(i)); 
end
X        = cell(1,numel(dm));
[X{:}]   = ndgrid(r{:});
clear r

% Transform
Y            = cell(size(X));
for i=1:numel(dm)
    Y{i}     = single(0);
    for j=1:numel(dm) 
        Y{i} = Y{i} + mat(i,j)*X{j}; 
    end
end
clear X

% Window function
f     = single(0);
for i=1:numel(dm) 
    f = f + Y{i}.^2; 
end    
f     = ((cos(min(f,pi^2/4)*4/pi) + 1)/2);

% Incorporate voxel size
for i=1:numel(dm)
    tmp                 = sin((vx(i))*Y{i})./(Y{i}.*cos(Y{i}/pi^(1/2)));
    tmp(~isfinite(tmp)) = vx(i);
    f                   = f.*tmp;
end
%==========================================================================

%==========================================================================
function thickslice_dir = get_thickslice_dir(Nii)
mat                = Nii.mat;
vx                 = sqrt(sum(mat(1:3,1:3).^2));   
[~,thickslice_dir] = max(vx);
%==========================================================================            