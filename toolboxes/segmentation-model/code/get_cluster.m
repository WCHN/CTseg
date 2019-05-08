function cluster = get_cluster(obs,bf,dm,GaussPrior,miss,in,varargin)
% FORMAT cluster = get_cluster(obs,bf,dm,GaussPrior,miss,in,['sort_pars',true/false])
% obs        - Observed image [Nx Ny Nz P] 
% bf         - Bias field [Nx Ny Nz P] (or 1 if no bias field)
% dm         - Image dimensions
% GaussPrior - GMM prior parameters {MU,b,V,n}
% miss       - Structure handling missing data (fields: C, L, nL)
% in         - Initial GMM posterior parameters {{MU,b},{V,n}}
% sort_pars  - If true, sort cluster by ascending magnitude.
% cluster    - Final GMM posterior parameters {{MU,b},{V,n}}
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Parse inputs
%--------------------------------------------------------------------------
p         = inputParser;
p.addParameter('sort_pars', true,  @islogical);
p.parse(varargin{:});
sort_pars = p.Results.sort_pars;

% Get GMM parameters (ML)
%--------------------------------------------------------------------------

if numel(bf) == 1
    bf = ones(size(obs),'single');
end

% Compute suffstats
SS0  = 0;
SS1  = 0;
SS2  = 0;
SS2b = 0;
for z=1:dm(3)            
    ix            = ix_slice(z,prod(dm(1:2)));   
    slice.obs     = double(bf(ix,:)).*double(obs(ix,:));
    slice.code    = miss.C(ix);
    slice.bin_var = (double(bf(ix,:)).^2)./12;

    if ~iscell(in)                    
        slice.Z   = double(in(ix,:));
    else        
        slice.template = double(in{1}(ix,:));
        if ~isempty(in{3})
            ix_l         = in{3}{1}(ix);
            CM           = double(in{3}{2});            
            slice.labels = CM(ix_l,:);
        else
            slice.labels = ones(1,size(in{1},2));
        end

        logPL = log(slice.labels);
        logPL = logPL(:,in{4}.lkp);

        logPI = log(spm_matcomp('softmax',slice.template,in{2}) + eps);
        logPI = logPI(:,in{4}.lkp);

        logPI(isnan(logPI)) = 0; 

        slice.Z = bsxfun(@plus, logPI, logPL);
        slice.Z = bsxfun(@plus, slice.Z,     log(in{4}.mg));

        slice.Z = bsxfun(@minus, slice.Z, max(slice.Z, [], 2));
        slice.Z = exp(slice.Z);
        slice.Z = bsxfun(@rdivide, slice.Z, sum(slice.Z, 2));
    end
    
    [dSS0,dSS1,dSS2] = spm_gmm_lib('SuffStat',slice.obs,slice.Z,1);
    
    SS2b = SS2b + spm_gmm_lib('SuffStat', 'bin', slice.bin_var, slice.Z, 1, {slice.code,miss.L});  
    
    SS0 = SS0 + dSS0;
    SS1 = SS1 + dSS1;
    SS2 = SS2 + dSS2;
end
SS2 = SS2 + SS2b;

% Compute ML GMM parameters
[MU,A] = spm_gmm_lib('UpdateClusters',SS0,SS1,SS2);

% Get GMM parameters (VB)
%--------------------------------------------------------------------------

lSS0 = cell(1,miss.nL); lSS0(:) = {0};
lSS1 = cell(1,miss.nL); lSS1(:) = {0};
lSS2 = cell(1,miss.nL); lSS2(:) = {0};
SS2b = 0;

for z=1:dm(3)        
    ix            = ix_slice(z,prod(dm(1:2)));   
    slice.obs     = double(bf(ix,:)).*double(obs(ix,:));
    slice.code    = miss.C(ix);
    slice.bin_var = (double(bf(ix,:)).^2)./12;
    
    if ~iscell(in)  
        slice.Z    = double(in(ix,:));        
    else        
        slice.template = double(in{1}(ix,:));
        if ~isempty(in{3})
            ix_l         = in{3}{1}(ix);
            CM           = double(in{3}{2});            
            slice.labels = CM(ix_l,:);
        else
            slice.labels = ones(1,size(in{1},2));
        end

        logPL = log(slice.labels);
        logPL = logPL(:,in{4}.lkp);

        logPI = log(spm_matcomp('softmax',slice.template,in{2}) + eps);
        logPI = logPI(:,in{4}.lkp);

        logPI(isnan(logPI)) = 0; 

        slice.Z = bsxfun(@plus, logPI, logPL);
        slice.Z = bsxfun(@plus, slice.Z,     log(in{4}.mg));

        slice.Z = bsxfun(@minus, slice.Z, max(slice.Z, [], 2));
        slice.Z = exp(slice.Z);
        slice.Z = bsxfun(@rdivide, slice.Z, sum(slice.Z, 2));
    end
    
    [dlSS0,dlSS1,dlSS2] = spm_gmm_lib('SuffStat', 'base', slice.obs, slice.Z, 1, {slice.code,miss.L});   
    
    SS2b = SS2b + spm_gmm_lib('SuffStat', 'bin', slice.bin_var, slice.Z, 1, {slice.code,miss.L});   
    
    for l=1:miss.nL
        if isempty(dlSS0{l}), continue; end

        lSS0{l} = lSS0{l} + dlSS0{l};
        lSS1{l} = lSS1{l} + dlSS1{l};
        lSS2{l} = lSS2{l} + dlSS2{l};
    end  
end

for i=1:32
    % Compute suffstats
    [SS0,SS1,SS2] = spm_gmm_lib('SuffStat', 'infer', lSS0, lSS1, lSS2, {MU,A}, miss.L);
    SS2           = SS2 + SS2b;
    
    % Compute VB GMM parameters
    [MU,A,b,V,n] = spm_gmm_lib('UpdateClusters',SS0, SS1, SS2, GaussPrior);
end

if sort_pars
    % Sort GMM parameters    
    [~,ix] = sort(vecnorm(MU,2,1));
    
    MU = MU(:,ix);
    b  = b(ix);
    V  = V(:,:,ix);
    n  = n(ix);
end
    
cluster = {{MU,b},{V,n}};
%==========================================================================     