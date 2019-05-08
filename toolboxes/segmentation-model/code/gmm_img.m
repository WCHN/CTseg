function varargout = gmm_img(varargin)
%__________________________________________________________________________
% GMM functions for images (2D and 3D).
%
% FORMAT [cluster,prop,lb,mom] = gmm_img('update_gmm_pars',X,bf,cluster,prop,Template,lb,dm,miss,part,varargin)
% FORMAT [lb,mom]              = gmm_img('img_lb_and_mom',obs,bf,BinWidth,template,labels,prop,mean,prec,miss,part,dm,mrf,varargin)
% FORMAT [slice,ix]            = gmm_img('get_slice_data',z,dm,obs,bf,template,code,labels,BinWidth)
% FORMAT [Z,lb,BX]             = gmm_img('slice_resp_and_lb',slice,mean,prec,prop,part,miss,const,lb,ix,varargin)
% FORMAT mom                   = gmm_img('slice_mom',mom,Z,slice,miss,BX)
% FORMAT [lb,mom]              = gmm_img('init_lb_and_mom',miss)
%
% FORMAT help gmm_img>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help gmm_img
    error('Not enough argument. Type ''help gmm_img'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'update'
        [varargout{1:nargout}] = update_gmm_pars(varargin{:});           
    case 'img_lb_and_mom'
        [varargout{1:nargout}] = img_lb_and_mom(varargin{:});  
    case 'slice_resp_and_lb'
        [varargout{1:nargout}] = slice_resp_and_lb(varargin{:});            
    case 'slice_mom'
        [varargout{1:nargout}] = slice_mom(varargin{:});            
    case 'getslice'
        [varargout{1:nargout}] = get_slice_data(varargin{:});            
    case 'init_lb_and_mom'
        [varargout{1:nargout}] = init_lb_and_mom(varargin{:});            
    otherwise
        help gmm_img
        error('Unknown function %s. Type ''help gmm_img'' for help.', id)
end
%==========================================================================

%==========================================================================
function [cluster,prop,lb,mom,mrf,part] = update_gmm_pars(X,bf,cluster,prop,Template,lb,dm,miss,part,mrf,ix_tiny,do_mg,varargin)
% FORMAT [cluster,prop,lb,mom,mrf] = gmm_img('update_gmm_pars',X,bf,cluster,prop,Template,lb,dm,miss,part,mrf,varargin)
%
% Short description
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% -------------------------------------------------------------------------
% Parse inputs
p = inputParser;
p.FunctionName = 'spm_gmm_loop';
p.addParameter('Resp',           [],    @isnumeric);
p.addParameter('GaussPrior',     {},    @iscell);
p.addParameter('PropPrior',      0,     @isnumeric);
p.addParameter('IterMax',        1024,  @(X) isscalar(X) && isnumeric(X));
p.addParameter('Tolerance',      1e-4,  @(X) isscalar(X) && isnumeric(X));
p.addParameter('SubIterMax',     32,   @(X) isscalar(X) && isnumeric(X));
p.addParameter('SubTolerance',   1e-4,  @(X) isscalar(X) && isnumeric(X));
p.addParameter('BinWidth',       0,     @isnumeric);
p.addParameter('Verbose',        0,     @(X) isscalar(X) && (isnumeric(X) || islogical(X)));
p.addParameter('Labels',         {},    @iscell);
p.addParameter('Constrained',    [],    @(X) isstruct(X) || isnumeric(X));
p.parse(varargin{:});
BinWidth     = p.Results.BinWidth;
GaussPrior   = p.Results.GaussPrior;
IterMax      = p.Results.IterMax;
Tolerance    = p.Results.Tolerance;
SubIterMax   = p.Results.SubIterMax;
SubTolerance = p.Results.SubTolerance;
Verbose      = p.Results.Verbose;
Labels       = p.Results.Labels;
Constrained  = p.Results.Constrained;

% -------------------------------------------------------------------------
% Unfold inputs
b   = 0;  % Mean degrees of freedom (posterior)
n   = 0;  % Precision degrees of freedom (posterior)
V   = []; % Scale matrix (posterior)
b0  = 0;  % Mean degrees of freedom (prior)
n0  = 0;  % Precision degrees of freedom (prior)
MU0 = []; % Mean (prior)
V0  = []; % Scale matrix (prior)

lkp = part.lkp;

% Get means + variances (ML) or posteriors over these (VB)
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

mean = {MU,b};
prec = {V,n};

% Get priors (VB)
if numel(GaussPrior) >= 1
    MU0 = GaussPrior{1};
    if numel(GaussPrior) >= 2
        b0 = GaussPrior{2};
        if numel(GaussPrior) >= 3
            V0 = GaussPrior{3};
            if numel(GaussPrior) >= 4
                n0 = GaussPrior{4};
            end
        end
    end
end

% Number of classes
K = size(MU,2); 

% Compute sufficient statistics
[dlb,mom,mrf] = img_lb_and_mom(X,bf,BinWidth,Template,Labels,prop,mean,prec,miss,part,dm,mrf,ix_tiny);
                   
% -------------------------------------------------------------------------
% EM loop
for em=1:IterMax    
                
    % ---------------------------------------------------------------------
    % sub-EM algorithm to update Mean/Precision with missing data
    % . Responsibilities (E[z]) are kept fixed
    % . Missing values (E[z*h], E[z*hh']) are updated
    % . Cluster parameters (MU,b,A/V,n) are updated

    LB         = NaN(1,SubIterMax);
    [lbMU,lbA] = spm_gmm_lib('KL', 'GaussWishart', {MU,b}, prec, {MU0,b0}, {V0,n0});
    LB(1)      = lbMU + lbA + dlb.X;

    for i=1:SubIterMax

        oMU  = MU;
        ob   = b;
        oV   = V;
        on   = n;
        oA   = A;
        olbX = dlb.X;

        % -------------------------------------------------------------
        % Infer missing suffstat
        % sum{E[z]}, sum{E[z*x]}, sum{E[z*xx']}
        [SS0,SS1,SS2] = spm_gmm_lib('SuffStat', 'infer', mom.SS0, mom.SS1, mom.SS2, {MU,A}, miss.L);
        SS2           = SS2 + mom.SS2b;

        % -------------------------------------------------------------
        % Update GMM
        [MU,A1,b,V1,n1] = spm_gmm_lib('UpdateClusters', ...
                                      SS0, SS1, SS2, {MU0,b0,V0,n0});
        for k=1:size(MU,2)
            [~,cholp] = chol(A1(:,:,k));
            if cholp == 0
                A(:,:,k) = A1(:,:,k);
                if sum(n1) > 0
                    V(:,:,k) = V1(:,:,k);
                    n(k)     = n1(k);
                end
            end
        end
        mean             = {MU,b};
        if ~sum(n), prec = {A};
        else,       prec = {V,n};   end

        % -------------------------------------------------------------
        % Marginal / Objective function
        [lbMU,lbA] = spm_gmm_lib('KL', 'GaussWishart', {MU,b}, prec, {MU0,b0}, {V0,n0});
        dlb.X      = spm_gmm_lib('MarginalSum', mom.SS0, mom.SS1, mom.SS2, mean, prec, miss.L, mom.SS2b);
        LB(i + 1)  = lbMU + lbA + dlb.X;            

        subgain = (LB(i + 1) - LB(i))/(max(LB(1:i + 1), [], 'omitnan') - min(LB(1:i + 1), [], 'omitnan'));
        if subgain < 0
            MU      = oMU;
            A       = oA;
            V       = oV;
            b       = ob;
            n       = on;
            dlb.X   = olbX;
            subgain = 0;

            mean             = {MU,b};
            if ~sum(n), prec = {A};
            else,       prec = {V,n};   end

            [lbMU,lbA] = spm_gmm_lib('KL', 'GaussWishart', {MU,b}, prec, {MU0,b0}, {V0,n0});
        end
        if Verbose >= 5
            switch sign(subgain)
                case 1,     incr = '(+)';
                case -1,    incr = '(-)';
                case 0,     incr = '(=)';
                otherwise,  incr = '';
            end
            fprintf('%-5s | %4d | lb = %-10.6f | gain = %-10.4g | %3s\n', 'sub', i, LB(i + 1), subgain, incr);
        end

        if subgain < SubTolerance && i > 1
            break
        end
    end      
    
    if do_mg        
        for k=1:numel(lkp)
            tmp   = SS0(lkp == lkp(k));
            mg(k) = (SS0(k) + eps)/sum(tmp + eps);
        end
        part.mg = mg;
    end
    
    % Compute sufficient statistics
    [dlb,mom,mrf] = img_lb_and_mom(X,bf,BinWidth,Template,Labels,prop,mean,prec,miss,part,dm,mrf,ix_tiny);
                      
    % ---------------------------------------------------------------------
    % Compute lower bound
    lb.Z(end+1)  = dlb.Z;       
    lb.MU(end+1) = lbMU;
    lb.A(end+1)  = lbA;
    lb.X(end+1)  = dlb.X;
    if numel(part.mg) > numel(prop)
        lb.mg(end+1) = dlb.mg;    
    end
    if ~isempty(Labels)
        lb.lab(end+1) = dlb.lab;
    end
    if mrf.do   
        lb.ZN(end+1) = dlb.ZN;
    end
    
    if any(Constrained.ElnDetV ~= 0)
        % Corrent lower bound for when a prior is placed on V0
        
        ElnDetV = Constrained.ElnDetV;
        
        for k=1:K
            lb.A(end) = lb.A(end) + 0.5*n0(k)*spm_matcomp('LogDet', V0(:,:,k)) ...
                                  - 0.5*n0(k)*ElnDetV(k);         
        end   
    end
    
    % Append lower bound
    [lb,gain] = check_convergence('gmm', lb, em, Verbose);

    % ---------------------------------------------------------------------
    % Check convergence    
    if gain < Tolerance && em > 1
        % Finished
        break;
    end            
end % <-- End loop em

% -------------------------------------------------------------------------
% Format output
cluster = {{MU,b},{V,n}};
%==========================================================================

%========================================================================== 
function [lb,mom,mrf] = img_lb_and_mom(obs,bf,BinWidth,template,labels,prop,mean,prec,miss,part,dm,mrf,ix_tiny,varargin)
% FORMAT [lb,mom,mrf] = gmm_img('img_lb_and_mom',obs,bf,BinWidth,template,labels,prop,mean,prec,miss,part,dm,mrf,varargin)
%
% Short description
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

K = numel(prop);

% Inits
[lb,mom] = init_lb_and_mom(miss);

%
const = spm_gmm_lib('Const', mean, prec, miss.L);

% Compute neighborhood part of responsibilities (from previous responsibilities)
lnPzN = gmm_mrf('apply',mrf);

for z=1:dm(3) % Loop over slices    
    
    % Get slice data
    [slice,ix] = get_slice_data(z,dm,obs,bf,template,miss.C,labels,BinWidth);

    if mrf.do        
        lnPzNz = double(lnPzN(ix,:));
    else
        lnPzNz = lnPzN;
    end

    % Compute responsibilities and lb
    [Z,lb,BX] = slice_resp_and_lb(slice,mean,prec,prop,part,miss,const,lnPzNz,ix_tiny,lb,ix,varargin{:});
       
    if mrf.do        
        mrf.oZ(:,:,z,:) = uint8((2^8)*reshape(cluster2template(Z,part),[dm(1:2) 1 K]));
    end
    
    % Compute sufficient statistics 
    mom = slice_mom(mom,Z,slice,miss,BX);
    
end % <-- End loop slices

lb.X = spm_gmm_lib('MarginalSum', mom.SS0, mom.SS1, mom.SS2, mean, prec, miss.L, mom.SS2b);
  
% Compute MRF part of lower bound (always requires most updated responsibilities)
lb.ZN = gmm_mrf('lowerbound',mrf); 
%========================================================================== 

%==========================================================================
function [Z,lb,BX] = slice_resp_and_lb(slice,mean,prec,prop,part,miss,const,lnPzN,ix_tiny,lb,ix_vox,varargin)
% FORMAT [Z,lb,BX] = gmm_img('slice_resp_and_lb',slice,mean,prec,prop,part,miss,const,lnPzN,lb,ix_vox,ix_tiny,varargin)
%
% Short description
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Parameters
MU     = mean{1};
lkp    = part.lkp;
mg     = part.mg;
K      = numel(prop);
lntiny = log(eps);

% Multiply bias field with observed data
BX = slice.bf.*slice.obs;

% Mask out where there are no observations and compute ln|bf|    
lnDetbf   = log(prod(slice.bf,2));    

%
if numel(varargin) == 1 && strcmpi(varargin{1}{1},'bf')                
    obf                   = double(varargin{1}{2}(ix_vox,:));        
    oBX                   = bsxfun(@times,slice.obs,obf);    
    obf(isnan(slice.obs)) = 1;
    olblnDetbf            = log(prod(obf,2));         
    olnpX                 = spm_gmm_lib('Marginal', oBX, [{MU} prec], const, {slice.code,miss.L}, (obf.^2)/12);
    olnpX(:,ix_tiny)      = lntiny;
elseif numel(varargin) == 1 && strcmpi(varargin{1}{1},'Template')        
    oTemplate             = double(varargin{1}{2}(ix_vox,:));
    ologPI                = log(spm_matcomp('softmax',oTemplate,prop) + eps);
elseif numel(varargin) == 1 && strcmpi(varargin{1}{1},'prop')            
    ologPI                = log(spm_matcomp('softmax',slice.template,varargin{1}{2}) + eps);        
end

% Compute ln(p(x))
lnpX = spm_gmm_lib('Marginal', BX, [{MU} prec], const, {slice.code,miss.L}, slice.bin_var);

% Compute ln(p(z))
if isempty(slice.template)   
    lnPI         = reshape(prop,[],numel(prop));
    slice.labels = ones(1,size(lnPI,2));
else
    lnPI = log(spm_matcomp('softmax',slice.template,prop) + eps);
    lnPI = lnPI(:,lkp);
end

% Get log of label part
lnPl = log(slice.labels);
lnPl = lnPl(:,lkp);

% MRF part
lnPzN = lnPzN(:,lkp);

% Force responsibilties to zero for ix_tiny classes
lnpX(:,ix_tiny) = lntiny;

% Compute responsibilities
%----------------------------------------------------------------------
if numel(varargin) == 0
    % Default
    Z = spm_gmm_lib('Responsibility', lnpX,  lnPI,         lnDetbf,   lnPl,log(mg),lnPzN);         
elseif strcmpi(varargin{1}{1},'bf')            
    % Using old bf parameters
    Z = spm_gmm_lib('Responsibility', olnpX,lnPI,         olblnDetbf,lnPl,log(mg),lnPzN); 
elseif strcmpi(varargin{1}{1},'Template') || strcmpi(varargin{1}{1},'prop')  
    % Using old template parameters
    Z = spm_gmm_lib('Responsibility', lnpX,  ologPI(:,lkp),lnDetbf,   lnPl,log(mg),lnPzN);         
end

% Set non-observed to zero
%----------------------------------------------------------------------
for k=1:numel(lkp)
    Z(slice.code == 0,k) = 0;
end

if nargout > 1 && nargin >= 10
    % Compute lower bound
    %----------------------------------------------------------------------
    lbZ  = spm_gmm_lib('KL', 'Categorical', Z, 1, lnPI);
    lbmg = sum(sum(Z .* log(mg), 2) .* 1, 1);
    if numel(lnPl) ~= K
        lbLab  = sum(sum(Z .* lnPl, 2) .* 1, 1);
        lb.lab = lb.lab + lbLab;
    end  

    lb.Z   = lb.Z   + lbZ;
    lb.mg  = lb.mg  + lbmg;
else
    lb     = struct;
end
%==========================================================================

%==========================================================================
function mom = slice_mom(mom,Z,slice,miss,BX)
% FORMAT mom = gmm_img('slice_mom',mom,Z,slice,miss,BX)
%
% Short description
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Compute sufficient statistics (bin uncertainty part)
dSS2b = spm_gmm_lib('SuffStat', 'bin', slice.bin_var, Z, 1, {slice.code,miss.L});                

if nargin < 4
    % Multiply bias field with observed data
    BX = slice.bf.*slice.obs;
end

% Compute fast sufficient statistics for each configuration of missing data
[dSS0,dSS1,dSS2] = spm_gmm_lib('SuffStat', 'base', BX, Z, 1, {slice.code,miss.L});

% Add up suffstats
mom.SS2b = mom.SS2b + dSS2b;    
for l=1:miss.nL
    if isempty(dSS0{l}), continue; end

    mom.SS0{l} = mom.SS0{l} + dSS0{l};
    mom.SS1{l} = mom.SS1{l} + dSS1{l};
    mom.SS2{l} = mom.SS2{l} + dSS2{l};
end  
%==========================================================================

%========================================================================== 
function [slice,ix] = get_slice_data(z,dm,obs,bf,template,code,labels,BinWidth)
% FORMAT [slice,ix] = gmm_img('get_slice_data',z,dm,obs,bf,template,code,labels,BinWidth)
%
% Short description
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 8, BinWidth = []; end

K  = size(template,2);
ix = ix_slice(z,prod(dm(1:2)));

% Store slice data in a struct
%--------------------------------------------------------------------------
             
% Image
slice.obs = double(obs(ix,:));

% Bias field
if numel(bf) == 1
    slice.bf = bsxfun(@times,BinWidth,ones(size(slice.obs)));
else
    slice.bf = double(bf(ix,:));    
end        
slice.bf(isnan(slice.obs)) = 1;

% Template
if isempty(template)
    slice.template = template;
else
    slice.template = double(template(ix,:));
end

% Code
slice.code = code(ix);

% Bin variance
if numel(bf) == 1
    slice.bin_var = (double(BinWidth).^2)./12;
else
    slice.bin_var = (slice.bf.^2)./12;
end

% Labels
if ~isempty(labels)
    ix_l         = labels{1}(ix);
    CM           = double(labels{2});            
    slice.labels = CM(ix_l,:);
else
    slice.labels = ones(1,K);
end
%========================================================================== 

%==========================================================================
function [lb,mom] = init_lb_and_mom(miss)
% FORMAT [lb,mom] = gmm_img('init_lb_and_mom',miss)
%
% Short description
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

lb.Z   = 0;
lb.mg  = 0;
lb.lab = 0;

mom      = struct;
mom.L    = miss.L;

SS0  = cell(1,miss.nL); SS0(:) = {0};
SS1  = cell(1,miss.nL); SS1(:) = {0};
SS2  = cell(1,miss.nL); SS2(:) = {0};
SS2b = 0;

mom.SS2b = SS2b;    
mom.SS0  = SS0;
mom.SS1  = SS1;
mom.SS2  = SS2;  
%==========================================================================