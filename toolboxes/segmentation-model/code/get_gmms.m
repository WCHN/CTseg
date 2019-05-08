function [dat,model] = get_gmms(obs,model,dat,opt)
S0      = numel(obs);
K       = opt.template.K;
init_ix = opt.gmm.hist.init_ix;
tiny    = get_par('tiny');
        
% Set posteriors
%--------------------------------------------------------------------------
parfor s=1:S0    
% for s=1:S0, fprintf('obs! for s=1:S0\n')
    
    modality    = dat{s}.modality{1}.name;   
    dm          = obs_info(dat{s}); 
    lkp         = get_par('lkp',modality,opt);
    [~,ix_tiny] = get_par('ix_tiny',dat{s}.population,lkp,opt);
    
    % Initialise mixing proportions    
    prop            = ones(1,K)./K;   
    dat{s}.gmm.prop = prop;        

    if strcmpi(modality,'ct')
        
        %------------------------------------------------------------------
        % Initilisation of GMMs when data is CT
        %------------------------------------------------------------------                

        if 1
            Verbose = 0;
            IterMax = 4;

            X  = double(obs{s}{1});
            B  = double(obs{s}{2});           

            [~,MU,~,~,b,V,n] = spm_gmm(X,K,B,'BinWidth',1,'GaussPrior',opt.ct.GaussPrior,'PropPrior',dat{s}.gmm.prop,'Start','prior','Verbose',Verbose,'IterMax',IterMax); 

            post = {{MU,b},{V,n}}; % Posteriors        
        else
            post = {{opt.ct.GaussPrior{1},opt.ct.GaussPrior{2}},{opt.ct.GaussPrior{3},opt.ct.GaussPrior{4}}}; % Posteriors        
        end
    else        
        
        %------------------------------------------------------------------
        % Initilisation of GMMs when data is MRI
        %------------------------------------------------------------------
        
        X    = get_obs(dat{s},'do_scl',true,'mskonlynan',opt.seg.mskonlynan); % Do not subsample!   
        C    = size(X,2);        
        miss = get_par('missing_struct',X);

%         if opt.gmm.labels.use && isfield(dat{s},'label') && opt.gmm.labels.cm.isKey(dat{s}.population)                        
%             
%             %--------------------------------------------------------------
%             % Labels provided
%             %--------------------------------------------------------------        
%             
%             sort_pars = false;
%             
%             ix = opt.gmm.labels.cm(dat{s}.population);
%                                     
%             % Get labels
%             %--------------------------------------------------------------
%             labels = get_labels(dat{s},opt);                        
%             labs   = ix(ix~=0);
%             K_l    = numel(labs); % Number of labelled classes
%             K_nl   = K - K_l;     % Number of non-labelled classes
%             ix_bg  = max(ix) + 1;
%                 
%             % Get resps for non-labelled (background) voxels
%             msk_nl = labels{1} == ix_bg;
%             % figure; imshow3D(squeeze(reshape(msk_nl ,[dm_s])))
%             X1     = X(msk_nl,:);
%             
%             Z_nl = resp_from_kmeans(X1,K_nl - sum(ix_tiny));                        
%             X1   = [];
%             
%             % Get resps for labelled voxels
%             cnt = 1;
%             Z   = zeros([size(X,1) K],'single');
%             for k=1:K
%                 if ix_tiny(k)
%                     continue; 
%                 end
%                 
%                 if ix(k)==0
%                     Z(msk_nl,k) = Z_nl(:,cnt);
%                     cnt         = cnt + 1;
%                 else    
%                     Z(:,k)      = labels{1}==ix(k);
%                 end
%             end
%             Z_nl   = []; 
%             msk_nl = [];
%         else                      
            %--------------------------------------------------------------
            % No labels provided
            % Get responsibilities from kmeans labels
            %--------------------------------------------------------------
            
            sort_pars = false;
            
            Z = resp_from_kmeans(X,K - sum(ix_tiny));
                        
            nZ  = zeros([size(X,1) K],'single');
            cnt = 1;
            for k=1:K
                if ix_tiny(k)
                    continue; 
                end
                
                nZ(:,k) = Z(:,cnt);
                cnt     = cnt + 1;
            end
            Z  = nZ;
            nZ = [];
%         end              
        
        % Add tiny value where resps are zero
        Z       = Z + eps;
        Z       = bsxfun(@rdivide,Z,sum(Z,2));                             
%         show_classes_in_Z(Z,dm,K);
        
        % Build prior
        [SS0,SS1] = spm_gmm_lib('SuffStat',X,Z,1);
        
        b  = ones(1,K);
        n  = C*ones(1,K);
        MU = double(bsxfun(@rdivide,SS1,SS0));                
        V  = eye(C);
        V  = repmat(V,[1 1 K]);
        
        % Set prior
        pr = {MU,b,V,n};

        % Compute suffstats from responsibilities, then GMM parameters
        post = get_cluster(X,ones(size(X),'single'),dm,pr,miss,Z,'sort_pars',sort_pars);                
    end       

    dat{s}.gmm.cluster = post; % Posteriors
end

% Set intensity prior as sample mean (for MRI)
%--------------------------------------------------------------------------
populations   = spm_json_manager('get_populations',dat);
P             = numel(populations);
GaussPrior    = containers.Map;
lb_pr         = struct;                                        
lb_pr.KL_qVpV = 0;
lb_pr.ElnDetV = zeros(1,K);
for p=1:P  
    population0 = populations{p}.name;
    modality    = populations{p}.type;
    names       = get_channel_names(dat,populations{p});
        
    if strcmpi(modality,'CT')
        GaussPrior(population0) = opt.ct.GaussPrior;
    else
        % Use sample mean
        MU  = 0;
        b   = 0;
        W   = 0;
        n   = 0;
        cnt = 0;
        for s=1:S0
            population = dat{s}.population;

            if strcmp(population0,population)
                MU = MU + dat{s}.gmm.cluster{1}{1};
                b  = b  + dat{s}.gmm.cluster{1}{2};
                W  = W  + dat{s}.gmm.cluster{2}{1};
                n  = n  + dat{s}.gmm.cluster{2}{2};

                cnt = cnt + 1;
            end
        end

        MU = MU./cnt;
        b  = b./cnt;
        W  = W./cnt;
        n  = n./cnt;
                
        C  = size(MU,1);
        b  = ones(size(b));
        n  = C*ones(size(n));
        MU = zeros(size(MU));
        W  = eye(C,C);
        W  = repmat(W,[1 1 K]);
        
        GaussPrior(population0) = {MU,b,W,n,names,lb_pr,1:K};            
    end
end
model.GaussPrior = GaussPrior;
clear GaussPrior

% Compute initial estimate of hyper-parameters of VB-GMM
%--------------------------------------------------------------------------
model = update_GaussPrior(dat,model,opt);

% Set tissue indices to match between different populations
%--------------------------------------------------------------------------
for p=1:P
    population0 = populations{p}.name;  
        
    if init_ix.isKey(population0) 
        ix = init_ix(population0);    
    else                          
        continue;
    end
    
    % Adjust posteriors
    for s=1:S0
        population = dat{s}.population;
        
        if strcmp(population0,population)
            dat{s}.gmm.cluster{1}{1} = dat{s}.gmm.cluster{1}{1}(:,ix);
            dat{s}.gmm.cluster{1}{2} = dat{s}.gmm.cluster{1}{2}(ix);
            dat{s}.gmm.cluster{2}{1} = dat{s}.gmm.cluster{2}{1}(:,:,ix);
            dat{s}.gmm.cluster{2}{2} = dat{s}.gmm.cluster{2}{2}(ix);            
        end
    end
    
    % Adjust priors
    pr    = model.GaussPrior(population0);
    pr{1} = pr{1}(:,ix);
    pr{2} = pr{2}(ix);
    pr{3} = pr{3}(:,:,ix);
    pr{4} = pr{4}(ix);
    model.GaussPrior(population0) = pr;    
end
%==========================================================================

%==========================================================================   
function Z = resp_from_kmeans(X,K)
% Get initial labels using kmeans
L = spm_kmeans(X,K,'Distance','sqeuclidian', ... % cityblock/sqeuclidian
                   'Start','plus', ...
                   'Order','magnitude', ...
                   'Missing',true);

% Compute responsibilities from kmeans labels
Z = zeros([numel(L) K],'single');
for k=1:K
    Z(:,k) = L(:) == k; 
end
%==========================================================================   

%==========================================================================   
function chn_names = get_channel_names(dat,population)
S0    = numel(dat);
C     = population.C; 
name  = population.name; 
for s=1:S0
   population0 = dat{s}.population;
   
   if strcmpi(population0,name)              
       [~,~,~,~,~,~,~,chn_names] = obs_info(dat{s});              
       
       return
   end
end
%==========================================================================    