function [dat,model,opt] = load_model(dat,opt)
% FORMAT [dat,model,opt] = load_model(dat,opt)
% dat   - Subjects data structure
% opt   - Options structure
% model - Model structure
%
% Load model from disk.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

S0    = numel(dat);
model = cell(1,S0);
for s=1:S0
    
    fprintf('---------------------------------------\n');
    fprintf('Initialising model\n');

    %--------------------------------------------------------------------------
    % Set template
    %--------------------------------------------------------------------------

    if isfield(opt.template,'pth_template') && exist(opt.template.pth_template,'file')

        model{s}.template.nii = nifti(opt.template.pth_template);
        d              = model{s}.template.nii.dat.dim(1:3);
        K              = model{s}.template.nii.dat.dim(4);
        opt.template.K = K; 

        % Get values to be used for FOV voxels when warping 
        model{s} = init_template_bg(model{s},opt);

        fprintf('1 | Loaded template with %i classes.\n',K);        
    else
        model{s} = init_uniform_template(dat,opt);
        K        = opt.template.K;

        % Uninformative template -> don't do registration
        opt.reg.do_aff = false;
        opt.reg.do_nl  = false;

        fprintf('1 | Initialised uniform template with %i classes.\n',K);    
    end

    %--------------------------------------------------------------------------
    % Set PropPrior
    %--------------------------------------------------------------------------

    if isfield(opt.gmm,'pth_PropPrior') && exist(opt.gmm.pth_PropPrior,'file')   
        fprintf('2 | Proportion prior for %s found...loaded.\n',dat{s}.population);
        
        load(opt.gmm.pth_PropPrior,'PropPrior');
        model{s}.PropPrior = PropPrior;    
    else
        fprintf('2 | Proportion prior for %s not found...using uninformative.\n',dat{s}.population);
        
        alpha = ones(1,K)*opt.prop.reg;

        model{s}.PropPrior.alpha = alpha;
        model{s}.PropPrior.norm  = 0;
    end
    
    %--------------------------------------------------------------------------
    % Set GaussPrior
    %--------------------------------------------------------------------------

    % Subject info
    [X,~,~,~,~,~,~,mn,mx]     = get_obs(dat{s},'do_scl',true,'mskonlynan',opt.seg.mskonlynan); % Do not subsample!
    [~,~,~,C,~,~,~,chn_names] = obs_info(dat{s});    

    % Uninformative VBGMM hyper-parameters
    b  = ones(1,K);
    n  = C*ones(1,K);
    MU = double(repmat(nanmedian(X)',[1 K]));        
    V  = zeros(C,C);
    for c=1:C
        vr     = ((mx(c) - mn(c))./(K)).^2;
        V(c,c) = 1./(n(1).*vr);
    end
    V  = repmat(V,[1 1 K]);

    if isfield(opt.gmm,'pth_GaussPrior') && exist(opt.gmm.pth_GaussPrior,'file')    

        load(opt.gmm.pth_GaussPrior,'GaussPrior');
        model{s}.GaussPrior = GaussPrior;
        clear GaussPrior

        if model{s}.GaussPrior.isKey(dat{s}.population)
            % Population exists in prior dictionary
            %----------------------------------------------------------------------

            % Adjust prior so that image channels maps to correct indices in the prior
            pr   = model{s}.GaussPrior(dat{s}.population);    
            nam1 = pr{5}; % Order of channels in intensity prior
            ix   = zeros(1,C);
            for c=1:C
                nam0 = chn_names{c}; % Order of channels in observed data

                % Figure out mapping from image to prior
                cnt = 1;
                for c1=1:C
                    if ischar(nam1)
                        nam11 = nam1;
                    else
                        nam11 = nam1{c1};
                    end

                    if strcmpi(nam11, nam0)
                        ix(c) = cnt;

                        continue
                    end

                    cnt = cnt + 1;        
                end
            end

            % Correct prior
            pr{1} = pr{1}(ix,:);
            pr{3} = pr{3}(ix,ix,:);
            pr{5} = pr{5}(ix);    

            no_GaussPrior = false(1,C);

            fprintf('3 | Intensity prior for %s found...loaded.\n',dat{s}.population);
        else
            % Population does not exist in prior dictionary, check if channel
            % exists in other population
            %----------------------------------------------------------------------

            fprintf('3 | Intensity prior for %s not found.\n',dat{s}.population);

            pr    = cell(1,7);
            pr{1} = zeros(C,K);
            pr{2} = zeros(1,K);
            pr{3} = zeros(C,C,K);
            pr{4} = zeros(1,K);
            pr{5} = cell(1,C);  
            pr{7} = get_par('lkp',dat{s}.modality{1}.name,opt);

            no_GaussPrior = true(1,C);
            keys          = model{s}.GaussPrior.keys;    
            for c1=1:C % loop over observed channel names
                nam0      = chn_names{c1};
                pr{5}{c1} = nam0;

                for k=1:numel(keys) % loop over available priors
                    pr1  = model{s}.GaussPrior(keys{k});
                    nam1 = pr1{5};

                    for c2=1:numel(nam1) % loop over channels in pulled out GaussPrior
                        if strcmpi(nam0,nam1{c2})
                           % Prior found

                           if size(pr{1},2) == size(pr1{1},2)
                               % Same number of tissues
                               pr{1}(c1,:)    = pr1{1}(c2,:);
                               pr{2}          = pr1{2};                           
                               pr{3}(c1,c1,:) = pr1{3}(c2,c2,:);
                               pr{4}          = pr1{4};

                               no_GaussPrior(c1) = false;

                               fprintf('2 | Channel %s will instead be loaded from population %s.\n',nam0,keys{k});

                               break;
                           end
                        end
                    end

                    if ~no_GaussPrior(c1)
                       break; 
                    end
                end    
            end
        end

        if any(no_GaussPrior == true)
            % Could not find a GaussPrior -> initialise as uniformative
            %----------------------------------------------------------------------        

            if sum(no_GaussPrior) == C
                % No channels have been initialised -> use uniformative b and n
                pr{2} = b;
                pr{4} = n;
            end

            for c=find(no_GaussPrior == true)
                nam0 = chn_names{c};

                if strcmpi(nam0,'CT')
                    pr(1:4) = opt.ct.GaussPrior(1:4);

                    fprintf('2 | Channel %s initialised from default CT prior.\n',nam0);
                else
                    pr{1}(c,:)   = MU(c,:);        
                    pr{3}(c,c,:) = V(c,c,:);        

                    fprintf('2 | Channel %s initialised as non-informative.\n',nam0);
                end                        
            end
        end
    else
        model{s}.GaussPrior = containers.Map;

        pr{1} = MU;
        pr{2} = b;
        pr{3} = V;
        pr{4} = n;    
        pr{5} = chn_names;    
        pr{7} = get_par('lkp',dat{s}.modality{1}.name,opt);

        if strcmpi(dat{s}.modality{1}.name,'CT')
            pr = opt.ct.GaussPrior;
        end

        fprintf('2 | Intensity prior initialised as uniformative.\n');
    end

    % Collapse prior
    [lkp0,mg] = get_par('lkp',dat{s}.modality{1}.name,opt);
    lkp       = pr{7};
    if ~isequal(lkp0,lkp)        
        [pr(1:4),mg] = spm_gmm_lib('extras','collapse_gmms',pr(1:4),lkp);
        lkp          = 1:numel(mg); 
    end

    % LB part when having prior on V0
    lb_pr         = struct;                                        
    lb_pr.KL_qVpV = 0;
    lb_pr.ElnDetV = zeros(1,numel(lkp));
    pr{6}         = lb_pr;
    pr{7}         = lkp;

    model{s}.GaussPrior(dat{s}.population) = pr;   
   
    fprintf('---------------------------------------\n');

    %--------------------------------------------------------------------------
    % Cluster proportions
    %--------------------------------------------------------------------------
    
    dat{s}.gmm.part.mg  = mg; 
    dat{s}.gmm.part.lkp = lkp;    
    
    %--------------------------------------------------------------------------
    % Tissue proportions
    %--------------------------------------------------------------------------

    dat{s}.gmm.prop = zeros([1 K]); 
end
%==========================================================================