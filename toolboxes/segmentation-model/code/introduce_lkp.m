function varargout = introduce_lkp(dat,model,opt,varargin)
% Introduce several GMM clusters (if needed) for each template class.
% This should only be done once.
%
% 1) 'Segment' case: only one subject
% FORMAT [dat,GaussPrior] = introduce_lkp(dat,model,opt,GaussPrior)
%
% 2) 'Train' case: several subjects
% FORMAT [dat,model] = introduce_lkp(dat,model,opt,iter)
%
% dat        - Subject data structure
% model      - Model structure (with field GaussPrior)
% opt        - Options structure
% GaussPrior - Cell describing the Gaussian mixture prior
%               {1} [PxKp]   Means
%               {2} [1xKp]   Mean d.f.
%               {3} [PxPxKp] Scale matrices
%               {4} [1xKp]   Precision d.f.
%               {5}
%               {6} [struct] Useful values (ElnDetV, ...)
%               {7} [1xKp]   Partitioning of the K template classes into Kp
%                            GMM classes
% iter       - Current EM iteration (only introduce if opt.start_it.do_mg)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

K  = opt.template.K;
S0 = numel(dat);   

if S0 == 1
    % ---------------------------------------------------------------------
    % Case 'Segment'
    GaussPrior = varargin{1};
    
    modality = dat.modality{1}.name; 
    lkp      = get_par('lkp',modality,opt);
        
    if numel(lkp) == K
        % As many GMM clusters as Template classes
        varargout{1} = dat;
        varargout{2} = model.GaussPrior(dat.population);
    
        return
    end
    
    % Modify posteriors
    gmm{1} = dat.gmm.cluster{1}{1}; % Mean
    gmm{2} = dat.gmm.cluster{1}{2}; % Mean df
    gmm{3} = dat.gmm.cluster{2}{1}; % Scale matrix
    gmm{4} = dat.gmm.cluster{2}{2}; % Precision df

    [gmm,mg] = spm_gmm_lib('extras', 'more_gmms', gmm, lkp);           

    dat.gmm.cluster{1}{1} = gmm{1};
    dat.gmm.cluster{1}{2} = gmm{2};
    dat.gmm.cluster{2}{1} = gmm{3};
    dat.gmm.cluster{2}{2} = gmm{4};

    dat.gmm.part.lkp = lkp; % New cluster to template mapping
    dat.gmm.part.mg  = mg;  % Within template weight per cluster
    
    % Modify GaussPrior   
    GaussPrior{7}         = lkp;
    GaussPrior{6}.ElnDetV = zeros(1,numel(lkp));                
    GaussPrior(1:4)       = spm_gmm_lib('extras', 'more_gmms', GaussPrior(1:4), lkp);        
    
    varargout{1} = dat;
    varargout{2} = GaussPrior;
else
    % ---------------------------------------------------------------------
    % Case 'Train'
    it_mod = varargin{1};
    
    if it_mod == opt.start_it.do_mg
        populations  = spm_json_manager('get_populations',dat);
        P            = numel(populations);

        for p=1:P  
            population0 = populations{p}.name;
            modality    = populations{p}.type;
            GaussPrior  = model.GaussPrior(population0);
            lkp         = get_par('lkp',modality,opt);

            if numel(lkp) == K
                continue
            end

            mg               = ones(1,numel(lkp));
            for k=1:max(lkp)
                kk           = sum(lkp == k);
                mg(lkp == k) = 1/kk;
            end   

            % Modify subject-specific posteriors
            for s=1:S0
                population = dat{s}.population;

                if strcmp(population0,population)
                    dat{s}.gmm.part.lkp = lkp;
                    dat{s}.gmm.part.mg  = mg;

                    gmm{1} = dat{s}.gmm.cluster{1}{1}; % Mean
                    gmm{2} = dat{s}.gmm.cluster{1}{2}; % Mean df
                    gmm{3} = dat{s}.gmm.cluster{2}{1}; % Scale matrix
                    gmm{4} = dat{s}.gmm.cluster{2}{2}; % Precision df

                    gmm = spm_gmm_lib('extras', 'more_gmms', gmm, lkp);           

                    dat{s}.gmm.cluster{1}{1} = gmm{1};
                    dat{s}.gmm.cluster{1}{2} = gmm{2};
                    dat{s}.gmm.cluster{2}{1} = gmm{3};
                    dat{s}.gmm.cluster{2}{2} = gmm{4};
                end
            end

            % Modify GaussPrior   
            GaussPrior{7}         = lkp;
            GaussPrior{6}.ElnDetV = zeros(1,numel(lkp));                                        
            GaussPrior(1:4)       = spm_gmm_lib('extras', 'more_gmms', GaussPrior(1:4), lkp);           

            model.GaussPrior(population0) = GaussPrior;
        end
    end
    
    varargout{1} = dat;
    varargout{2} = model;
end
%==========================================================================