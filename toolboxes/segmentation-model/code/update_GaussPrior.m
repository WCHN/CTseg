function model = update_GaussPrior(dat,model,opt)
% FORMAT model = update_GaussPrior(dat,model,opt)
% dat       - Subjects data structure
% model     - Model structure
% opt       - Options structure
%
% Update Gaussian Mixture priors
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging


dir_model   = opt.dir_model;                    % Directory were model file sare stored
verbose     = opt.gmm.hist.verbose;             % Verbosity
constrained = opt.gmm.GaussPrior.constrained;   % Constrain variances to be similar
GaussPrior  = model.GaussPrior;                 % Prior parameters dictionary {mean df scale df}
ModSharing  = opt.gmm.GaussPrior.mods;

S0          = numel(dat);                              % Number of subjects
populations = spm_json_manager('get_populations',dat); % list of populations
P           = numel(populations);                      % Number of populations

if ~isempty(ModSharing)
    %----------------------------------------------------------------------
    % First we update populations which share the prior
    %----------------------------------------------------------------------

    for m0=1:numel(ModSharing)
        
        NameMod = ModSharing{m0};
        
        clusters = {};    
        cnt      = 1;     
        has_mod  = false; 

        for p=1:P % Iterate over populations
            population0 = populations{p}.name;
            modality    = get_modality_name(dat,population0);

            if ~strcmpi(modality,NameMod)
                continue
            end

            has_mod = true;

            pr = GaussPrior(population0);

            % Get lkp
            for s=1:S0
                population = dat{s}.population;

                if strcmp(population0,population)
                    lkp = dat{s}.gmm.part.lkp; % Map GMM cluster to Template classes
                    break
                end
            end

            for s=1:S0 % Iterate over subjects      
                population = dat{s}.population; 

                if strcmpi(population0,population)
                    clusters{cnt} = dat{s}.gmm.cluster;
                    cnt           = cnt + 1;
               end
            end
        end

        if has_mod
            [pr(1:4),extras] = spm_gmm_lib('UpdateHyperPars',clusters,pr(1:4), ...
                                           'constrained',constrained,...
                                           'figname',modality,...
                                           'verbose',verbose,...
                                           'lkp',lkp);    


            lb_pr                  = pr{6};                                        
            lb_pr.KL_qVpV(end + 1) = extras.lb;  % KL(q(V0) || p(V0))
            lb_pr.ElnDetV          = extras.ldW; % E[ln|V0|]    
            pr{6}                  = lb_pr;    

            for p=1:P % Iterate over populations
                population0 = populations{p}.name;
                modality    = get_modality_name(dat,population0);

                if ~strcmpi(modality,NameMod)
                    continue
                end

                GaussPrior(population0) = pr; % Update population prior
            end
        end
    end
end

%--------------------------------------------------------------------------
% Then we update individual GaussPriors
%--------------------------------------------------------------------------

for p=1:P % Iterate over populations
    population0 = populations{p}.name;
    modality    = get_modality_name(dat,population0);
        
    if any(strcmp(ModSharing,modality)) && ~isempty(modality) 
        % Skip this modality
        continue
    end
    
    pr           = GaussPrior(population0);
    clusters     = {};  
    cnt          = 1;

    % Get lkp
    for s=1:S0
        population = dat{s}.population;

        if strcmp(population0,population)
            lkp = dat{s}.gmm.part.lkp;
            break
        end
    end

    for s=1:S0 % Iterate over subjects      
        population = dat{s}.population; 

        if strcmpi(population0,population)
            clusters{cnt} = dat{s}.gmm.cluster;
            cnt           = cnt + 1;
       end
    end

    [pr(1:4),extras] = spm_gmm_lib('UpdateHyperPars',clusters,pr(1:4), ...
                                   'constrained',constrained,...
                                   'figname',population0,...
                                   'verbose',verbose,...
                                   'lkp',lkp);    


    lb_pr                  = pr{6};                                        
    lb_pr.KL_qVpV(end + 1) = extras.lb;  % KL(q(V0) || p(V0))
    lb_pr.ElnDetV          = extras.ldW; % E[ln|V0|]    
    pr{6}                  = lb_pr;    

    GaussPrior(population0) = pr;
end

%--------------------------------------------------------------------------
% Save prior parameters dictionary
%--------------------------------------------------------------------------

fname = fullfile(dir_model,'GaussPrior.mat');
save(fname,'GaussPrior');

model.GaussPrior = GaussPrior;

deal_figs(model);
%==========================================================================