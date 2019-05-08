function dat = meancorrect_bf(dat,GaussPrior,opt)
% FORMAT dat = meancorrect_bf(dat,GaussPrior,opt)
% dat        - Subjects data structure
% GaussPrior - Cell of GMM prameters
% opt        - Options structure
%
% Zero-centre bias field mean-value across subjects.
% Since GMM parameters are stongly impacted by the bias mean value, they
% are updated after zero-centering.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if opt.bf.mc_bf && opt.bf.do
    % Parameters
    Verbose     = opt.bf.mc_bf_verbose;
    populations = spm_json_manager('get_populations',dat);
    P           = numel(populations);
    S0          = numel(dat);
    ModSharing  = opt.gmm.GaussPrior.mods;

    if ~isempty(ModSharing)
        
        for m0=1:numel(ModSharing)
            
            NameMod = ModSharing{m0};
            
            % Iterate over populations
            sm_dc_ln  = 0;
            sm_dc_int = 0;
            S_p       = 0;
            for p=1:P
                population0 = populations{p}.name;
                modality    = get_modality_name(dat,population0);

                if ~strcmpi(modality,NameMod) || strcmpi(modality,'CT') 
                    % Skip this modality
                    continue
                end

                population_p = populations{p}.name;      % Get population

                % For the current population, compute the mean DC component of the
                % bias field (in log and intensity space)
                
                for s=1:S0
                    population = dat{s}.population;

                    if strcmp(population_p,population)
                        sm_dc_ln  = sm_dc_ln  + dat{s}.bf.dc.ln;
                        sm_dc_int = sm_dc_int + dat{s}.bf.dc.int;

                        S_p = S_p + 1;
                    end
                end
            end
            
            mn_dc_ln  = sm_dc_ln./S_p;
            mn_dc_int = sm_dc_int./S_p;

            % Number of channels
            C = numel(mn_dc_ln); 

            if Verbose
                % Some verbose
                fprintf('mn_dc_int (%s) = [',modality);
                for c=1:C - 1
                    fprintf('%4.2f ',mn_dc_int(c));
                end
                fprintf('%4.2f',mn_dc_int(C))
                fprintf(']\n');
            end

            for p=1:P
                population0 = populations{p}.name;
                modality    = get_modality_name(dat,population0);

                if ~strcmpi(modality,NameMod) || strcmpi(modality,'CT') 
                    % Skip this modality
                    continue
                end

                population_p = populations{p}.name;      % Get population
                pr           = GaussPrior(population_p); % Get prior
            
                % For the current population, subtract the mean (log) DC component to
                % each subject's bias field DC component. Then correct the GMM
                % parameters via each subject's suffstats.
                for s=1:S0
                    population = dat{s}.population;

                    if strcmp(population_p,population)  
                        % Subtract mean DC components
                        for c=1:C
                            dat{s}.bf.chan(c).T(1,1,1) = dat{s}.bf.chan(c).T(1,1,1) - mn_dc_ln(c);
                        end

                        scl = (1./mn_dc_int)';
        %                 scl = (mn_dc_int)';

                        % Correct GMM parameters via suffstats
                        dSS0 = dat{s}.gmm.mom.SS0;
                        dSS1 = dat{s}.gmm.mom.SS1;
                        dSS2 = dat{s}.gmm.mom.SS2;

                        Missing = iscell(dSS0);

                        if Missing
                            SS2b    = dat{s}.gmm.mom.SS2b;
                            L       = dat{s}.gmm.mom.L;
                            cluster = dat{s}.gmm.cluster;
                            MU      = cluster{1}{1};
                            K       = size(MU,2);
                            A       = bsxfun(@times,cluster{2}{1},reshape(cluster{2}{2},[1 1 K]));

                            for i=2:numel(L)
                                c        = L(i);
                                observed = code2bin(c, C);

                                scl1 = scl(observed);

                                dSS1{i} = bsxfun(@times,dSS1{i},scl1);
                                for k=1:size(dSS2{i},3)
                                    dSS2{i}(:,:,k) = (scl1*scl1').*dSS2{i}(:,:,k);
                                end     
                            end

                            [SS0,SS1,SS2] = spm_gmm_lib('SuffStat', 'infer', dSS0, dSS1, dSS2, {MU,A}, L);
                            SS2           = SS2 + SS2b;
                        else
                            dSS1 = bsxfun(@times,dSS1,scl);
                            for k=1:size(dSS2,3)
                                dSS2(:,:,k) = (scl*scl').*dSS2(:,:,k);
                            end         
                            SS0 = dSS0;
                            SS1 = dSS1;
                            SS2 = dSS2;
                        end   

                        dat{s}.gmm.mom.SS1 = dSS1;
                        dat{s}.gmm.mom.SS2 = dSS2;

                        [MU,~,b,V,n] = spm_gmm_lib('UpdateClusters', ...
                                                   SS0, SS1, SS2, pr);    

                        dat{s}.gmm.cluster = {{MU,b},{V,n}};
                    end
                end
            end
    
        end
    end
       
    % Iterate over populations
    for p=1:P
        population0 = populations{p}.name;
        modality    = get_modality_name(dat,population0);

        if (any(strcmp(ModSharing,modality)) && ~isempty(modality)) || strcmpi(modality,'CT') 
            % Skip this modality
            continue
        end

        population_p = populations{p}.name;      % Get population
        pr           = GaussPrior(population_p); % Get prior
        
        % For the current population, compute the mean DC component of the
        % bias field (in log and intensity space)
        sm_dc_ln  = 0;
        sm_dc_int = 0;
        S_p       = 0;
        for s=1:S0
            population = dat{s}.population;

            if strcmp(population_p,population)
                sm_dc_ln  = sm_dc_ln  + dat{s}.bf.dc.ln;
                sm_dc_int = sm_dc_int + dat{s}.bf.dc.int;

                S_p = S_p + 1;
            end
        end
        mn_dc_ln  = sm_dc_ln./S_p;
        mn_dc_int = sm_dc_int./S_p;

        % Number of channels
        C = numel(mn_dc_ln); 

        if Verbose
            % Some verbose
            fprintf('mn_dc_int (%s) = [',population_p);
            for c=1:C - 1
                fprintf('%4.2f ',mn_dc_int(c));
            end
            fprintf('%4.2f',mn_dc_int(C))
            fprintf(']\n');
        end

        % For the current population, subtract the mean (log) DC component to
        % each subject's bias field DC component. Then correct the GMM
        % parameters via each subject's suffstats.
        for s=1:S0
            population = dat{s}.population;

            if strcmp(population_p,population)  
                % Subtract mean DC components
                for c=1:C
                    dat{s}.bf.chan(c).T(1,1,1) = dat{s}.bf.chan(c).T(1,1,1) - mn_dc_ln(c);
                end

                scl = (1./mn_dc_int)';
%                 scl = (mn_dc_int)';

                % Correct GMM parameters via suffstats
                dSS0 = dat{s}.gmm.mom.SS0;
                dSS1 = dat{s}.gmm.mom.SS1;
                dSS2 = dat{s}.gmm.mom.SS2;

                Missing = iscell(dSS0);

                if Missing
                    SS2b    = dat{s}.gmm.mom.SS2b;
                    L       = dat{s}.gmm.mom.L;
                    cluster = dat{s}.gmm.cluster;
                    MU      = cluster{1}{1};
                    K       = size(MU,2);
                    A       = bsxfun(@times,cluster{2}{1},reshape(cluster{2}{2},[1 1 K]));

                    for i=2:numel(L)
                        c        = L(i);
                        observed = code2bin(c, C);

                        scl1 = scl(observed);

                        dSS1{i} = bsxfun(@times,dSS1{i},scl1);
                        for k=1:size(dSS2{i},3)
                            dSS2{i}(:,:,k) = (scl1*scl1').*dSS2{i}(:,:,k);
                        end     
                    end

                    [SS0,SS1,SS2] = spm_gmm_lib('SuffStat', 'infer', dSS0, dSS1, dSS2, {MU,A}, L);
                    SS2           = SS2 + SS2b;
                else
                    dSS1 = bsxfun(@times,dSS1,scl);
                    for k=1:size(dSS2,3)
                        dSS2(:,:,k) = (scl*scl').*dSS2(:,:,k);
                    end         
                    SS0 = dSS0;
                    SS1 = dSS1;
                    SS2 = dSS2;
                end   
                
                dat{s}.gmm.mom.SS1 = dSS1;
                dat{s}.gmm.mom.SS2 = dSS2;

                [MU,~,b,V,n] = spm_gmm_lib('UpdateClusters', ...
                                           SS0, SS1, SS2, pr);    

                dat{s}.gmm.cluster = {{MU,b},{V,n}};
            end
        end
    end
end
%==========================================================================

% =========================================================================
function bin = code2bin(code, length)
% FORMAT bin = spm_gmm_lib('code2bin', code, length)

bin = dec2bin(code,length) == '1';
bin = bin(end:-1:1);
% =========================================================================