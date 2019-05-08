function [dat,bf] = update_bf(dat,obs,bf,template,dm,labels,miss,opt)
% FORMAT [dat,bf] = update_bf(dat,obs,bf,template,dm,labels,miss,opt)
% dat      - Subject's data structure (one subject)
% obs      - Observed image
% bf       - Bias fields
% template - Warped + softmaxed template
% dm       - Image dimensions
% labels   - Manual labels
% miss     - Missing data structure
% opt      - Options structure
%
% Update bias fields, one at a time, by Gauss-Newton optimisation.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Parse input
chan    = dat.bf.chan;      % Bias field encoding (C,B1,B2,B3,T)
cluster = dat.gmm.cluster;  % GMM posterior parameters
prop    = dat.gmm.prop;     % Tissue proportions
armijo  = dat.armijo.bf;    % Previous armijo factor
verbose = opt.verbose.bf;   % Verbosity level
part    = dat.gmm.part;     % Clusters to template mapping (lkp,mg)

% Parameters
kron             = @(a,b) spm_krutil(a,b); % Redefine kronecker product
cluster_scale    = cluster{2}{1};
cluster_scale_df = cluster{2}{2};
cluster_mean     = cluster{1}{1};
nb_clusters = size(cluster_mean,2);   % Number of clusters
nb_channels = numel(chan);            % Number of channels
const       = spm_gmm_lib('Const', cluster{1}, cluster{2}, miss.L); % GMM likelihood norm
ix_tiny     = get_par('ix_tiny',dat.population,part.lkp,opt);      % Labels to template mapping

for c=1:nb_channels % Loop over channels
    
    % Neighborhood part
    lnPzN = gmm_mrf('apply',dat.mrf);

    %----------------------------------------------------------------------
    % Compute objective function and its first and second derivatives
    %----------------------------------------------------------------------
    
    d3 = numel(chan(c).T); % Number of DCT parameters
    H  = zeros(d3,d3);     % Second derivatives w.r.t. DC coefficients
    gr = zeros(d3,1);      % First derivatives w.r.t. DC coefficients 
    
    for z=1:dm(3) % Loop over slices
        
        % Get slices
        [slice,ix] = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels);
        % slice.obs
        % slice.bf
        % slice.template
        % slice.code
        % slice.bin_var
        % slice.labels

        if dat.mrf.do && numel(lnPzN) > nb_clusters           
            lnPzNz = double(lnPzN(ix,:));
        else
            lnPzNz = lnPzN;
        end
    
        % Compute responsibilities and lb
        [Z,~,BX] = gmm_img('slice_resp_and_lb',slice,cluster{1},cluster{2},prop,part,miss,const,lnPzNz,ix_tiny);
        
        if dat.mrf.do        
            dat.mrf.oZ(:,:,z,:) = reshape(uint8((2^8)*cluster2template(Z,part)),[dm(1:2) 1 max(part.lkp)]);
        end
        
        slice_grad = zeros(dm(1:2));
        slice_hess = zeros(dm(1:2));
        
        % -----------------------------------------------------------------
        % For each combination of missing voxels
        for i=1:miss.nL

            % -------------------------------------------------------------
            % Get mask of missing modalities (with this particular code)
            code              = miss.L(i);
            observed_channels = spm_gmm_lib('code2bin', code, nb_channels);
            missing_channels  = ~observed_channels;
            if missing_channels(c), continue; end
            if isempty(slice.code), selected_voxels = ones(size(slice.code), 'logical');
            else,                   selected_voxels = (slice.code == code);
            end
            nb_channels_missing  = sum(missing_channels);
            nb_voxels_coded      = sum(selected_voxels);
            if nb_voxels_coded == 0, continue; end

            % -------------------------------------------------------------
            % Convert channel indices to observed indices
            mapped_c     = 1:nb_channels;
            mapped_c     = mapped_c(observed_channels);
            mapped_c     = find(mapped_c == c);
            cc           = mapped_c; % short alias
            
            selected_obs = BX(selected_voxels,observed_channels);
            gi = 0; % Gradient accumulated accross clusters
            Hi = 0; % Hessian accumulated accross clusters
            for k=1:nb_clusters

                % ---------------------------------------------------------
                % Compute expected precision (see GMM+missing data)
                Voo = cluster_scale(observed_channels,observed_channels,k);
                Vom = cluster_scale(observed_channels,missing_channels,k);
                Vmm = cluster_scale(missing_channels,missing_channels,k);
                Vmo = cluster_scale(missing_channels,observed_channels,k);
                Ao  = Voo - Vom*(Vmm\Vmo);
                Ao  = (cluster_scale_df(k)-nb_channels_missing) * Ao;
                MUo = cluster_mean(observed_channels,k);

                % ---------------------------------------------------------
                % Compute statistics
                gk = bsxfun(@minus, selected_obs, MUo.') * Ao(cc,:).';
                Hk = Ao(cc,cc);
                
                selected_resp = Z(selected_voxels,k);
                gk = bsxfun(@times, gk, selected_resp);
                Hk = bsxfun(@times, Hk, selected_resp);
                clear selected_resp

                % ---------------------------------------------------------
                % Accumulate across clusters
                gi = gi + gk;
                Hi = Hi + Hk;
                clear sk1x sk2x
            end
            
            % -------------------------------------------------------------
            % Deal with binning uncertainty
            binvar = 0;
            if numel(slice.bin_var) > 1
                binvar = slice.bin_var(selected_voxels,c);
            end
                
            % -------------------------------------------------------------
            % Multiply with bias corrected value (chain rule)
            gi = gi .* selected_obs(:,cc);
            Hi = Hi .* (selected_obs(:,cc).^2 + binvar);
            clear selected_obs binvar
            
            % -------------------------------------------------------------
            % Normalisation term
            gi = gi - 1;
            
            % -------------------------------------------------------------
            % Accumulate across missing codes
            slice_grad(selected_voxels) = slice_grad(selected_voxels) + gi;
            slice_hess(selected_voxels) = slice_hess(selected_voxels) + Hi;
            clear selected_voxels
        end

        b3 = chan(c).B3(z,:)';
        gr = gr + kron(b3,spm_krutil(slice_grad,chan(c).B1,chan(c).B2,0));
        H  = H  + kron(b3*b3',spm_krutil(slice_hess,chan(c).B1,chan(c).B2,1));
        clear slice_grad slice_hess b3 Z msk

    end % <-- Loop over slices
    
    % Inverse covariance of priors
    ICO = chan(c).C;         

    % Gauss-Newton update of bias field parameters
    Update = reshape((H + ICO)\(gr + ICO*chan(c).T(:)),size(chan(c).T));
    clear H gr

    % Line-search
    %------------------------------------------------------------------
    obf          = bf;
    oT           = chan(c).T;
    bf_reg       = dat.lb.bf_reg(end,:);
    nline_search = opt.nline_search.bf;
    for line_search=1:nline_search

        % Update bias-field parameters
        chan(c).T = chan(c).T - armijo(c)*Update;

        % Compute new bias-field (only for channel c)
        [bf,bf_reg] = get_bf(chan,dm,bf,c,bf_reg);                                               

        % Compute new lower bound                      
        lblnDetbf                   = bf;
        lblnDetbf(isnan(lblnDetbf)) = 1;        
        lblnDetbf                   = log(prod(lblnDetbf,2));  
        lblnDetbf                   = sum(lblnDetbf);

        [dlb,~,dat.mrf] = gmm_img('img_lb_and_mom',obs,bf,[],template,labels,prop,cluster{1},cluster{2},miss,part,dm,dat.mrf,ix_tiny,{'bf',obf});                      

        % Check new lower bound
        if (dlb.X + lblnDetbf + sum(bf_reg)) > (dat.lb.X(end) + dat.lb.lnDetbf(end) + sum(dat.lb.bf_reg(end,:)))              
            armijo(c) = min(armijo(c)*1.25,1);

            dat.lb.X(end + 1)        = dlb.X;                    
            dat.lb.lnDetbf(end + 1)  = lblnDetbf;
            dat.lb.bf_reg(end + 1,:) = bf_reg;                                                                                                         
            
            break;
        else                                
            armijo(c) = armijo(c)*0.5;
            chan(c).T = oT;

            if line_search == nline_search                                        
                bf = obf;                                         
            end
        end
    end
    clear oT Update obf
    
    dat.lb = check_convergence('bf',dat.lb,c,verbose,armijo(c));
end % <-- Loop over channels

if verbose >= 3
    show_bf_and_ivel(obs,dm,bf);
end
            
% Get DC component
dc.ln = zeros(1,nb_channels);
for c=1:nb_channels
    dc.ln(c) = chan(c).T(1,1,1);
end
dc.int = scl_from_bf(chan);

% Set output
dat.bf.chan   = chan;
dat.bf.dc     = dc;
dat.armijo.bf = armijo;
%==========================================================================