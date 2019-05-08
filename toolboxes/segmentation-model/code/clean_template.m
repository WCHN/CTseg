function model = clean_template(model,opt,Verbose)
% FORMAT model = clean_template(model,opt)
% model   - Model structure
% opt     - Options structure
% Verbose - Show masking results
%
% Ad-hoc template clean-up routine, to remove non-brain voxels from classes
% that should only contain brain and non-air voxels from classes that
% should only contain air.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 3, Verbose = false; end

% Parameters
BrainClass  = opt.template.clean.brain;
AirClass    = opt.template.clean.air;
LesionClass = opt.template.clean.les;
Tiny        = eps('single');
ThrldBrain  = opt.template.clean.val_brain;
ThrldAir    = opt.template.clean.val_air;
FigNum      = 666;
dm          = model.template.nii.dat.dim;
K           = dm(4);
dil_er      = opt.template.clean.dil_er;
it_dil_er   = opt.template.clean.it_dil_er;
SmoothSig   = 1;

% Get soft-maxed template
template = single(model.template.nii.dat(:,:,:,:));
template = spm_matcomp('softmax',template);  

if ~isempty(BrainClass)
    %----------------------------------------------------------------------
    % First attempt to clean up within brain classes
    %----------------------------------------------------------------------

    nfigs_sp1 = 2*numel(BrainClass) + 4;
    nr        = floor(sqrt(nfigs_sp1));
    nc        = ceil(nfigs_sp1/nr);      
    fig_cnt   = 1;            
    
    % Get brain class(es)
    cls = template(:,:,:,BrainClass);
    msk = sum(cls,4) > ThrldBrain;

    if Verbose
        % Some verbose
        figure(FigNum)
        for k=1:size(cls,4)
            subplot(nr,nc,fig_cnt); fig_cnt = fig_cnt + 1;            
            img = cls(:,:,:,k);
            imagesc3d(img); axis off; drawnow
        end
        subplot(nr,nc,fig_cnt); fig_cnt = fig_cnt + 1;
        imagesc3d(msk); axis off; drawnow
    end

    % Find largest connected component in brain mask
    C   = bwlabeln(msk);
    unqC  = unique(C);
    vls = zeros([1 numel(unqC)]);
    for i=1:numel(unqC)
        vls(i) = sum(sum(sum(C == unqC(i))));
    end
    [~,ix] = sort(vls);
    ix     = (ix(end - 1) - 1);
    msk    = C == ix;

    if Verbose
        % Some verbose
        subplot(nr,nc,fig_cnt); fig_cnt = fig_cnt + 1;
        imagesc3d(msk); axis off; drawnow
    end

    % Fill holes
    msk = imfill(msk,'holes');

    if Verbose
        % Some verbose
        subplot(nr,nc,fig_cnt); fig_cnt = fig_cnt + 1;
        imagesc3d(msk); axis off; drawnow
    end

    if Verbose
        % Some verbose
        subplot(nr,nc,fig_cnt); fig_cnt = fig_cnt + 1;    
    end
    
    if dil_er
        % Dilate/erode
        se  = strel('sphere',1);
        for i=1:it_dil_er
            msk = imdilate(msk, se);    

            if Verbose
                % Some verbose        
                imagesc3d(msk); axis off; drawnow
            end
        end        

        for i=1:it_dil_er
            msk = imerode(msk, se);

            if Verbose
                % Some verbose        
                imagesc3d(msk); axis off; drawnow
            end
        end       
    end
    
    if Verbose
        % Some verbose
        figure(FigNum)
        for k=1:size(cls,4)
            subplot(nr,nc,fig_cnt); fig_cnt = fig_cnt + 1;
            img = cls(:,:,:,k);
            img = msk.*img;
            imagesc3d(img); axis off; drawnow
        end
    end

    % Update template with mask information
    template = adjust_template_brain(template,msk,BrainClass,LesionClass,Tiny,SmoothSig);

    if Verbose
        % Some verbose
        figure(FigNum + 1)
        nr = floor(sqrt(K));
        nc = ceil(K/nr);      
        for k=1:K   
            subplot(nr,nc,k)
            imagesc3d(template(:,:,:,k)); axis off; drawnow
        end
    end
end

if ~isempty(AirClass)
    %----------------------------------------------------------------------
    % Then attempt to clean up air class
    %----------------------------------------------------------------------

    nfigs_sp1 = 2*numel(AirClass) + 3;
    nr        = floor(sqrt(nfigs_sp1));
    nc        = ceil(nfigs_sp1/nr);      
    fig_cnt   = 1;      
    
    % Get air class(es)
    cls = template(:,:,:,AirClass);
    msk = sum(cls,4) > ThrldAir;

    if Verbose
        % Some verbose
        figure(FigNum + 2)
        for k=1:size(cls,4)
            subplot(nr,nc,fig_cnt); fig_cnt = fig_cnt + 1;  
            img = cls(:,:,:,k);
            imagesc3d(img); axis off; drawnow
        end
        subplot(nr,nc,fig_cnt); fig_cnt = fig_cnt + 1;  
        imagesc3d(msk); axis off; drawnow
    end

    % Do connected component analysis
    C    = bwlabeln(msk);
    unqC = unique(C(:));
    unqC = unqC(:)';
    nC   = numel(unqC);    
        
    % Find the largest volume that is not air
    mx_c      = 5;            % Max number of volumes to search through    
    mx_c      = min(nC,mx_c);
    unqC      = unqC(1:mx_c);
    vols_brd  = [];
    for c=unqC        
        msk_c     = C == c;                   % Get a component mask
        vol_brd_c = get_border_volume(msk_c); % Compute volumes of its borders        
        vols_brd  = [vols_brd vol_brd_c];
    end
    
    % Pick the largest volume from vols_brd, this should be the air class
    [~,ix_air] = max(vols_brd);
    msk        = C == unqC(ix_air);
    
    % Fill holes
    msk = ~msk; % Invert
    msk = imfill(msk,'holes');
    msk = ~msk; % Invert
    
    % Dilate
    se  = strel('sphere',1);
    msk = imdilate(msk, se);    
    
    if Verbose
        % Some verbose
        subplot(nr,nc,fig_cnt); fig_cnt = fig_cnt + 1;  
        imagesc3d(msk); axis off; drawnow
    end

    % Update template with mask information
    template = adjust_template_air(template,msk,AirClass,Tiny,SmoothSig);

    if Verbose
        % Some verbose
        figure(FigNum + 3)
        nr = floor(sqrt(K));
        nc = ceil(K/nr);      
        for k=1:K   
            subplot(nr,nc,k)
            imagesc3d(template(:,:,:,k)); axis off; drawnow
        end
    end
end

if nargout > 0 && (~isempty(BrainClass) || ~isempty(AirClass))
    % Take log and write to template NIfTI
    model.template.nii.dat(:,:,:,:) = log(template);
end
%==========================================================================

%==========================================================================
function template = adjust_template_brain(template,msk,BrainClasses,LesionClass,Tiny,SmoothSig)
dm       = size(template);
K        = dm(4);
msk      = reshape(msk,[prod(dm(1:3)) 1]);
template = reshape(template,[prod(dm(1:3)) K]);
for k=1:K    
    if ismember(k,BrainClasses) || ismember(k,LesionClass)
        template(~msk,k) = Tiny;
        
        if SmoothSig > 0
            % Smooth mask a little bit
            template          = reshape(template,[dm(1:3) K]);
            template(:,:,:,k) = imgaussfilt3(template(:,:,:,k),SmoothSig);
            template          = reshape(template,[prod(dm(1:3)) K]);
        end
    end    
end

% Renormalise
template = reshape(template,[dm(1:3) K]);
template = bsxfun(@rdivide,template,sum(template,4) + eps);
%==========================================================================

%==========================================================================
function template = adjust_template_air(template,msk,AirClass,Tiny,SmoothSig)
dm       = size(template);
K        = dm(4);
msk      = reshape(msk,[prod(dm(1:3)) 1]);
template = reshape(template,[prod(dm(1:3)) K]);
for k=1:K    
    if ismember(k,AirClass)
        template(msk,k) = 1 - Tiny;

        if SmoothSig > 0
            % Smooth mask a little bit
            template          = reshape(template,[dm(1:3) K]);
            template(:,:,:,k) = imgaussfilt3(template(:,:,:,k),SmoothSig);
            template          = reshape(template,[prod(dm(1:3)) K]);
        end
    end    
end

% Renormalise
template = reshape(template,[dm(1:3) K]);
template = bsxfun(@rdivide,template,sum(template,4) + eps);
%==========================================================================

%==========================================================================
function val = get_border_volume(msk)

dm = size(msk);
dm = [dm 1];

if dm(3) ==1
    % 2D
    N       = 4;
    vals    = zeros([1 N]);   
    vals(1) = sum(sum(msk(:,1)));
    vals(2) = sum(sum(msk(:,end)));
    vals(3) = sum(sum(msk(1,:)));
    vals(4) = sum(sum(msk(end,:)));
else
    % 3D
    N       = 6;
    vals    = zeros([1 N]);   
    vals(1) = sum(sum(sum(msk(:,:,1))));
    vals(2) = sum(sum(sum(msk(:,:,end))));
    vals(3) = sum(sum(sum(msk(1,:,:))));
    vals(4) = sum(sum(sum(msk(end,:,:))));
    vals(5) = sum(sum(sum(msk(:,1,:))));
    vals(6) = sum(sum(sum(msk(:,end,:))));
end

val = sum(vals);

return
%==========================================================================