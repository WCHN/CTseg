function Z = clean_brain(Z,model,opt,y,Verbose)
% Clean-up resulting (soft) segmentations using a series of educated
% heuristics.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 5, Verbose = false; end

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

template = single(model.template.nii.dat(:,:,:,:));
template = spm_matcomp('softmax',template);  

% Build a 3x3x3 seperable smoothing kernel (krn)
krn.x=[0.75 1 0.75];
krn.y=[0.75 1 0.75];
krn.z=[0.75 1 0.75];
sm=sum(kron(kron(krn.z,krn.y),krn.x))^(1/3);
krn.x=krn.x/sm; krn.y=krn.y/sm; krn.z=krn.z/sm;
clear sm

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
        % Close holes
        msk = binmorph('close',msk,it_dil_er);

        if Verbose
            % Some verbose        
            imagesc3d(msk); axis off; drawnow
        end  
        
        msk = imfill(msk,'holes');
        
        if Verbose
            % Some verbose        
            imagesc3d(msk); axis off; drawnow
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
    
    if Verbose
        % Show template before clean
        figure(FigNum + 1)
        nr = floor(sqrt(K));
        nc = ceil(K/nr);      
        for k=1:K   
            subplot(nr,nc,k)
            imagesc3d(template(:,:,:,k)); axis off; drawnow
        end
    end
    
    % Update template with mask information
    template = adjust_Z_brain(template,msk,BrainClass,LesionClass,Tiny,SmoothSig);
    
    if Verbose
        % Show template after clean
        figure(FigNum + 2)
        nr = floor(sqrt(K));
        nc = ceil(K/nr);      
        for k=1:K   
            subplot(nr,nc,k)
            imagesc3d(template(:,:,:,k)); axis off; drawnow
        end
    end        

    if Verbose
        % Show Z before clean
        figure(FigNum + 3)
        nr = floor(sqrt(K));
        nc = ceil(K/nr);      
        for k=1:K   
            subplot(nr,nc,k)
            imagesc3d(Z(:,:,:,k)); axis off; drawnow
        end
    end
    
    % Warp mask to subject space
    msk                 = single(msk);
    msk                 = imgaussfilt3(msk,SmoothSig);
    msk                 = spm_diffeo('bsplins',msk,y,[4 4 4  0 0 0]);
    msk(~isfinite(msk)) = 0;    
    msk                 = msk > 0;
    
    % Fill holes
    msk = imfill(msk,'holes');
    
    % Erode
    msk = binmorph('erode',msk,it_dil_er);
    
    if Verbose
        % Some warped mask    
        figure(FigNum + 4)
        imagesc3d(msk); axis off; drawnow
    end          
    
    % Update Z with mask information    
    Z = adjust_Z_brain(Z,msk,BrainClass,LesionClass,Tiny,SmoothSig);

    if Verbose
        % Show Z after clean
        figure(FigNum + 5)
        nr = floor(sqrt(K));
        nc = ceil(K/nr);      
        for k=1:K   
            subplot(nr,nc,k)
            imagesc3d(Z(:,:,:,k)); axis off; drawnow
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

    % Warp mask to subject space
    msk                 = spm_diffeo('bsplins',single(msk),y,[0 0 0  0 0 0]);
    msk(~isfinite(msk)) = 0;

    % Update template and Z with mask information    
    Z = adjust_Z_air(Z,msk,BrainClass,LesionClass,Tiny,SmoothSig);

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
%==========================================================================

%==========================================================================
function template = adjust_Z_brain(template,msk,BrainClasses,LesionClass,Tiny,SmoothSig)
dm       = size(template);
K        = dm(4);
msk      = reshape(msk,[prod(dm(1:3)) 1]);
template = reshape(template,[prod(dm(1:3)) K]);
for k=1:K    
    if ismember(k,BrainClasses) || ismember(k,LesionClass)
        template(~msk,k) = Tiny;
    end    
end

% Renormalise
template = reshape(template,[dm(1:3) K]);
template = bsxfun(@rdivide,template,sum(template,4) + eps);
%==========================================================================

%==========================================================================
function template = adjust_Z_air(template,msk,AirClass,Tiny,SmoothSig)
dm       = size(template);
K        = dm(4);
msk      = reshape(msk,[prod(dm(1:3)) 1]);
template = reshape(template,[prod(dm(1:3)) K]);
for k=1:K    
    if ismember(k,AirClass)
        template(msk,k) = 1 - Tiny;
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

%==========================================================================
function vol = binmorph(action,bim,n)

if nargin < 3, n = 1; end % Iterations
    
vol = uint8(bim);

kx = [1 1 1];
ky = [1 1 1];
kz = [1 1 1];

order = sum(kx(:) ~= 0)*sum(ky(:) ~= 0);

switch lower(action)
    
	case 'dilate'
	
        sz = size(vol);
        vol2 = zeros(sz(1)+(2*n),sz(2)+(2*n),sz(3)+(2*n),'uint8');
        vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n) = vol;
        for i = 1:n
            spm_conv_vol(vol2,vol2,kx,ky,kz,-[1 1 1]);
            vol2 = uint8(vol2~=0);
            
%             imagesc3d(vol2); axis off; drawnow
        end
        vol = vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n);

	case 'erode'

        for i = 1:n
            spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
            vol = uint8(vol>=order);
            
%             imagesc3d(vol); axis off; drawnow
        end
        
	case 'close'
        
        sz = size(vol);
        vol2 = zeros(sz(1)+(2*n),sz(2)+(2*n),sz(3)+(2*n),'uint8');
        vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n) = vol;
        for i = 1:n
            spm_conv_vol(vol2,vol2,kx,ky,kz,-[1 1 1]);
            vol2 = uint8(vol2~=0);
            
%             imagesc3d(vol2); axis off; drawnow
        end                
        
        for i = 1:n
            spm_conv_vol(vol2,vol2,kx,ky,kz,-[1 1 1]);
            vol2 = uint8(vol2>=order);
            
%             imagesc3d(vol2); axis off; drawnow
        end
        vol = vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n);
        clear vol2                        
end

vol = logical(vol);
%==========================================================================