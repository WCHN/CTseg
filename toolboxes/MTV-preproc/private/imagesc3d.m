function imagesc3d(img,varargin)
% FORMAT imagesc3d(img,varargint)
% img      - A 2D or 3D image
% varargin - Options for MATLAB's imagesc() function
%
% Same as imagesc() if 2D, otherwise displays a montage of slices from a 3D
% image.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

dm        = size(img);
if numel(dm) == 2
    dm(3) = 1;
end
is_rbg    = dm(3) == 3 || (numel(dm) == 4 && dm(4) == 3);

if numel(dm) > 3 && ~is_rbg
    error('imagesc3d does not support >3D arrays, unless 4D array has three channels in fourth dimension')
elseif dm(3) > 1 && dm(3) ~= 3
    
    % Montage parameters: this is the spacing for picking slices in the
    % z-axis
    Spacing = 8;
    
    % Set up montage
    z  = 1:Spacing:dm(3);
    N  = numel(z);    
    nc = floor(sqrt(N));
    nr = ceil(N/nc);  
    
    % Create montage
    if is_rbg
        mtg = zeros([nr*dm(1) nc*dm(2) 3],'single');
    else
        mtg = zeros([nr*dm(1) nc*dm(2)],'single');
    end
    cnt = 1;
    for r=1:nr
        for c=1:nc
            if cnt > numel(z)
                break
            end
            
            if is_rbg
                mtg(1 + (r - 1)*dm(1):r*dm(1),1 + (c - 1)*dm(2):c*dm(2),:) = squeeze(img(:,:,z(cnt),:));
            else
                mtg(1 + (r - 1)*dm(1):r*dm(1),1 + (c - 1)*dm(2):c*dm(2)) = img(:,:,z(cnt));
            end
            
            cnt = cnt + 1;
        end
    end   
    img = mtg;
    clear mtg
end

if is_rbg
    img = permute(img,[2 1 3]);
else
    img = img';
end

% Show image with imagesc()
imagesc(img,varargin{:});
%==========================================================================