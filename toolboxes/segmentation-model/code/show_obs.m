function show_obs(img,dm,title_nam,figname)

if nargin < 4
    figname = '(SPM) Observed data';
end

C = size(img,2);
if nargin < 3, title_nam = cell(1,C); end

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   

img = reshape(img,[dm C]);

if dm(3)==1
    % 2d 
    %----------------------------------------------------------------------
    
    for c=1:C    
        subplot(1,C,c)
        
        slice = img(:,:,1,c);
        slice = permute(slice,[2 1 3]);
        imagesc(slice); axis off xy; colorbar;  
        
        if ~isempty(title_nam{c})
            title(title_nam{c})
        end
    end
    
    colormap(gray);
else
    % 3d
    %----------------------------------------------------------------------
    
    for c=1:C  
        subplot(C,3,3*(c - 1) + 1)
        
        slice = img(:,:,floor(dm(3)/2) + 1,c);
        slice = permute(slice,[2 1 3]);
        imagesc(slice); axis off xy; colorbar;  

        if ~isempty(title_nam{c})
            title(title_nam{c})
        end
        
        subplot(C,3,3*(c - 1) + 2)
        
        slice = permute(img(:,floor(dm(2)/2) + 1,:,c),[3 1 2 4]);  
        imagesc(slice); axis off xy; colorbar;  

        if ~isempty(title_nam{c})
            title(title_nam{c})
        end
        
        subplot(C,3,3*(c - 1) + 3)
        
        slice = permute(img(floor(dm(1)/2) + 1,:,:,c),[2 3 1 4]);
        slice = permute(slice,[2 1 3]);
        imagesc(slice); axis off xy; colorbar;  
        
        if ~isempty(title_nam{c})
            title(title_nam{c})
        end
        
        colormap(gray);          
    end    
end

drawnow
%==========================================================================