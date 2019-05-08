function make_animations(operation,opt,iter)

figHandles = findobj('Type', 'figure');

if strcmp(operation,'savefigs')    
    for i=1:numel(figHandles)
        save_figs(figHandles(i).Name,opt.dir_animations,iter);
    end
elseif strcmp(operation,'clear')  
    for i=1:numel(figHandles)
        clf(figHandles(i));
    end    
elseif strcmp(operation,'make')
    write2gif(figHandles,opt.dir_animations,iter);
else
    error('Unknown command!')
end
%==========================================================================

%==========================================================================
function save_figs(figname,dir_animations,n)
f = findobj('Type', 'Figure', 'Name', figname);
if ~isempty(f)
    set(0, 'CurrentFigure', f);      
    filename = fullfile(dir_animations,[figname '_' num2str(n) '.png']);
    saveas(gcf,filename)
end
%==========================================================================

%==========================================================================
function write2gif(figHandles,dir_animations,iter)
for i=1:numel(figHandles)
    figname = figHandles(i).Name;
    gifname = fullfile(dir_animations,[figname '.gif']);
    
    first_im = true;
    for n=1:iter
        filename = fullfile(dir_animations,[figname '_' num2str(n) '.png']);
        
        if ~isfile(filename), continue; end
        
        im         = imread(filename);
        [imind,cm] = rgb2ind(im,256); 

        if first_im
            imwrite(imind,cm,gifname,'gif', 'Loopcount',inf); 
            first_im = false;
        else 
            imwrite(imind,cm,gifname,'gif','WriteMode','append'); 
        end 
    end    
end
%==========================================================================      