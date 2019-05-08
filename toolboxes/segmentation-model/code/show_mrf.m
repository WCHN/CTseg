function show_mrf(mrf)

figname = '(SPM) MRF';

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);  

G = mrf.G;

imagesc(G); axis off; colorbar;

drawnow;


%==========================================================================