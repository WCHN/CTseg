function show_reg(Z,Template,figname)
              
if nargin < 3
    figname = '(SPM) Registration';
end

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   

spm_gmm_lib('plot','showcatimg',{Z,Template},{'Z','Template'},{'bg','st','gm','bas','wm','wmh','csf','ven','cer','spn','skl'});
%==========================================================================