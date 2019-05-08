function deal_figs(model)
% FORMAT deal_figs(model)
%
% Distribute all plotting windows in a nice way.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

screen_size = get(0,'ScreenSize');
w           = screen_size(3);
h           = screen_size(4);
top_offset  = 90;
dm          = model.template.nii.dat.dim;
is3d        = dm(3) > 1;

figname = '(SPM) Tissue proportions';

f = findobj('Type', 'Figure', 'Name', figname);
if ~isempty(f) && ~strcmpi(get(f, 'WindowStyle'), 'docked')
    set(f, 'Position', [w - w/6, 0, w/6, h/2]) % [x y w h]
end

figname = '(SPM) Template';

f = findobj('Type', 'Figure', 'Name', figname);
if ~isempty(f) && ~strcmpi(get(f, 'WindowStyle'), 'docked')
    if is3d
        set(f, 'Position', [w - (w/6 + w/2), 0, w/2, h/2]) % [x y w h]    
    else
        set(f, 'Position', [w - (w/6 + w/3), 0, w/3, h/2])   % [x y w h]
    end
end

figname = '(SPM) Model lower bound';

f = findobj('Type', 'Figure', 'Name', figname);
if ~isempty(f) && ~strcmpi(get(f, 'WindowStyle'), 'docked')
    set(f, 'Position', [0, h - (h/2 + top_offset), w/2, h/2]) % [x y w h]
end

figname = '(SPM) GaussPrior';

figHandles = findobj('Type', 'figure');
for i=1:numel(figHandles)
    nam = figHandles(i).Name;
    
    if numel(nam) >= 16 && strcmp(figname,nam(1:16)) ...
            && ~strcmpi(get(figHandles(i), 'WindowStyle'), 'docked')
        set(figHandles(i), 'Position', [w - w/2, h - (h/2 + top_offset), w/2, h/2]) % [x y w h]
    end
end

figname = '(SPM) Bias field and initial velocity';

f = findobj('Type', 'Figure', 'Name', figname);
if ~isempty(f) && ~strcmpi(get(f, 'WindowStyle'), 'docked')
    set(f, 'Position', [w - w/3, 0, w/3, h]) % [x y w h]
end

figname = '(SPM) Subject Lower Bound';

f = findobj('Type', 'Figure', 'Name', figname);
if ~isempty(f) && ~strcmpi(get(f, 'WindowStyle'), 'docked')
    set(f, 'Position', [w/2, h - (h/2 + top_offset), w/2, h/2]) % [x y w h]    
end

figname = '(SPM) Plot GMM';

f = findobj('Type', 'Figure', 'Name', figname);
if ~isempty(f) && ~strcmpi(get(f, 'WindowStyle'), 'docked')
    set(f, 'Position', [0, h - (h/2 + top_offset), w/2, h/2]) % [x y w h]
end

figname = '(SPM) Observed, template and responsibilities';

f = findobj('Type', 'Figure', 'Name', figname);
if ~isempty(f) && ~strcmpi(get(f, 'WindowStyle'), 'docked')
    if is3d
        set(f, 'Position', [0, 0, w/1.5, h/2]) % [x y w h]
    else
        set(f, 'Position', [0, 0, w/1.5, h/3]) % [x y w h]
    end
end

figname = '(SPM) MRF';

f = findobj('Type', 'Figure', 'Name', figname);
if ~isempty(f) && ~strcmpi(get(f, 'WindowStyle'), 'docked')
    set(f, 'Position', [w - w/5, h - (h/5 + top_offset), w/5, h/5]) % [x y w h]   
end

drawnow;
%==========================================================================