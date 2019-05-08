function dat = init_load_a_der(dat,opt)
if opt.template.load_a_der                    
    S0 = numel(dat);
    for s=1:S0           
        if ~isfield(dat{s},'template')
            dat{s}.template = struct;
        end
        
        if ~isfield(dat{s}.template,'pth_gr') && ~isfield(dat{s}.template,'pth_H')
            dat{s}.template.pth_gr = fullfile(opt.dir_a_der,['gr-' num2str(s) '.nii']);
            dat{s}.template.pth_H  = fullfile(opt.dir_a_der,['H-' num2str(s) '.nii']);               
        end
        dat{s}.template.bb = zeros(3,2); % For saving bounding-box        
    end        
end
%==========================================================================  