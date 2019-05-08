function spm_check_path(varargin)
% Check that SPM is on the MATLAB path. Can also check to see if some other
% tools are avaiable.
%
% EXAMPLE USAGE
%
% Call as spm_check_path('Shoot','Longitudinal','pull') to check if these
% toolboxes and/or functions are available. If you just want to check if,
% for example, the Shoot toolbox is available, just do spm_check_path('Shoot').
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
   
% Check that SPM is on the MATLAB path
if ~(exist('spm','file') == 2)
    error('SPM is not on the MATLAB path, get the latest SPM version from: https://www.fil.ion.ucl.ac.uk/spm/software/ and add it to the MATLAB path'); 
end

if ~isempty(varargin)
    
    if any(strcmpi(varargin,'Shoot'))
        % Check that the Shoot toolbox is on the MATLAB path
        
        try
            spm_shoot3d;
        catch e
            if strcmp(e.message,'Undefined function or variable ''spm_shoot3d''.')
                error('Add the Shoot toolbox, from /spm/toolbox/Shoot, to the MATLAB path')
            end
        end
    end

    if any(strcmpi(varargin,'Longitudinal'))
        % Check that the Longitudinal toolbox is on the MATLAB path
        
        try
            spm_groupwise_ls;
        catch e
            if strcmp(e.message,'Undefined function or variable ''spm_groupwise_ls''.')
                error('Add the Longitudinal toolbox, from /spm/toolbox/Longitudinal, to the MATLAB path')
            end
        end
    end
    
    if any(strcmpi(varargin,'pull'))
        % Check that spm_diffeo('pull') is available
        
        try
            spm_diffeo('pull');
        catch e
            if strcmp(e.message,'Option not recognised.')
                error('The function spm_diffeo(''pull'') is not available, update to the latest SPM version from: https://www.fil.ion.ucl.ac.uk/spm/software/')
            end
        end
    end
end
%==========================================================================