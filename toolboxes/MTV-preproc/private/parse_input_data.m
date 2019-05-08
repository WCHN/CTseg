function [Nii,C,is3d] = parse_input_data(Nii,InputImages,use_projmat)
% Parse input data and return Nii cell array
%
% InputImages  - [1 x C cell array] or empty (let user select)
% method       - String either 'superres' or 'denoise'
% Nii          - [1 x C cell array] Output array with observed data
% C            - Number of image channels
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

OrdinalNumbers = {'1st', '2nd', '3rd', '4th', '5th', '6th', '7th', '8th', '9th', '10th', '11th', ...
                  '12th', '13th', '14th', '15th', '16th', '17th', '18th', '19th', '20th', '21st', ...
                  '22nd', '23rd', '24th', '25th', '26th', '27th', '28th', '29th', '30th', '31st', '32nd'};
              
Nii.x = {};
if isempty(InputImages)
    % Empty input, let user select images from disk using spm_select
    C = input('Please specify number of image channels: ');        
    
    for c=1:C
        fprintf('Select all images of the %s channel.\n',OrdinalNumbers{c});
        Nii.x{c} = nifti(spm_select(Inf,'nifti','Select image'));    
    end
else
    isNii = true;
    C     = numel(InputImages);
    for c=1:C
        if iscell(InputImages)
            if ~isa(InputImages{c},'nifti')
                Nii.x{c} = nifti(InputImages{c}); 
                isNii    = false;
            end
        else
            if ~isa(InputImages(c),'nifti')
                Nii.x{c} = nifti(InputImages(c)); 
                isNii    = false;
            end
        end
    end        
    
    if isNii
        if isa(InputImages(c),'nifti')
            for c=1:C
                Nii.x{c} = InputImages(c);
            end
        else
            Nii.x = InputImages;
        end
    end
end

% Sanity check input (for denoising)
dm0 = Nii.x{1}(1).dat.dim;
for c=1:C
    N = numel(Nii.x{c});
    for n=1:N
        dm = Nii.x{c}(n).dat.dim; 
        if ~use_projmat && c > 1 && (~isequal(dm,dm0) || ~isequal(dm,dm0))
            error('Images are not all the same size!')
        end
    end
end
dm0  = [dm0 1];
is3d = dm0(3) > 1;
%==========================================================================