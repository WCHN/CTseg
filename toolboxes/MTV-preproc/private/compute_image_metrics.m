function [avg_psnr,avg_ssim] = compute_image_metrics(Nii,Nii_ref)
% Compute image quality metrics
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimagin

C     = numel(Nii_ref);
psnrs = zeros(1,C);
ssims = zeros(1,C);

for c=1:C % Loop over channels
    if iscell(Nii(c))
        img = get_nii(Nii{c});
    else
        img = get_nii(Nii(c));
    end
    
    ref = get_nii(Nii_ref(c));   
    
    psnrs(c) = get_psnr(img(:),ref(:));
    ssims(c) = ssim(img,ref);    
end

avg_psnr = mean(psnrs);
avg_ssim = mean(ssims);
%==========================================================================

%==========================================================================
function val = get_psnr(hat,ref)
% Calculates the peak signal-to-noise ratio (PSNR) for the image in array 
% hat, with the image in array ref as the reference.
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

hat = double(hat(:)); 
ref = double(ref(:));

peak = max(hat);
peak = max(peak,max(ref));

RMSE = sqrt(mean((ref - hat).^2)); 

val = 20*log10(peak/RMSE);
%==========================================================================