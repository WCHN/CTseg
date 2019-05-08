function Nii_H = approx_hessian(Nii_H,dat)
% Compute approximation to the diagonal of the Hessian 
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C = numel(dat);

for c=1:C        
    H        = At(A(ones(dat(c).dm,'single'),dat(c)),dat(c));   
%     H        = max(H(:))*ones(size(H),'single'); % Infinity norm
    Nii_H(c) = put_nii(Nii_H(c),H);
end    
%==========================================================================