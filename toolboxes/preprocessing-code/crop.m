function Nii = crop(Nii,rem_neck)
if nargin < 2, rem_neck = true; end

fprintf('Cropping...')
N        = numel(Nii{1});
prefix   = 'cr';
for n=1:N
    f = Nii{1}(n).dat.fname;
    
    atlas_crop(f,prefix,rem_neck);
    
    [pth,nam,ext] = fileparts(f);
    nf            = fullfile(pth,[prefix nam ext]);
    
    delete(f);
    Nii{1}(n) = nifti(nf);
end
fprintf('done!\n')
%==========================================================================

%==========================================================================
function [Affine,bb] = atlas_crop(P,prefix,rem_neck)
% Removes air outside of head
% FORMAT [Affine,bb] = atlas_crop(P,Affine,prefix,rem_neck)
% P        - Path to NIfTI file
% prefix   - File prefix (if empty -> overwrites) ['']
% rem_neck - Remove neck/spine [false]
% bb - Computed bounding box
%
% This function rigidly registers the SPM atlas to an image and then
% removes image data outside of the head.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<2, prefix   = 'cr'; end
if nargin<3, rem_neck = 0;  end

% Locate TPM.nii in SPM
pth_tpm = fullfile(spm('dir'),'tpm','TPM.nii,');

Vin  = spm_vol(P);
Vtpm = spm_vol(pth_tpm);

mat    = Vin.mat;
mattpm = Vtpm.mat;

tpm = spm_load_priors8(Vtpm);    

V = spm_vol(P);

M               = V(1).mat;
c               = (V(1).dim+1)/2;
V(1).mat(1:3,4) = -M(1:3,1:3)*c(:);
[Affine1,ll1]   = spm_maff8(V(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid
Affine1         = Affine1*(V(1).mat/M);

% Run using the origin from the header
V(1).mat      = M;
[Affine2,ll2] = spm_maff8(V(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid

% Pick the result with the best fit
if ll1>ll2, Affine  = Affine1; else Affine  = Affine2; end

Affine = spm_maff8(P,4,32,tpm,Affine,'mni');
Affine = spm_maff8(P,4,1,tpm,Affine,'mni');

% Voxel locations in TPM.nii
Ltpm1 = [120 72.2 37.3 1]'; Ltpm2 = [120 72.2 75.9 1]';
Rtpm1 = [3  72.2 37.3 1]'; Rtpm2 = [3  72.2 75.9 1]';

Stpm1 = [58.6 42.6 119 1]'; Stpm2 = [58.60 99.4 119 1]';
if rem_neck==2
    Itpm1 = [61 42 31   1]'; Itpm2 = [61 107 29 1]';   
elseif rem_neck==1
    Itpm1 = [58.6 39.4 2.5   1]'; Itpm2 = [58.60 99.4 2.5  1]';    
else
    Itpm1 = [58.6 39.4 -200  1]'; Itpm2 = [58.60 99.4 -200 1]';
end
Atpm1 = [58.6 144 28.4 1]'; Atpm2 = [58.60 144 82.3 1]';
Ptpm1 = [58.6 3.5  28.4 1]'; Ptpm2 = [58.60 3.5 82.3 1]'; 

% Voxel locations in input
T  = mat\(Affine\mattpm);
L1 = T*Ltpm1; L2 = T*Ltpm2;
R1 = T*Rtpm1; R2 = T*Rtpm2;
U1 = T*Stpm1; U2 = T*Stpm2;
D1 = T*Itpm1; D2 = T*Itpm2;
A1 = T*Atpm1; A2 = T*Atpm2;
P1 = T*Ptpm1; P2 = T*Ptpm2;

% Bounding-box
bb = zeros(2,3);
for i=1:3
    X       = [L1(i) R1(i) U1(i) D1(i) A1(i) P1(i)...
               L2(i) R2(i) U2(i) D2(i) A2(i) P2(i)];
    bb(1,i) = max(X);
    bb(2,i) = min(X);
end

% Do cropping
subvol(Vin,bb,prefix);      
%==========================================================================