function Nii = alloc_aux_vars(Nii,do_readwrite,dm,mat,use_projmat,p)
% Allocate MTV model auxiliary variables
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C            = numel(Nii.x);
dir_tmp      = p.Results.TemporaryDirectory;
EstimateBias = p.Results.EstimateBias;
ApplyBias    = p.Results.ApplyBias;
IsMPM        = p.Results.IsMPM;

if do_readwrite
    % Read/write temporary variables from disk (stored as NIfTIs)
    Nii.y = nifti;
    Nii.u = nifti;
    Nii.w = nifti;
    Nii.H = nifti;
else
    % Keep temporary variables in memory
    Nii.y = struct;
    Nii.u = struct;
    Nii.w = struct;
    Nii.H = struct;   
end
Nii.b = cell(1,C);

bf_list = cell(1,2);
for c=1:C
    if do_readwrite
        fname_y = fullfile(dir_tmp,['y' num2str(c) '.nii']);
        fname_u = fullfile(dir_tmp,['u' num2str(c) '.nii']); 
        fname_w = fullfile(dir_tmp,['w' num2str(c) '.nii']);
                       
        create_nii(fname_y,zeros(dm,      'single'),mat,[spm_type('float32') spm_platform('bigend')],'y');
        create_nii(fname_u,zeros([dm 3 2],'single'),mat,[spm_type('float32') spm_platform('bigend')],'u');
        create_nii(fname_w,zeros([dm 3 2],'single'),mat,[spm_type('float32') spm_platform('bigend')],'w');
        
        Nii.y(c) = nifti(fname_y);
        Nii.u(c) = nifti(fname_u);
        Nii.w(c) = nifti(fname_w);
        
        if use_projmat        
            % Approximate Hessian when using A and At
            fname_H = fullfile(dir_tmp,['H' num2str(c) '.nii']);        
            create_nii(fname_H,zeros(dm,'single'),mat,[spm_type('float32') spm_platform('bigend')],'H');
            Nii.H(c) = nifti(fname_H);
        else
            Nii.H(c) = nifti;
        end          
        
        % Bias field parameters
        Nii.b{c} = nifti;
        for n=1:numel(Nii.x{c})
            
            fname_img = Nii.x{c}(n).dat.fname;

            if ApplyBias || EstimateBias
                
                fname_b = fullfile(dir_tmp,['bf' num2str(c) num2str(n) '.nii']);  
                                
                if ApplyBias
                    
                    if IsMPM
                        % MPM data, use the same bf over echoes
                        [bf,bf_list] = get_bf_mpm(fname_img,bf_list,fname_b);
                    else
                        bf           = get_bf_spm_preproc8(fname_img);
                    end
                                        
                elseif EstimateBias
                    bf = zeros(dm,'single');
                end
                
                create_nii(fname_b,bf,mat,[spm_type('float32') spm_platform('bigend')],'bf');
                Nii.b{c}(n) = nifti(fname_b);
            else
                Nii.b{c}(n) = nifti;
            end
        end
    else
        Nii.y(c).dat = zeros(dm,      'single');
        Nii.u(c).dat = zeros([dm 3 2],'single');
        Nii.w(c).dat = zeros([dm 3 2],'single');
        
        if use_projmat
            % Approximate Hessian when using A and At
            Nii.H(c).dat = zeros(dm,'single');
        else
            Nii.H(c).dat = 0;
        end
        
        % Bias field parameters
        for n=1:numel(Nii.x{c})
    
            fname_img = Nii.x{c}(n).dat.fname;
            
            if ApplyBias || EstimateBias
                
                if ApplyBias                    
                    bf = get_bf_spm_preproc8(fname_img);                    
                elseif EstimateBias
                    bf = zeros(dm,'single');
                end
                
                Nii.b{c}(n).dat = bf;
            else
                Nii.b{c}(n).dat = 0;
            end
        end
    end
end

return
%==========================================================================

%==========================================================================
function [bf,bf_list] = get_bf_mpm(fname_img,bf_list,fname_b)

[pth,nam]    = fileparts(fname_img);
ix           = strfind(nam,'MP');
nam_echo_ser = nam(1:14);
echoes       = spm_select('FPList',pth,['^' nam_echo_ser '.*\.nii$']);
nam_echo_ser = nam(ix:14);
for i=1:size(echoes,1)
    fname_echo1 = strtrim(echoes(i,:));
    num         = fname_echo1(end - 5:end - 4);
    
    if str2double(num) == 1
        break
    end
end

if any(strcmp(bf_list{1},nam_echo_ser))
    ix      = strcmp(bf_list{1},nam_echo_ser);
    fname_b = bf_list{2}{ix};
    Nii     = nifti(fname_b);
    bf      = single(Nii.dat(:,:,:));
else
    bf                  = get_bf_spm_preproc8(fname_echo1);
    bf_list{1}{end + 1} = nam_echo_ser;
    bf_list{2}{end + 1} = fname_b;
end
%==========================================================================

%==========================================================================
function bfc = get_bf_spm_preproc8(fname_img)

[~,nam,ext] = fileparts(fname_img);
dir_bf      = '.';
fname_bf    = fullfile(dir_bf,['BiasField_' nam ext]);

obj          = struct;
obj.bb       =  NaN(2,3);
obj.bb       = [-90 -126 -72; 90 90 108];
obj.vox      =  NaN;
obj.cleanup  =  1;
obj.mrf      =  2;
obj.affreg   = 'mni';
obj.reg      = [0 0.001 0.5 0.05 0.2]*0.1;
obj.fwhm     = 1;
obj.samp     = 3;
obj.biasreg  = 0.001*(1/5);
obj.biasfwhm = 60;

tpmname      = fullfile(spm('dir'),'tpm','TPM.nii');
obj.lkp      = [1 1 2 2 3 3 4 4 5 5 5 6 6];
obj.tpm      = spm_load_priors8(tpmname);
obj.image    = spm_vol(char(fname_img));


M = obj.image(1).mat;
c = (obj.image(1).dim+1)/2;
obj.image(1).mat(1:3,4) = -M(1:3,1:3)*c(:);
[Affine1,ll1]    = spm_maff8(obj.image(1),8,(0+1)*16,obj.tpm,[],obj.affreg); % Closer to rigid
Affine1          = Affine1*(obj.image(1).mat/M);
obj.image(1).mat = M;

% Run using the origin from the header
[Affine2,ll2]    = spm_maff8(obj.image(1),8,(0+1)*16,obj.tpm,[],obj.affreg); % Closer to rigid

% Pick the result with the best fit
if ll1>ll2
    obj.Affine  = Affine1;
else
    obj.Affine  = Affine2;
end

% Initial affine registration.
obj.Affine     = spm_maff8(obj.image(1),obj.samp*2,(obj.fwhm+1)*16,obj.tpm, obj.Affine, obj.affreg); % Closer to rigid
obj.Affine     = spm_maff8(obj.image(1),obj.samp*2, obj.fwhm,      obj.tpm, obj.Affine, obj.affreg);

% Run the actual segmentation
res = spm_preproc8(obj);

% Final iteration, so write out the required data.
required    = false(1,2);
required(1) = true;
spm_preproc_write8(res,false(max(obj.lkp),4),required,false(1,2),obj.mrf,obj.cleanup,obj.bb,obj.vox,dir_bf);

% Get bias field
Nii                 = nifti(fname_bf);
bfc                 = single(-log(Nii.dat(:,:,:)));
bfc(~isfinite(bfc)) = 0;

delete(fname_bf);
%==========================================================================