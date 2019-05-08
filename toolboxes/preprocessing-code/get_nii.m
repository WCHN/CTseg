function Nii = get_nii(in,DirOut)
if isa(in,'nifti')
    Nii = {in};
elseif isfolder(in)
    files = spm_select('FPListRec',in,'^.*\.dcm$');
    if ~isempty(files)    
        if ~(exist(DirOut,'dir') == 7)  
            mkdir(DirOut);  
        end
        hdr = spm_dicom_headers(files);
        out = spm_dicom_convert(hdr,'all','flat',spm_get_defaults('images.format'),DirOut);
        Nii = {nifti(out.files{1})};
    else
        error('No DICOMs in folder!')
    end        
else        
    Nii = {nifti(in)};
end
%==========================================================================
