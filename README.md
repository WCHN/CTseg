# CTseg: Brain CT segmentation and normalisation

<img style="float: right;" src="https://github.com/WCHN/CTseg/blob/master/demo.png" width="80%" height="80%">

This is a MATLAB implementation of a model for segmenting and spatially normalising computed tomography (CT) brain scans. The model is an extension of the popular unified segmentation routine (part of the SPM12 software) with: improved registration, priors on the Gaussian mixture model parameters, an atlas learned from both MRIs and CTs (with more classes). These improvements leads to a more **robust** segmentation routine that can better handle image with lots of noise and/or large morphological variability (see figure above). The algorithm can produce native|warped|modulated space segmentations of:

1. Gray matter (GM)
2. White matter (WM)
3. Cerebrospinal fluid (CSF)
4. Dura (DUR)
5. Skull, calcifications and hyper-intensities (BONE)
6. Soft tissue (ST)
7. Background (BG)

The input should be provided as NIfTI files (.nii). The resulting tissue segmentations are in the same format as the output of the SPM12 segmentation routine (```c*```, ```wc*```, ```mwc*```). Note that **the atlas is encoded in log-space**, not probabilisticly, a softmax operation is therefore needed to give voxel values that sum to one (see Example section). The normalised segmentations are in MNI space.

The code can be used either as: **(1)** an SPM12 extension, by adding it to the toolbox folder of SPM and using the batch interface (SPM -> Tools -> CT Segmentation); or **(2)** by interfacing with the code directly (example below).

The orientation matrix in the NIfTI header of CT scans could be messed up, this means that the atlas will not align with the image data. This is here fixed this by a preprocessing step. Note that this operation requires reslicing of the image data and therefore creates a copy of the original input data (as ```r*.nii```). Setting the ```correct_header``` option of CTseg to ```false``` disables this preprocessing step.

For converting DICOM CT to NIfTI, we recommend using SPM12's ```spm_dicom_convert```. This DICOM converter can deal with the fact that many CT images often are acquired with variable slice thickness. If this is not accounted for when reconstructing the NIfTI file, the head shape can be deformed.

If you find the code useful, please consider citing one of the publications in the *References* section.

## Dependencies

The algorithm requires that the following packages are on the MATLAB path:
* **SPM12:** Download from https://www.fil.ion.ucl.ac.uk/spm/software/spm12/.
* **Shoot toolbox:** Add Shoot folder from the toolbox directory of the SPM source code.
* **Longitudinal toolbox:** Add Longitudinal folder from the toolbox directory of the SPM source code.
* **Multi-Brain toolbox:** Download (or clone) from https://github.com/WTCN-computational-anatomy-group/diffeo-segment.

## Example use case

Below are two MATLAB snippets. The first takes as input a CT image (as ```*.nii```) and produces native space GM, WM, CSF tissue segmentations (```c[1-3]*.nii```), as well as template space (MNI) non-modulated (```wc[1-3]*.nii```) and modulated (```mwc[1-3]*.nii```) ones. The forward deformation that warps the atlas to the native space CT is also written to disk (as ```y_*.nii```). The second snippet uses the forward deformation to warp: (1) the CT image to the template space; and (2), the template to the space of the CT image. Note that the template is here softmaxed to make it probabalistic. Results are written to ```dir_out```.

### 1. CT segmentation and normalisation

```
% Path to a CT image
pth_ct = 'CT.nii';  

% Output directory (if empty, uses same directory as input data)
dir_out = ''; 

% What tissue classes to write to disk 
% (column: native, warped, modulated | row: GM, WM, CSF, DUR, BONE, ST, BG)
tc = [true(3, 3); false(4, 3)];  

% Write forward deformation to disk?
def = true;  

% Correct orientation matrix?
correct_header = false;  

% Run segmentation+normalisation
CTseg(pth_ct, dir_out, tc, def, correct_header)
```

### 2. Warping with the generated deformation

```
% Path to tissue template (this should be located in the CTseg folder)
pth_mu = 'mu_CTseg.nii';

% Path to forward deformation (in dir_out)
pth_y = 'y_*.nii';

% Write softmaxed template to dir_out
Nii = nifti(pth_mu);
mu  = Nii.dat();
mu  = spm_mb_shape('template_k1',mu);
mu  = exp(mu);
[pth,nam,ext] = fileparts(pth_mu);
nam           = ['softmax' nam(3:end)];
pth_mu        = fullfile(dir_out,[nam ext]);
fa            = file_array(pth_mu,size(mu),'float32',0);
Nmu           = nifti;
Nmu.dat       = fa;
Nmu.mat       = Nii.mat;
Nmu.mat0      = Nii.mat;
Nmu.descrip   = 'Template (softmax)';
create(Nmu);
Nmu.dat(:,:,:,:) = mu;

% Use forward deformation (pth_y) to warp CT image (pth_ct) to softmaxed 
% template space (pth_mu)
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {pth_y};
matlabbatch{1}.spm.util.defs.comp{1}.inv.space       = {pth_mu};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {pth_ct};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {dir_out};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = 'w';  % Output prefix
spm_jobman('run',matlabbatch); % Run job

% Use forward deformation (pth_y) to warp softmaxed template (pth_mu) to native 
% CT image space (pth_ct) 
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.comp{1}.def = {pth_y};
matlabbatch{1}.spm.util.defs.comp{1}.space       = {pth_ct};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {pth_mu};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {dir_out};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = 'w';  % Output prefix
spm_jobman('run',matlabbatch); % Run job
```

## Troubleshooting

* **Segmentation results not as expected:** The model file could not have been found. Make sure that the files ```prior_CTseg.mat``` and ```mu_CTseg.nii``` exist in the directory of CTseg. They are in the ```model.zip``` file, which should get automatically downloaded and unzipped when the code is executed for the first time.

* **Error related to spm_diffeo:** This code uses a recent version of SPM12; therefore, if your SPM12 version is quite old, the function ```spm_diffeo``` might break. Updating to the latest version of SPM12 will resolve this issue.

## Improved runtime (Linux and Mac)

For a faster algorithm, consider compiling SPM with OpenMP support. Just go to the *src* folder of SPM and do:
```
make distclean
make USE_OPENMP=1
make install
```

## References

1. Brudfors, M., Balbastre, Y., Flandin, G., Nachev, P., & Ashburner, J. (2020).
Flexible Bayesian Modelling for Nonlinear Image Registration. 
International Conference on Medical Image Computing and Computer Assisted Intervention.

2. Brudfors, M. (2020). 
Generative Models for Preprocessing of Hospital Brain Scans.
Doctoral dissertation, UCL (University College London).

## License

This software is released under the [GNU General Public License version 3](LICENSE) (GPL v3). As a result, you may copy, distribute and modify the software as long as you track changes/dates in source files. Any modifications to or software including (via compiler) GPL-licensed code must also be made available under the GPL along with build & install instructions.

[TL;DR: GPL v3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))
