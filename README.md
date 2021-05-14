# CTseg: Brain CT segmentation, normalisation, skull-stripping and total brain/intracranial volume computation

## Overview

<img style="float: right;" src="https://github.com/WCHN/CTseg/blob/master/demo.png" width="80%" height="80%">

This is an algorithm for segmenting and spatially normalising computed tomography (CT) brain scans. The model is an extension of the popular unified segmentation routine (part of the SPM12 software) with: improved registration, priors on the Gaussian mixture model parameters, an atlas learned from both MRIs and CTs (with more classes). These improvements leads to a more **robust** segmentation routine that can better handle image with lots of noise and/or large morphological variability (see figure above). The algorithm can produce native|warped|modulated space segmentations of:

1. Gray matter (GM)
2. White matter (WM)
3. Cerebrospinal fluid (CSF)
4. Skull (BONE)
5. Soft tissue (ST)
6. Background (BG)

The implementation is done in MATLAB and depends on the SPM12 package (and its MB toolbox). If you find the code useful, please consider citing the publications in the *References* section.

## Further details

The input to CTseg should be provided as NIfTI files (```.nii```). The resulting tissue segmentations are in the same format as the output of the SPM12 segmentation routine (```c*```, ```wc*```, ```mwc*```). The normalised segmentations (```wc*```, ```mwc*```) are in MNI space.

A **skull-stripped** version of the input image is produced by default (prefixed ```ss_``` to the original filename). **Total brain volume** (TBV) and **intercranial volume** (TIV) are also computed by the algorithm and returned as the second argument of the CTseg function. Note that both of these routines uses only the GM, WM and CSF segmentations of the algorithm. The skull-stripped volume will therefore not include the meninges, the sinuses or any calcifications; the TIV might therefore also be slighly underestimated.

For converting **DICOM** CT to NIfTI, we recommend using SPM12's ```spm_dicom_convert```. This DICOM converter can deal with the fact that many CT images are often acquired with variable slice thickness. If this is not accounted for when reconstructing the NIfTI file, the head shape can be deformed.

## Dependencies

The algorithm is developed using MATLAB and relies on external functionality from the SPM12 software. The following are therefore required downloads and need to be placed on the MATLAB search path (using `addpath`):
* **SPM12:** Download from https://www.fil.ion.ucl.ac.uk/spm/software/download/.
* **Shoot toolbox:** The Shoot folder from the toolbox directory of SPM12.
* **Longitudinal toolbox:** The Longitudinal folder from the toolbox directory of SPM12.
* **Multi-Brain toolbox:** Download/clone https://github.com/WTCN-computational-anatomy-group/mb and follow the *Installation instructions*.

## Example use case

Below are two MATLAB snippets. The first takes as input a CT image (as ```*.nii```) and produces native space GM, WM, CSF tissue segmentations (```c[1-3]*.nii```), as well as template space (MNI) non-modulated (```wc[1-3]*.nii```) and modulated (```mwc[1-3]*.nii```) ones. The forward deformation that warps the atlas to the native space CT is also written to disk (as ```y_*.nii```). The second snippet uses the forward deformation to warp: (1) the CT image to the template space; and (2), the template to the space of the CT image. Note that the template is here softmaxed to make it probabilistic. Results are written to ```dir_out```.

### 1. CT segmentation and normalisation

``` matlab
% Path to a CT image
pth_ct = 'CT.nii';  

% Output directory (if empty, uses same directory as input data)
dir_out = ''; 

% Write all types of segmentations
tc = true;  

% Write forward deformation to disk?
def = true;  

% Correct orientation matrix?
correct_header = true;  

% Do skull-stripping?
ss = true;

% Template space voxel size
vox = 1.0;

% Run segmentation+normalisation
[res,vol] = spm_CTseg(pth_ct, dir_out, tc, def, correct_header, ss, vox)
% res: a struct with paths to result niftis
% vol: a struct containing TBV and TIV
```

### 2. Warping with the generated deformation

``` matlab
% Path to tissue template (this should be located in the CTseg folder)
pth_mu = 'mu_CTseg.nii';

% Path to forward deformation (in dir_out)
pth_y = 'y_*.nii';

% Use forward deformation (pth_y) to warp CT image (pth_ct) to 
% template space (pth_mu)
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def     = {pth_y};
matlabbatch{1}.spm.util.defs.comp{1}.inv.space           = {pth_mu};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {pth_ct};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {dir_out};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = 'w';
spm_jobman('run',matlabbatch);

% Use forward deformation (pth_y) to warp template (pth_mu) to native 
% CT image space (pth_ct) 
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.comp{1}.def         = {pth_y};
matlabbatch{1}.spm.util.defs.comp{1}.space               = {pth_ct};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {pth_mu};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {dir_out};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = 'w';
spm_jobman('run',matlabbatch);
```

## Troubleshooting

* **Error using spm_gmm_lib>loop:** If you get the error message "At least one of Prop, LogProp or Dir must be provided.", ensure the option `correct_header` is set to true.

* **Out of memory error:** Some CT scans can have quite large file size, as they might have large coverage and small voxels (submillimetric), which could lead to memory issues when running CTseg. Two solutions to this problem is to either subsample the CT image, or find a computer with more RAM...

* **Segmentation results not as expected:** The model file could not have been found. Make sure that the files ```prior_CTseg.mat``` and ```mu_CTseg.nii``` exist in the directory of CTseg. They are in the ```model.zip``` file, which should get automatically downloaded and unzipped when the code is executed for the first time.

* **Error related to spm_diffeo:** This code uses a recent version of SPM12; therefore, if your SPM12 version is quite old, the function ```spm_diffeo``` might break. Updating to the latest version of SPM12 will resolve this issue.

## Improved runtime (Linux and Mac)

For a faster algorithm, consider compiling SPM with OpenMP support. Just go to the *src* folder of SPM and do:
``` bash
make distclean
make USE_OPENMP=1 && make install
```

## References

``` latex
@inproceedings{brudfors2020flexible,
  title={Flexible Bayesian Modelling for Nonlinear Image Registration},
  author={Brudfors, Mikael and Balbastre, Ya{\"e}l and Flandin, Guillaume and Nachev, Parashkev and Ashburner, John},
  booktitle={International Conference on Medical Image Computing and Computer-Assisted Intervention},
  pages={253--263},
  year={2020},
  organization={Springer}
}

@phdthesis{brudfors2020generative,
  title={Generative Models for Preprocessing of Hospital Brain Scans},
  author={Brudfors, Mikael},
  year={2020},
  school={UCL (University College London)}
}
```
## License

CTseg is free but copyright software, distributed under the terms of the GNU General Public Licence as published by the Free Software Foundation (either version 2, or at your option, any later version).
