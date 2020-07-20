# CTseg: Brain CT Segmentation and Normalization

<img style="float: right;" src="https://github.com/WCHN/CTseg/blob/master/demo.png" width="80%" height="80%">

This is a MATLAB implementation of a model for segmenting and spatially normalising computed tomography (CT) brain scans. The model is an extension of the popular unified segmentation routine (part of the SPM12 software) with: improved registration, priors on the Gaussian mixture model parameters, an atlas learned from both MRIs and CTs (with more classes), and more. These improvements leads to a more **robust** segmentation routine that can better handle image with lots of noise and/or large morphological variability (see figure above). The algorithm can produce native|warped|modulated space segmentations of:
1. Gray matter (GM)
2. White matter (WM)
3. Cerebrospinal fluid (CSF)
4. Dural venous sinuses (SIN)
5. Skull, calcifications and hyper-intensities (BONE)
6. Soft tissue (ST)
7. Background (BG)

The input should be provided as nifti files (.nii). The resulting tissue segmentations are in the same format as the output of the SPM12 segmentation routine (c*, wc*, mwc*).

The code can be used either as: **(1)** an SPM12 extension, by adding it to the toolbox folder of SPM and using the batch interface (SPM -> Tools -> CT Segmentation); or **(2)** by interfacing with the code directly (example below).

The orientation matrix in the nifti header of CT scans could be messed up, this means that the atlas will not align with the image data. This is here fixed this by a preprocessing step. Note that this operation requires reslicing of the image data and therefore creates a copy of the original input data (prefixed *r\*.nii*). Setting the ```correct_header``` option of CTseg to ```false``` disables this preprocessing step.

If you find the code useful, please consider citing one of the publications in the *References* section.

## Dependencies

The algorithm requires that the following packages are on the MATLAB path:
* **SPM12:** Download from https://www.fil.ion.ucl.ac.uk/spm/software/spm12/.
* **Shoot toolbox:** Add Shoot folder from the toolbox directory of the SPM source code.
* **Longitudinal toolbox:** Add Longitudinal folder from the toolbox directory of the SPM source code.
* **Multi-Brain toolbox:** Download (or clone) from https://github.com/WTCN-computational-anatomy-group/diffeo-segment.

## Example

Below is a MATLAB snippet that takes as input a CT image (as *.nii*) and produces native space GM, WM, CSF tissue segmentations (*c[1-3]\*.nii*), as well as template space non-modulated (*wc[1-3]\*.nii*) and modulated (*mwc[1-3]\*.nii*) ones. The forward deformation that warps the atlas to the native space CT is also written to disk (*y_\*.nii*).
```
% Set algorithm input
pth_ct = 'CT.nii';  % Path to a CT image
odir = '';  % Output directory (if empty, uses same directory as input data)
% What tissue classes to write to disk (column: native, warped, modulated, 
% row: GM, WM, CSF, SIN, BONE, ST, BG)
tc = [true(3, 3); false(4, 3)];  
def = true;  % Write forward deformation to disk?
correct_header = true;  % Correct orientation matrix?

% Run segmentation+normalisation
CTseg(pth_ct, odir, tc, def, correct_header)
```

## Troubleshooting

* **Segmentation results not as expected:** The model file could not have been found. Make sure that the files *spm_mb_model.mat* and *spm_mb_mu.nii* exist in the directory of CTseg. They are in the *model.zip* file, which should get automatically downloaded and unzipped when the code is executed for the first time.

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
