# CTseg

<img style="float: right;" src="https://github.com/WCHN/CTseg/blob/master/demo.png" width="60%" height="60%">

This is a MATLAB implementation of a model for segmenting and spatially normalising computed tomography (CT) brain scans. The model is an extension of the popular unified segmentation routine (part of the SPM12 software) with: improved registration, priors on the Gaussian mixture model parameters, an atlas learned from both MRIs and CTs (with more classes), and more. 

The segmentation results are **grey matter (GM)**, **white matter (WM)** and **cerebrospinal fluid (CSF)**, in native and template (normalised) space. The input should be provided as nifti files (*.nii*), the resulting tissue segmentations are in the same format as the output of the SPM12 segmentation routine. 

The code can be used either as: (1) an SPM12 extension, by adding it to the toolbox folder of SPM and using the batch interface (SPM -> Tools -> CT Segmentation); or (2) by interfacing with the code directly (example below).

If you find the code useful, please consider citing one of the publications in the *References* section.

## Dependencies

The algorithm requires that the following packages are on the MATLAB path:
* **SPM12:** Download from https://www.fil.ion.ucl.ac.uk/spm/software/spm12/.
* **diffeo-segment:** Download (or clone) from https://github.com/WTCN-computational-anatomy-group/diffeo-segment.
* **auxiliary-functions:** Download (or clone) from https://github.com/WTCN-computational-anatomy-group/auxiliary-functions.

## Example

Below is a MATLAB snippet that takes as input a CT image (as *.nii*) and produces native space GM, WM, CSF tissue segmentations (*c[1-3]\*.nii*), as well as template space non-modulated (*wc[1-3]\*.nii*) and modulated (*mwc[1-3]\*.nii*) ones. The forward deformation that warps the atlas to the native space CT is also written to disk (*y_\*.nii*).
```
# Set algorithm input
pth_ct = 'CT.nii';  % Path to a CT image
odir = '';  % Output directory (if empty, use same as input CT)
tc = [1, 1, 1];  % Tissue classes to write to disk [native, unmodulated, modulated]
def = true;  % Write forward deformation to disk?
correct_header = false;  % Correct orientation matrix? (CT images can have messed up header information)

# Run segmentation+normalisation
CTseg(pth_ct, odir, tc, def, correct_header)
```

## Improved runtime (Linux only)

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
