# CTseg: A Tool for Brain CT Segmentation, Spatial Normalisation, and Volumetrics

 [![arXiv](https://img.shields.io/badge/arXiv-2605.05154-b31b1b.svg)](https://arxiv.org/abs/2605.05154)

> **New:** CTseg paper has been published on arXiv — *CTseg: A Tool for Brain CT Segmentation, Spatial Normalisation, and Volumetrics* — [arXiv:2605.05154](https://arxiv.org/abs/2605.05154) (2026).

<img style="float: right;" src="https://github.com/WCHN/CTseg/blob/master/example_1.png" width="80%" height="80%">
<img style="float: right;" src="https://github.com/WCHN/CTseg/blob/master/example_2.png" width="80%" height="80%">

## Table of Contents

- [Overview](#overview)
- [Installation Instructions](#installation-instructions)
- [Docker (no MATLAB needed)](#docker-no-matlab-needed)
- [Available Atlases](#available-atlases)
- [Hemisphere Segmentation](#hemisphere-segmentation)
- [Example use cases](#example-use-cases)
- [Troubleshooting](#troubleshooting)
- [Improved runtime (Linux and Mac)](#improved-runtime-linux-and-mac)
- [References](#references)
- [License](#license)

## Overview

This is an algorithm for segmenting and spatially normalising computed tomography (CT) brain scans. The model is an extension of the popular unified segmentation routine (part of the SPM12 software) with: improved registration, priors on the Gaussian mixture model parameters, an atlas learned from both MRIs and CTs (with more classes). These improvements leads to a more **robust** segmentation routine that can better handle image with lots of noise and/or large morphological variability (see figure above). The algorithm can produce native|warped|modulated space segmentations of:

1. Gray matter (GM)
2. White matter (WM)
3. Cerebrospinal fluid (CSF)
4. Skull (BONE)
5. Soft tissue (ST)
6. Background (BG)

The implementation is done in MATLAB and depends on the SPM12 package (and its MB toolbox), but can be run without MATLAB using Docker. The dependencies are packaged in the latest release. If you find the code useful, please consider citing the publications in the *References* section.

#### Further details

The input to CTseg should be provided as NIfTI files (```.nii```). The resulting tissue segmentations are in the same format as the output of the SPM12 segmentation routine (```c*```, ```wc*```, ```mwc*```). The normalised segmentations (```wc*```, ```mwc*```) are in atlas space.

A **skull-stripped** version of the input image is produced by default (prefixed ```ss_``` to the original filename). **Total brain volume** (TBV) and **intercranial volume** (TIV) are also computed by the algorithm and returned as the second argument of the CTseg function. Note that both of these routines uses only the GM, WM and CSF segmentations of the algorithm. The skull-stripped volume will therefore not include the meninges, the sinuses or any calcifications; the TIV might therefore also be slighly underestimated.

For converting **DICOM** CT to NIfTI, we recommend using SPM12's ```spm_dicom_convert```. This DICOM converter can deal with the fact that many CT images are often acquired with variable slice thickness. If this is not accounted for when reconstructing the NIfTI file, the head shape can be deformed.

## Installation Instructions

The algorithm is developed using MATLAB and relies on external functionality from the SPM12 software. The following are required:

* **SPM12:** Download from https://www.fil.ion.ucl.ac.uk/spm/software/download/ and add to the MATLAB path.
* **Shoot toolbox:** The Shoot folder from the toolbox directory of SPM12 (add to path).
* **Longitudinal toolbox:** The Longitudinal folder from the toolbox directory of SPM12 (add to path).
* **Multi-Brain toolbox:** Included as a git submodule. After cloning CTseg, initialise and compile it:

``` bash
git clone --recursive https://github.com/WCHN/CTseg
cd CTseg/mb
make
```

If you have already cloned without `--recursive`, run:

``` bash
git submodule update --init
cd mb
make
```

On Windows, if `make` is not available, compile from MATLAB:

``` matlab
cd mb
mex -O -largeArrayDims spm_gmmlib.c gmmlib.c
```

Note that CTseg will attempt to compile automatically on first use if the MEX file is not found.

## Docker (no MATLAB needed)

CTseg can be run with Docker, which does *not* require you to have MATLAB installed on your computer. First build the image:

```bash
docker build -t ubuntu:ctseg .
```

CTseg can then be run, for example, by:

```bash
docker run --rm -it -v dir_host:/data ubuntu:ctseg function spm_CTseg '/data/CT.nii'
```
 or 
 
```bash
docker run --rm -it -v dir_host:/data ubuntu:ctseg eval "spm_CTseg('/data/CT.nii')"
```

where `dir_host` is the absolute path to a folder on your local machine that contains a `CT.nii` image. After CTseg has finished running, its output can be found in the `dir_host` folder.

**Atlas availability inside Docker**: only the default atlas (`spm15`, 66 MB) is bundled in the image — it is downloaded to `/opt/spm12/spm12_mcr/.../toolbox/CTseg/models/` at image build time. If you pass `'spm10'` or `'ctseg'` as the `mu` argument, `spm_CTseg` will try to download that atlas on first use, which only works if the container has outbound network access. The simplest workaround is to pre-download the atlas on the host and bind-mount it into the container's models directory.

Maintainers: see [`build/README.md`](build/README.md) for how to rebuild the standalone binary after changes to CTseg code or models.

## Available Atlases

The CTseg atlas, when learned, is in the group-wise optimal space of the training data; however, for downstream analysis it could be important to have the deformations map to the space of the SPM Tissue Probability Map (TPM). Therefore, CTseg provides three atlases: the `ctseg` atlas, mapping to original group-wise optimal space; and two atlases mapping to the space of the SPM atlas `spm15`/`spm10`, produced by directly registering the CTseg template to SPM's `TPM.nii` with the Multi-Brain toolbox (categorical, tissue-probability-to-tissue-probability alignment). Warped segmentations (`wc*`, `mwc*`) for `spm15`/`spm10` therefore land in SPM space directly. Atlases are downloaded automatically on first use to the `models/` directory. The `bb` parameter allows to set the field-of-view (FOV) of the template-space segmentations, where `full` gives the original FOV of the input atlas (larger, includes spinal cord, etc) and `spm` gives the FOV of the SPM atlas (default). A smaller resolution should give more accurate segmentations, but with increased runtime and memory use.

| Shorthand   | Space                | Resolution | Size   |
|-------------|----------------------|------------|--------|
| `'spm15'` (default)  | SPM TPM  | 1.5mm      | 66 MB  |
| `'spm10'`   | SPM TPM  | 1.0mm      | 224 MB |
| `'ctseg'`   | CTseg     | 1.0mm      | 224 MB |

Pass the shorthand as the `mu` parameter (10th argument) of `spm_CTseg`:

``` matlab
% Default: SPM-aligned 1.5 mm atlas (downloads mu_CTseg_spm15.nii on first use)
res = spm_CTseg('CT.nii');

% SPM-aligned 1.0 mm atlas
res = spm_CTseg('CT.nii', '', true, true, true, false, NaN, [], [], 'spm10');

% Legacy groupwise optimal atlas
res = spm_CTseg('CT.nii', '', true, true, true, false, NaN, [], [], 'ctseg');

% Use the full input-atlas FOV (disable default SPM TPM cropping)
res = spm_CTseg('CT.nii', '', true, true, true, false, NaN, [], [], 'spm15', false, 'full');

% Use a custom atlas file path
res = spm_CTseg('CT.nii', '', true, true, true, false, NaN, [], [], '/path/to/my_atlas.nii');
```

## Hemisphere Segmentation

CTseg can optionally separate GM and WM into left and right hemisphere segmentations (8 tissue classes instead of 6). Enable this by setting the `hemisphere` parameter (11th argument) to `true`:

``` matlab
% Standard segmentation (6 classes: GM, WM, CSF, Bone, ST, BG)
res = spm_CTseg('CT.nii');

% Hemisphere segmentation (8 classes: GM-L, GM-R, WM-L, WM-R, CSF, Bone, ST, BG)
res = spm_CTseg('CT.nii', '', true, true, true, false, NaN, [], [], '', true);
```

This can be combined with any atlas.

## Example use cases

Below are two MATLAB snippets. The first takes as input a CT image (as ```*.nii```) and produces native space GM, WM, CSF tissue segmentations (```c[1-3]*.nii```), as well as template space, non-modulated (```wc[1-3]*.nii```) and modulated (```mwc[1-3]*.nii```) segmentations. The forward deformation that warps the atlas to the native space CT is also written to disk (as ```y_*.nii```). The second snippet uses the forward deformation to warp: (1) the CT image to the template space; and (2), the template to the space of the CT image. Results are written to ```dir_out```.

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
correct_header = false;  

% Do skull-stripping?
ss = true;

% Template space voxel size
vox = 1.0;

% Spatial regularisation (empty = default, scalar = multiplier)
v_settings = [0.0001 0 0.4 0.1 0.4] * 3;

% Stopping tolerance
tol = 0.001;

% Atlas: shorthand name or path
mu = 'spm15';

% Hemisphere segmentation (split GM/WM into left/right)?
hemisphere = false;

% Bounding box for template-space outputs
% 'spm' (default) = SPM TPM FOV; 'full' = full atlas FOV; or a 2x3 numeric matrix
bb = 'spm';

% Run segmentation+normalisation
[res,vol] = spm_CTseg(pth_ct, dir_out, tc, def, correct_header, ss, vox, ...
                       v_settings, tol, mu, hemisphere, bb)
% res: a struct with paths to result niftis
% vol: a struct containing TBV and TIV
```

### 2. Warping with the generated deformations

``` matlab
% Path to tissue template (located in the CTseg folder)
pth_mu = 'mu_CTseg.nii';

% Path to forward deformation (in dir_out)
pth_y = 'y_*.nii';

% image-to-template (pull)
% Use forward deformation (pth_y) to warp CT image (pth_ct) to 
% template space (pth_mu) by pulling
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def     = {pth_y};
matlabbatch{1}.spm.util.defs.comp{1}.inv.space           = {pth_mu};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames          = {pth_ct};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {dir_out};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp          = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask            = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm            = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix          = 'wpull';
spm_jobman('run',matlabbatch);

% image-to-template (push)
% Use forward deformation (pth_y) to warp CT image (pth_ct) to 
% template space (pth_mu) by pushing (with modulation and smoothing)
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.def                 = {pth_y};
matlabbatch{1}.spm.util.defs.out{1}.push.fnames          = {pth_ct};
matlabbatch{1}.spm.util.defs.out{1}.push.weight          = {''};
matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {dir_out};
matlabbatch{1}.spm.util.defs.out{1}.push.fov.file        = {pth_mu};
matlabbatch{1}.spm.util.defs.out{1}.push.preserve        = 1;
matlabbatch{1}.spm.util.defs.out{1}.push.fwhm            = [10 10 10];
matlabbatch{1}.spm.util.defs.out{1}.push.prefix          = 'wpush';
spm_jobman('run',matlabbatch);

% template-to-image
% Use forward deformation (pth_y) to warp template (pth_mu) to native 
% CT image space (pth_ct) by pulling
matlabbatch = {};
matlabbatch{1}.spm.util.defs.comp{1}.def                 = {pth_y};
matlabbatch{1}.spm.util.defs.comp{2}.id.space            = {pth_ct};
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

* **Segmentation results not as expected:** The model file could not have been found. Make sure that ```models/prior_CTseg.mat``` exists in the CTseg directory (included in the repository). Atlas files (```mu_CTseg*.nii```) are downloaded automatically to the ```models/``` directory on first use.

* **Error related to spm_diffeo:** This code uses a recent version of SPM12; therefore, if your SPM12 version is quite old, the function ```spm_diffeo``` might break. Updating to the latest version of SPM12 will resolve this issue.

## Improved runtime (Linux and Mac)

For a faster algorithm, consider compiling SPM with OpenMP support. Just go to the *src* folder of SPM and do:

``` bash
make distclean
make USE_OPENMP=1 && make install
```

## References

``` latex
@misc{brudfors2026ctsegtoolbrainct,
      title={CTseg: A Tool for Brain CT Segmentation, Spatial Normalisation, and Volumetrics}, 
      author={Mikael Brudfors},
      year={2026},
      eprint={2605.05154},
      archivePrefix={arXiv},
      primaryClass={eess.IV},
      url={https://arxiv.org/abs/2605.05154}, 
}
```

## License

CTseg is free but copyright software, distributed under the terms of the GNU General Public Licence as published by the Free Software Foundation (either version 2, or at your option, any later version).
