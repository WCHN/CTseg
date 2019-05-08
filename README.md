<img style="float: right;" src="https://github.com/WCHN/CTseg/blob/master/figures/im.png" width="75%" height="75%">
<img style="float: right;" src="https://github.com/WCHN/CTseg/blob/master/figures/seg.png" width="75%" height="75%">

# CTseg

This is a tool for segmenting and spatially normalising routine clinical, computed tomography (CT) scans. It is simillar to the widely used SPM software (see https://www.fil.ion.ucl.ac.uk/spm/software/spm12/), but with a number of key differences, for example:

1. The probabalistic atlas, which is deformed to the subject scan, has been learnt from a combination of MRIs and CTs and contain eight classes, representative of the tissue distribution of routine CT scans (see Fig. TODO).
2. The deformation model combines affine with diffemorphic registration.
3. The gaussian mixture model now has priors on its means and covariances. The prior hyper-parameters has been learnt from a large number of CT scans, capturing their Hounsfield-based intensity distribution, for a more robust segmentation (see Fig. TODO).

## Usage

This tool runs in MATLAB and builds on the SPM12 software. SPM12 will therefore have to be downloaded (from https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and put on the MATLAB path. Furthermore, make sure that both the Shoot and the Longitudinal toolbox of SPM12 (available in `spm/toolbox/{Shoot,Longitudinal}`) is also on the MATLAB path.

Having installed SPM12, add CTseg to SPM12 by copying it into the toolbox folder (`spm/toolbox/CTseg`). CTseg will now be available as a batch job in SPM12:

1. Launch SPM12 from the MATLAB command window: `spm pet`
2. Press `Batch`
3. Press `SPM -> Tools -> CT Segmentation`
4. Double click the `CT scans` option and pick one or more CT images.
5. Press the green play button to run with default options (should work OK).

Depending what output you chose to write to disk, you will have this written to the `Output directory` option. By default this is `./CTseg-Results/image-name`.

## Publication

A publication is currently being written up and will be added here once it has been published, for potential users of the software to cite.
