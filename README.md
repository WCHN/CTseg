# CTseg
...

## Usage

This tool runs in MATLAB and builds on the SPM12 software. SPM12 will therefore have to be downloaded ([www.fil.ion.ucl.ac.uk/spm](from http://www.fil.ion.ucl.ac.uk/spm/)) and put on the MATLAB path. Furthermore, make sure that both the Shoot and the Longitudinal toolbox of SPM12 (available in `spm/toolbox/{Shoot,Longitudinal}`) is also on the MATLAB path.

Having installed SPM12, add CTseg to SPM12 by copying it into the toolbox folder (`spm/toolbox/CTseg`). CTseg will now be available as a batch job in SPM12:

1. Launch SPM12 from the MATLAB command window: `spm pet`
2. Press `Batch`
3. Press `SPM -> Tools -> CT Segmentation`
4. Double click the `CT scans` option and pick one or more CT images.
5. Press the green play button to run with default options (should work OK).

Depending what output you chose to write to disk, you will have this written to the `Output directory` option. By default this is `.CTseg-Results/`.
