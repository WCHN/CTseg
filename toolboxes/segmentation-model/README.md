# A generative model for learning from neuroimaging data

This code implements a generative model for learning from neuroimaging data. Once learnt, the model can be used to segment images of the brain.

## Example

To familiarise a user with the framework, the demo folder contains the script:

**train_MRBrainS18_2d**: Demonstrates using the segmentation-model code to learn from
MR scans part of the MRBrainS18 segmentation challenge (http://mrbrains18.isi.uu.nl/).

## Dependencies

This project has strong dependencies to SPM12 and its `Shoot` toolbox. Both of them should be added to Matlab's path. SPM can be downloaded at [www.fil.ion.ucl.ac.uk/spm](http://www.fil.ion.ucl.ac.uk/spm/).

Core functions also depend on our [`auxiliary-functions` toolbox](https://github.com/WTCN-computational-anatomy-group/auxiliary-functions), which gathers lots of low-level functions.

Furthermore, executable scripts depend on our [`distributed-computing` toolbox](https://github.com/WTCN-computational-anatomy-group/distributed-computing), which helps parallelising parts of the code either on the local workstation (using Matlab's parallel processing toolbox) or on a computing cluster (see the toolbox help file for use cases and limitations).
