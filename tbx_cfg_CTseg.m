function ctseg = tbx_cfg_CTseg
% Simple configuration for Mikael's CT segmentation

%--------------------------------------------------------------------------
% data CT Volumes
%--------------------------------------------------------------------------
data          = cfg_files;
data.tag      = 'data';
data.name     = 'CT scans';
data.filter   = 'image';
data.ufilter  = '.*';
data.num      = [0 Inf];
data.help     = {'Select the CT images for segmentation.'};
data.preview  = @(f) spm_check_registration(char(f));

%--------------------------------------------------------------------------
% outdir Output directory
%-------------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.val{1}  = {'./CTseg-Results'};
outdir.help    = {'Segmentation results will be written into this output directory.',...
                  ' If no directory is specified, output is written to a folder CTseg-Results in the working directory.'};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];

%--------------------------------------------------------------------------
% cleanbrain Clean-up brain
%-------------------------------------------------------------------------
cleanbrain        = cfg_menu;
cleanbrain.tag    = 'cleanbrain';
cleanbrain.name   = 'Clean up any partitions';
cleanbrain.help   = {'This uses a crude routine for extracting the brain from segmented images.'}';
cleanbrain.labels = {'Dont do clean-up'
                     'Do clean-up'}';
cleanbrain.values = {0 1};
cleanbrain.val    = {0};

%--------------------------------------------------------------------------
% mrf MRF Parameter
%--------------------------------------------------------------------------
mrf         = cfg_entry;
mrf.tag     = 'mrf';
mrf.name    = 'MRF Parameter';
mrf.help    = {'When tissue class images are written out, a few iterations of a simple Markov Random Field (MRF) cleanup procedure are run.  This parameter controls the strength of the MRF. Setting the value to zero will disable the cleanup.'};
mrf.strtype = 'r';
mrf.num     = [1 1];
mrf.val     = {1};

%--------------------------------------------------------------------------
% samp Sampling distance
%--------------------------------------------------------------------------
samp         = cfg_entry;
samp.tag     = 'samp';
samp.name    = 'Sampling distance';
samp.help    = {
    'This encodes the approximate distance between sampled points when estimating the model parameters.'
    'Smaller values use more of the data, but the procedure is slower and needs more memory. Determining the ``best'''' setting involves a compromise between speed and accuracy.'
    }';
samp.strtype = 'r';
samp.num     = [1  1];
samp.val     = {2};

%--------------------------------------------------------------------------
% writeim Save Images
%--------------------------------------------------------------------------
writeim         = cfg_menu;
writeim.tag     = 'image';
writeim.name    = 'Images';
writeim.help    = {
    'Save native and/or template space image.'
    }';
writeim.labels = {
                'None'
                'Native'
                'Template'
                'Native + Template'
                }';
writeim.values = {
                [0 0]
                [1 0]
                [0 1]
                [1 1]
                }';
writeim.val    = {[1 0]};

%--------------------------------------------------------------------------
% native Native Tissue
%--------------------------------------------------------------------------
native         = cfg_menu;
native.tag     = 'native';
native.name    = 'Native Tissue';
native.help    = {'The native space option allows you to produce a tissue class image (c*) that is in alignment with the original/* (see Figure \ref{seg1})*/. It can also be used for ``importing'''' into a form that can be used with the Dartel toolbox (rc*).'};
native.labels = {
    'None'
    'Native Space'
    'Dartel Imported'
    'Native + Dartel Imported'
    }';
native.values = {
                 [0 0]
                 [1 0]
                 [0 1]
                 [1 1]
                 }';
native.val    = {[1 0]};

%--------------------------------------------------------------------------
% warped Warped Tissue
%--------------------------------------------------------------------------
warped         = cfg_menu;
warped.tag     = 'warped';
warped.name    = 'Warped Tissue';
warped.help    = {
    'You can produce spatially normalised versions of the tissue class - both with (mwc*) and without (wc*) modulation (see below). These can be used for voxel-based morphometry. All you need to do is smooth them and do the stats.'
    ''
    '``Modulation'''' is to compensate for the effect of spatial normalisation.  When warping a series of images to match a template, it is inevitable that volumetric differences will be introduced into the warped images.  For example, if one subject''s temporal lobe has half the volume of that of the template, then its volume will be doubled during spatial normalisation. This will also result in a doubling of the voxels labelled grey matter.  In order to remove this confound, the spatially normalised grey matter (or other tissue class) is adjusted by multiplying by its relative volume before and after warping.  If warping results in a region doubling its volume, then the correction will halve the intensity of the tissue label. This whole procedure has the effect of preserving the total amount of grey matter signal in the normalised partitions.  Actually, in this version of SPM the warped data are not scaled by the Jacobian determinants when generating the "modulated" data.  Instead, the original voxels are projected into their new location in the warped images.  This exactly preserves the tissue count, but has the effect of introducing aliasing artifacts - especially if the original data are at a lower resolution than the warped images.  Smoothing should reduce this artifact though.'
    'Note also that the "unmodulated" data are generated slightly differently in this version of SPM. In this version, the projected data are corrected using a kind of smoothing procedure. This is not done exactly as it should be done (to save computational time), but it does a reasonable job. It also has the effect of extrapolating the warped tissue class images beyond the range of the original data.  This extrapolation is not perfect, as it is only an estimate, but it may still be a good thing to do.'
    }';
warped.labels = {
                 'None'
                 'Modulated'
                 'Unmodulated'
                 'Modulated + Unmodulated'
                 }';
warped.values = {
                 [0 0]
                 [1 0]
                 [0 1]
                 [1 1]
                 }';
warped.val    = {[0 0]};

%--------------------------------------------------------------------------
% ctseg Segment
%-------------------------------------------------------------------------
ctseg        = cfg_exbranch;
ctseg.tag    = 'ctseg';
ctseg.name   = 'CT Segmentation';
ctseg.val    = {data outdir writeim native warped cleanbrain samp mrf};
ctseg.prog   = @ctseg_run;
ctseg.vout   = @vout;
ctseg.help   = {
'This is Mikael''s CT segmentation algorithm. /* See his PhD thesis \cite{brudfors2019} for more information.*/',...
'',...
'The segmention results are:',...
'    Grey matter.',...
'    White matter.',...
'    More stuff.'
};


function dep = vout(job)
for i=1:8
    dep(i) = cfg_dep;
    dep(i).sname      = sprintf('Tissue %d',i);
    dep(i).src_output = substruct('.','tiss','()',{i},'.','c','()',{':'});
    dep(i).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

