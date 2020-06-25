function ctseg = tbx_cfg_CTseg
% Configuration for CTSeg

%--------------------------------------------------------------------------
% data CT Volumes
%--------------------------------------------------------------------------
data          = cfg_files;
data.tag      = 'data';
data.name     = 'CT scans';
data.filter   = 'image';
data.ufilter  = '.*';
data.num      = [0 Inf];
data.help     = {'Select the CT images to segment.'};
data.preview  = @(f) spm_check_registration(char(f));

%--------------------------------------------------------------------------
% odir
%-------------------------------------------------------------------------
odir         = cfg_files;
odir.tag     = 'odir';
odir.name    = 'Output directory';
odir.val{1}  = {''};
odir.help    = {'Segmentations will be written into this directory. If no directory is specified, output is written to same directory as the input images.'};
odir.filter  = 'dir';
odir.ufilter = '.*';
odir.num     = [0 1];

%--------------------------------------------------------------------------
% tc
%--------------------------------------------------------------------------
tc         = cfg_menu;
tc.tag     = 'tc';
tc.name    = 'Tissues';
tc.help    = {'The native space option allows you to produce a tissue class image (c*) that is in alignment with the original. You can also produce spatially normalised versions of the tissue class - both with (mwc*) and without (wc*) modulation.'};
tc.labels = {
    'Native'
    'Unmodulated'
    'Modulated'
    'Native + Unmodulated'
    'Native + Modulated'
    'Unmodulated + Modulated'    
    'Native + Unmodulated + Modulated'
    }';
tc.values = {[1 0 0]
             [0 1 0]
             [0 0 1]
             [1 1 0]
             [1 0 1]
             [0 1 1]
             [1 1 1]}';
tc.val = {[1 1 1]};

%--------------------------------------------------------------------------
% def
%-------------------------------------------------------------------------
def        = cfg_menu;
def.tag    = 'def';
def.name   = 'Deformation Fields';
def.help   = {'Deformation fields can be saved to disk, and used by the Deformations Utility to spatially normalise images to MNI space. Deformations are saved as .nii files, which contain three volumes to encode the x, y and z coordinates.'};
def.labels = {'None'
              'Forward'}';
def.values = {0 1};
def.val    = {1};

%--------------------------------------------------------------------------
% correct_header
%-------------------------------------------------------------------------
correct_header        = cfg_menu;
correct_header.tag    = 'correct_header';
correct_header.name   = 'Correct Orientation Matrix';
correct_header.help   = {'CT images can have messed up orientation matrices in their headers. This means that the atlas will not be able to align with the image data. This resets the orientation matrixes and reslices the image data. OBS: This will create a copy of the input image data and reslice it (prefixed r*).'};
correct_header.labels = {'No'
                         'Yes'}';
correct_header.values = {0 1};
correct_header.val    = {0};

%--------------------------------------------------------------------------
% ctseg
%-------------------------------------------------------------------------
ctseg        = cfg_exbranch;
ctseg.tag    = 'ctseg';
ctseg.name   = 'CT Segmentation';
ctseg.val    = {data odir tc def correct_header};
ctseg.prog   = @ctseg_run;
ctseg.vout   = @vout;
ctseg.help   = {
'This is a CT segmentation algorithm.',...
'',...
'The segmention results are:',...
'    Grey matter',...
'    White matter',...
'    Cerebrospinal fluid'...
'in native and template (normalised) space.',...
'The resulting tissue segmentations are in the',...
'same format as the default SPM12 segmentation routine'
};
%==========================================================================

%==========================================================================
function dep = vout(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated.function dep = vout(job)

cdep = cfg_dep;

for i=1:3
    if job.tc(1)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('c%d Images',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','c','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if job.tc(2)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('wc%d Images',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','wc','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if job.tc(3)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('mwc%d Images',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','mwc','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
end

if job.def
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Forward Deformations';
    cdep(end).src_output = substruct('.','fordef','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

dep = cdep;