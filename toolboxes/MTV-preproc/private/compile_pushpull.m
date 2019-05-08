function compile_pushpull(path_to_spm)
% FORMAT compile_pushpull(path_to_spm)

if nargin == 0
    path_to_spm = '';
end
if isempty(path_to_spm)
    path_to_spm = fileparts(which('spm'));
end
if isempty(path_to_spm)
    error('SPM not found');
end
path_to_spm = fullfile(path_to_spm, 'src');
path_to_src = fileparts(which(mfilename));

mex(['-I' path_to_spm], ...
    fullfile(path_to_src, 'pushpull.c'), ...
    fullfile(path_to_spm, 'shoot_boundary.c'));
