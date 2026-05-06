function copy_to_manuscript()
% Copy figures and tables from experiments/results to manuscript/figures.
%
% Run this independently after experiments are complete and results
% have been reviewed. Avoids overwriting manuscript material during
% experiment reruns.
%
% Usage:
%   >> copy_to_manuscript

cfg     = config();
src_dir = cfg.fig_dir;               % experiments/results/
dst_dir = fullfile(cfg.exp_dir, '..', 'manuscript', 'figures');

if ~exist(dst_dir, 'dir')
    mkdir(dst_dir);
end

% Copy figures (.pdf, .png) and tables (.tex)
patterns = {'*.pdf', '*.png', '*.tex'};
n_copied = 0;
for p = 1:numel(patterns)
    files = dir(fullfile(src_dir, patterns{p}));
    for f = 1:numel(files)
        src = fullfile(src_dir, files(f).name);
        dst = fullfile(dst_dir, files(f).name);
        copyfile(src, dst);
        fprintf('  %s\n', files(f).name);
        n_copied = n_copied + 1;
    end
end

fprintf('Copied %d files from %s to %s\n', n_copied, src_dir, dst_dir);
