% clean_pp.m — Delete all pp_ files from the HN data folder.
addpath(fileparts(fileparts(mfilename('fullpath'))));  % experiments/
cfg = config();
files = dir(fullfile(cfg.data_dir, '**', 'pp_*'));
fprintf('Found %d pp_ files\n', numel(files));
for i = 1:numel(files)
    delete(fullfile(files(i).folder, files(i).name));
end
fprintf('Deleted %d files\n', numel(files));
