function test_atlas_download()
% test_atlas_download — Verify atlas download and resolution logic.
%
% Tests the multi-atlas infrastructure without running full segmentation.
% Requires internet access to download atlases from Dropbox.
%
% Usage:
%   >> test_atlas_download

dir_ctseg  = fileparts(mfilename('fullpath'));
dir_models = fullfile(dir_ctseg, 'models');
results    = struct('passed', 0, 'failed', 0);

% Test 1: models/ directory exists
fprintf('Test 1: models/ directory exists... ');
if exist(dir_models, 'dir') == 7
    results = pass(results); else results = fail(results, 'models/ directory not found'); end

% Test 2: prior_CTseg.mat exists in models/
fprintf('Test 2: prior_CTseg.mat in models/... ');
if exist(fullfile(dir_models, 'prior_CTseg.mat'), 'file') == 2
    results = pass(results); else results = fail(results, 'prior_CTseg.mat not found in models/'); end

% Test 3: Mmni.mat exists in models/
fprintf('Test 3: Mmni.mat in models/... ');
if exist(fullfile(dir_models, 'Mmni.mat'), 'file') == 2
    results = pass(results); else results = fail(results, 'Mmni.mat not found in models/'); end

% Test 4: get_atlas_registry returns all 7 entries
fprintf('Test 4: Atlas registry has 7 entries... ');
registry = get_atlas_registry();
names = fieldnames(registry);
if numel(names) == 7
    results = pass(results); else results = fail(results, sprintf('Expected 7 entries, got %d', numel(names))); end

% Test 5: Each registry entry has 'file' and 'url' fields
fprintf('Test 5: Registry entries have file and url fields... ');
all_ok = true;
for i = 1:numel(names)
    entry = registry.(names{i});
    if ~isfield(entry, 'file') || ~isfield(entry, 'url')
        all_ok = false; break;
    end
end
if all_ok
    results = pass(results); else results = fail(results, sprintf('Entry ''%s'' missing file or url', names{i})); end

% Test 6: Download and gunzip each atlas
for i = 1:numel(names)
    name = names{i};
    fprintf('Test 6.%d: Download atlas ''%s''... ', i, name);
    try
        pth = download_atlas(name, dir_models);
        if exist(pth, 'file') == 2
            info = dir(pth);
            if info.bytes > 0
                results = pass(results);
            else
                results = fail(results, 'File is empty');
            end
        else
            results = fail(results, 'File not found after download');
        end
    catch ME
        results = fail(results, ME.message);
    end
end

% Test 7: Re-download skips existing files (no re-download)
fprintf('Test 7: Re-download skips existing... ');
try
    pth = download_atlas('spm15', dir_models);
    if exist(pth, 'file') == 2
        results = pass(results); else results = fail(results, 'File not found'); end
catch ME
    results = fail(results, ME.message);
end

% Test 8: Unknown shorthand produces error
fprintf('Test 8: Unknown shorthand error... ');
try
    download_atlas('nonexistent_atlas', dir_models);
    results = fail(results, 'No error thrown for unknown atlas');
catch ME
    if contains(ME.message, 'Unknown atlas')
        results = pass(results); else results = fail(results, sprintf('Wrong error: %s', ME.message)); end
end

% Test 9: Shorthand detection (isfield check)
fprintf('Test 9: Shorthand vs file path detection... ');
if isfield(registry, 'spm15') && ~isfield(registry, '/path/to/atlas.nii')
    results = pass(results); else results = fail(results, 'Shorthand detection logic wrong'); end

% Summary
fprintf('\n=== Results: %d passed, %d failed ===\n', results.passed, results.failed);
if results.failed > 0
    error('%d test(s) failed.', results.failed);
end

%==========================================================================
function results = pass(results)
fprintf('PASSED\n');
results.passed = results.passed + 1;
%==========================================================================

%==========================================================================
function results = fail(results, msg)
fprintf('FAILED: %s\n', msg);
results.failed = results.failed + 1;
%==========================================================================

%==========================================================================
function registry = get_atlas_registry()
% Copy of the registry from spm_CTseg.m for standalone testing.
registry.default    = struct('file','mu_CTseg.nii',          'url','https://www.dropbox.com/scl/fi/k9v7yccfcknb97860ci8s/mu_CTseg.nii.gz?rlkey=3ps0epza7yv8l18wz7dj1kyjo&st=796t7fyw&dl=1');
registry.spm10      = struct('file','mu_CTseg_spm10.nii',    'url','https://www.dropbox.com/scl/fi/w6t2hzih7ve06p3pki7a3/mu_CTseg_spm10.nii.gz?rlkey=rig0kqp1i280v8dr7e9ggcfqi&st=0zyv9wgf&dl=1');
registry.spm15      = struct('file','mu_CTseg_spm15.nii',    'url','https://www.dropbox.com/scl/fi/4pi8nfmj4uv3zqxcxyvkr/mu_CTseg_spm15.nii.gz?rlkey=tei8ouu7zfamk6ckw8vzyuc0x&st=9t3is2k5&dl=1');
registry.icbm10asym = struct('file','mu_CTseg_icbm10asym.nii','url','https://www.dropbox.com/scl/fi/t2pa9t60rhy8gke4uc9i1/mu_CTseg_icbm10asym.nii.gz?rlkey=4ebid783llfv7zsug44dotrke&st=grj8bu2s&dl=1');
registry.icbm10sym  = struct('file','mu_CTseg_icbm10sym.nii', 'url','https://www.dropbox.com/scl/fi/wlhrppe8hnfd7j6pzdguo/mu_CTseg_icbm10sym.nii.gz?rlkey=fnbns7mao2xopm53snnc7ng2c&st=e3vskk97&dl=1');
registry.icbm15asym = struct('file','mu_CTseg_icbm15asym.nii','url','https://www.dropbox.com/scl/fi/1v9cfdcx15cya25qu6mum/mu_CTseg_icbm15asym.nii.gz?rlkey=bvnte7utnh0u6gj2arg4ajsxl&st=ozerm1vg&dl=1');
registry.icbm15sym  = struct('file','mu_CTseg_icbm15sym.nii', 'url','https://www.dropbox.com/scl/fi/b4w4vo4f3gmouqguit8ux/mu_CTseg_icbm15sym.nii.gz?rlkey=irswa7izexmsumt8mgngv3qa2&st=tlcv24zr&dl=1');
%==========================================================================

%==========================================================================
function pth = download_atlas(name, dir_models)
% Copy of the download function from spm_CTseg.m for standalone testing.
registry = get_atlas_registry();
if ~isfield(registry, name)
    error('Unknown atlas ''%s''. Valid names: %s', ...
        name, strjoin(fieldnames(registry), ', '));
end
entry = registry.(name);
pth   = fullfile(dir_models, entry.file);
if exist(pth, 'file') == 2
    return
end
pth_gz = fullfile(dir_models, [entry.file '.gz']);
fprintf('Downloading atlas ''%s'' (first use only)... ', name)
websave(pth_gz, entry.url);
fprintf('done.\n')
fprintf('Extracting... ')
gunzip(pth_gz, dir_models);
delete(pth_gz);
fprintf('done.\n')
%==========================================================================
