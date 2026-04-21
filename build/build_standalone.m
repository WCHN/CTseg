function build_standalone(spm_dir, ctseg_dir, out_zip)
% BUILD_STANDALONE  Compile SPM + CTseg into a Linux standalone ZIP.
%
%   build_standalone(SPM_DIR, CTSEG_DIR, OUT_ZIP)
%
% Produces the ZIP consumed by the CTseg Dockerfile. Run from MATLAB on
% Linux (the standalone binary is platform-specific to the host OS).
%
% INPUTS
%   SPM_DIR    Path to an SPM12 source tree (e.g. '~/spm12').
%   CTSEG_DIR  Path to this CTseg repo (e.g. '/mnt/c/.../CTseg').
%   OUT_ZIP    Output ZIP path, named to match the Dockerfile's ADD URL
%              (e.g. '~/spm12_r7771_BI_Linux_R2026a.zip').
%
% PREREQUISITES
%   * MATLAB + MATLAB Compiler installed on Linux
%   * Multi-Brain MEX built for Linux:
%       cd(fullfile(CTSEG_DIR,'mb')); system('make MEXBIN=$(which mex)');
%
% SIDE EFFECTS
%   * Stashes CTSEG_DIR/models/mu_CTseg*.nii to CTSEG_DIR/.build_stash/
%     during compile (atlases are downloaded by the Dockerfile at image
%     build time, not baked into the standalone) and restores them after.
%   * Overwrites SPM_DIR/toolbox/CTseg with a pruned copy of CTSEG_DIR,
%     then removes it when done.
%
% See also: BUILD.md at the repo root.

assert(nargin == 3, 'Usage: build_standalone(spm_dir, ctseg_dir, out_zip)');
assert(isunix && ~ismac, 'Linux only (mcc targets the host OS).');

spm_dir   = char(spm_dir);
ctseg_dir = char(ctseg_dir);
out_zip   = char(out_zip);

assert(exist(fullfile(spm_dir,'spm.m'),'file') == 2, ...
    'SPM not found at %s', spm_dir);
assert(exist(fullfile(ctseg_dir,'spm_CTseg.m'),'file') == 2, ...
    'CTseg not found at %s', ctseg_dir);
assert(exist(fullfile(ctseg_dir,'mb','spm_gmmlib.mexa64'),'file') == 2, ...
    'Multi-Brain MEX not built. Run: cd %s/mb && make', ctseg_dir);

stash_dir  = fullfile(ctseg_dir, '.build_stash');
toolbox_cp = fullfile(spm_dir, 'toolbox', 'CTseg');
build_dir  = tempname;
cleanup    = onCleanup(@() local_cleanup(stash_dir, toolbox_cp, build_dir, ctseg_dir));

% 1. Stash large atlases so they aren't baked into the CTF.
if ~exist(stash_dir,'dir'), mkdir(stash_dir); end
atlases = dir(fullfile(ctseg_dir,'models','mu_CTseg*.nii'));
for i = 1:numel(atlases)
    movefile(fullfile(atlases(i).folder, atlases(i).name), stash_dir);
end

% 2. Replace any existing toolbox/CTseg with a real copy (mcc -a does NOT
%    follow symlinks, so a symlink alone causes silent exclusion).
if exist(toolbox_cp,'dir') || ~isempty(dir(toolbox_cp))
    rmdir_force(toolbox_cp);
elseif exist(toolbox_cp,'file')
    delete(toolbox_cp); % stale symlink
end
system_check(sprintf('cp -rL %s %s', escape(ctseg_dir), escape(toolbox_cp)));

% 2b. Also place mb as a top-level SPM toolbox. spm_make_standalone's config
%     scan is non-recursive: it only sees *_cfg_*.m directly under toolbox/*,
%     so mb's tbx_cfg_mb.m would be missed if mb is only nested under CTseg.
%     Keep the CTseg/mb copy too -- spm_CTseg.m resolves mb via a path
%     relative to its own location.
toolbox_mb = fullfile(spm_dir, 'toolbox', 'mb');
if exist(toolbox_mb,'dir'), rmdir_force(toolbox_mb); end
system_check(sprintf('cp -rL %s %s', ...
    escape(fullfile(ctseg_dir,'mb')), escape(toolbox_mb)));

% 3. Prune non-runtime files to keep the CTF lean.
drop = {'experiments','manuscript','.git','.gitmodules','.gitignore', ...
        '.build_stash','out','test','examples','matlab_workspace.mat', ...
        'Dockerfile','README.md','LICENSE','BUILD.md','scripts','desktop.ini'};
for i = 1:numel(drop)
    p = fullfile(toolbox_cp, drop{i});
    if exist(p,'dir'),  rmdir_force(p); end
    if exist(p,'file'), delete(p); end
end
for pat = {'example_*.png','example_*.svg'}
    d = dir(fullfile(toolbox_cp, pat{1}));
    for i = 1:numel(d), delete(fullfile(d(i).folder, d(i).name)); end
end

% 4. Compile.
mkdir(build_dir);
addpath(spm_dir);
spm_jobman('initcfg');
fprintf('Compiling standalone (this takes ~5 min)...\n');
spm_make_standalone(build_dir);

% 5. Package: ZIP must have a top-level `spm12/` dir so Dockerfile's
%    `unzip -d /opt` lands the binary at `/opt/spm12/spm12`.
stage = fullfile(tempname);
stage_spm12 = fullfile(stage, 'spm12');
mkdir(stage_spm12);
for f = {'spm12','spm12.ctf','run_spm12.sh'}
    copyfile(fullfile(build_dir, f{1}), stage_spm12);
end
if exist(out_zip,'file'), delete(out_zip); end
system_check(sprintf('cd %s && zip -qr %s spm12', escape(stage), escape(out_zip)));
rmdir_force(stage);

info = dir(out_zip);
fprintf('Built: %s (%.1f MB)\n', out_zip, info.bytes / 1024 / 1024);
end

% -------------------------------------------------------------------------
function local_cleanup(stash_dir, toolbox_cp, build_dir, ctseg_dir)
% Restore stashed atlases and remove scratch dirs even on error.
if exist(stash_dir,'dir')
    d = dir(fullfile(stash_dir,'*.nii'));
    for i = 1:numel(d)
        movefile(fullfile(d(i).folder, d(i).name), fullfile(ctseg_dir,'models'));
    end
    try rmdir(stash_dir); catch, end
end
if exist(toolbox_cp,'dir'), rmdir_force(toolbox_cp); end
toolbox_mb = fullfile(fileparts(toolbox_cp), 'mb');
if exist(toolbox_mb,'dir'), rmdir_force(toolbox_mb); end
if exist(build_dir,'dir'),  rmdir_force(build_dir);  end
end

function rmdir_force(p)
[status, msg] = rmdir(p, 's');
if ~status && exist(p,'dir')
    system_check(sprintf('rm -rf %s', escape(p)));
end
end

function system_check(cmd)
[status, out] = system(cmd);
if status ~= 0
    error('Command failed (%d): %s\n%s', status, cmd, out);
end
end

function s = escape(p)
s = ['''' strrep(p, '''', '''\''''') ''''];
end
