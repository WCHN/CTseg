# Building the CTseg Docker image

The Docker image in `Dockerfile` downloads a precompiled SPM + CTseg standalone binary from GitHub Releases. This document is for the maintainer who needs to **rebuild that binary** after changes to CTseg code or models.

End users who just want to run CTseg in Docker should follow the "Docker" section of the main `README.md`; they do not need anything here.

## Overview

The standalone is SPM12 with CTseg bundled as a toolbox, compiled via MATLAB Compiler (`mcc`) and run against MATLAB Runtime (MCR) inside the container. Producing it requires:

1. A Linux machine with MATLAB + MATLAB Compiler installed (the binary is platform-specific — Windows MATLAB cannot build the Linux binary).
2. A local clone of SPM12 source.
3. This CTseg repo, with the Multi-Brain MEX (`spm_gmmlib.mexa64`) built for Linux.

Output: a ZIP (e.g. `spm12_r7771_BI_Linux_R2026a.zip`) uploaded to a GitHub Release, plus a Dockerfile update pointing to it.

## 1. Set up MATLAB on Linux

These instructions assume WSL2 Ubuntu on Windows 11. On a native Linux box, skip the WSL-specific bits.

Install OS prerequisites (Ubuntu 24.04 ships `t64` variants; adjust for older Ubuntu):

```bash
sudo apt update
sudo apt install -y unzip build-essential \
  libxt6t64 libxmu6 libxcomposite1 libasound2t64 \
  libgtk-3-0t64 libnss3 libgbm1 libxkbcommon0 libxrandr2 \
  libxdamage1 libxfixes3 libegl1 libglu1-mesa \
  libxshmfence1 xdg-utils
```

From Windows, download the MATLAB Linux installer from mathworks.com. Copy into WSL, extract, run:

```bash
mkdir -p ~/matlab_install && cd ~/matlab_install
cp /mnt/c/Users/<you>/Downloads/matlab_*_glnxa64.zip .
unzip -q matlab_*_glnxa64.zip -d installer
cd installer && sudo ./install
```

Select **MATLAB + MATLAB Compiler** at minimum. Image Processing and Statistics toolboxes are useful for other CTseg workflows but aren't required for the standalone build.

Verify:

```bash
/usr/local/MATLAB/<release>/bin/matlab -batch "ver; disp(license('test','Compiler')); disp(which('mcc'))"
```

## 2. Prepare the build tree

Clone SPM12 (the public GitHub mirror ends at r7771 — that is the last public SPM12 release; r8168 in the original Docker came from an internal FIL build no longer available):

```bash
cd ~
git clone https://github.com/spm/spm12.git
```

Build the Multi-Brain MEX for Linux. The Makefile hardcodes an old MATLAB path — override it:

```bash
cd /path/to/CTseg/mb
make MEXBIN=/usr/local/MATLAB/<release>/bin/mex
# produces spm_gmmlib.mexa64
```

## 3. Compile the standalone

From MATLAB on Linux:

```matlab
addpath('/path/to/CTseg/build');
build_standalone('/home/<you>/spm12', ...
                 '/path/to/CTseg', ...
                 '/home/<you>/spm12_r7771_BI_Linux_R2026a.zip');
```

The script (see `build/build_standalone.m`):

- Stashes `models/mu_CTseg*.nii` so the large atlases are *not* baked into the CTF — the Dockerfile downloads them at image build time instead.
- Copies (not symlinks — `mcc -a` silently skips symlinks) a pruned CTseg tree into `<spm12>/toolbox/CTseg`.
- Runs `spm_jobman('initcfg')` + `spm_make_standalone`.
- Packages the three output files (`spm12`, `spm12.ctf`, `run_spm12.sh`) under a `spm12/` prefix in the ZIP so the Dockerfile's `unzip -d /opt` lands them at `/opt/spm12/`.
- Restores the stashed atlases and cleans up the scratch copy in SPM's toolbox dir on success or error.

Expect ~5 minutes and a ~250 MB ZIP.

## 4. Host the ZIP

Two options:

**GitHub Releases (recommended long-term)** — stable versioned URL, no expiring tokens. Cut a release on `github.com/WCHN/CTseg` (e.g. tag `standalone-<spm_rev>-<matlab_release>`), upload the ZIP as an asset. URL pattern:

```
https://github.com/WCHN/CTseg/releases/download/<tag>/spm12_<spm_rev>_BI_Linux_<matlab_release>.zip
```

Asset size limit is 2 GB — our ZIP is well under.

**Dropbox (interim)** — fine for feature branches before a release is cut. Upload via the Dropbox web UI, grab the share link, and append `&dl=1` so `ADD` / `wget` gets the raw file rather than the preview page. Links include a per-file `rlkey` and `st` token; these have broken before (the original FIL-hosted URL rotted in July 2024, prompting the first Dropbox switch), so treat this as temporary.

## 5. Update the Dockerfile

Bump the version variables and the download URL. Four lines change:

```dockerfile
ENV MATLAB_VERSION=R2026a
ENV MCR_DIR=R2026a

# Filename is now `MATLAB_Runtime_R2026a_glnxa64.zip` (no `_Update_N_` suffix).
RUN wget ... https://ssd.mathworks.com/supportfiles/downloads/${MATLAB_VERSION}/Release/0/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_${MATLAB_VERSION}_glnxa64.zip ...

ENV SPM_REVISION=r7771

ADD <url-to-zip> /opt
```

Where `<url-to-zip>` is either the GitHub Releases URL or a Dropbox share link with `&dl=1` appended.

The MCR subdirectory name inside `/opt/mcr/` matches the MATLAB release (e.g. `/opt/mcr/R2026a/`) in recent MathWorks installers. Older releases (through R2023b) used a numeric versioned name like `v911` instead. `LD_LIBRARY_PATH` uses `${MCR_DIR}` — check the installer output or `docker run --rm --entrypoint=ls ubuntu:ctseg /opt/mcr/` and adjust if needed.

## 6. Build and test the Docker image

```bash
docker build -t ubuntu:ctseg -f Dockerfile .
```

The build runs `spm_CTseg(1)` near the end, which triggers a download of the default atlas (`mu_CTseg_spm15.nii`, ~66 MB) into the image.

End-to-end test (replace the volume mount and filename with your own CT):

```bash
docker run --rm -it \
  -v /mnt/c/Users/mbrudfors/Data/CT-brain:/data \
  ubuntu:ctseg \
  eval "spm_CTseg('/data/CT.nii', '', true, true, false, false, NaN, [], [], 'spm15', false, 'spm')"
```

Confirm that tissue maps (`c1..c6`, `wc*`, `mwc*`) appear in the mounted directory and match a MATLAB-native reference run on the same input.

## Gotchas

- **`mcc -a` does not follow symlinks.** A symlinked `spm12/toolbox/CTseg` compiles without errors but silently omits the toolbox. `build_standalone.m` uses `cp -rL` to sidestep this.
- **Build must run on Linux.** `mcc` produces binaries for the host OS only. A Windows MATLAB can't produce the Linux binary.
- **MATLAB Update version vs. MCR.** MathWorks may lag in publishing MCR updates. We compile with MATLAB R2026a Update 1 but the public MCR for R2026a is `Release/0`; the runtime is forward-compatible with minor updates within a release.
- **Large atlases.** `models/mu_CTseg*.nii` (up to 224 MB each) must be excluded from the CTF or the standalone ZIP balloons. `build_standalone.m` stashes and restores them automatically.
- **`spm_make_standalone` toolbox scan is non-recursive.** It looks for `*_cfg_*.m` only in immediate children of `spm12/toolbox/`, so Multi-Brain's `tbx_cfg_mb.m` (normally nested at `CTseg/mb/`) is missed unless `mb` is also placed at `spm12/toolbox/mb/`. `build_standalone.m` bundles `mb` at both locations.
- **`addpath`, `rmpath`, and `mex` are forbidden in deployed code.** `spm_CTseg.m` must guard all such calls with `if ~isdeployed` — in the standalone, paths are baked in at compile time and MEX auto-compilation is unavailable. The batch config registration that the runtime symlink dance is meant to achieve is already done by `build_standalone.m` at compile time.
