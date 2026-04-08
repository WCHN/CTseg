# CTseg Validation Experiments

Reproduces all results for the CTseg validation paper.

## Prerequisites

**MATLAB:**
- MATLAB (R2020a+)
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/download/) on MATLAB path
  - Shoot toolbox (`spm12/toolbox/Shoot`)
  - Longitudinal toolbox (`spm12/toolbox/Longitudinal`)
  - [Multi-Brain toolbox](https://github.com/WTCN-computational-anatomy-group/mb) in `spm12/toolbox/mb`
- [CTseg](https://github.com/WCHN/CTseg) on MATLAB path
- [spm-hospital-preproc](https://github.com/WTCN-computational-anatomy-group/spm-hospital-preproc) — set path in `config.m`
- [PredictPRoNTo](https://github.com/brudfors/PredictPRoNTo) — set path in `config.m`; includes bundled PRoNTo v2
- Image Processing Toolbox (for `bwperim`, `bwdist`)
- Statistics Toolbox (for `signrank`, `prctile`)

**Python (preprocessing only):**
- Python 3 with SimpleITK (`pip install SimpleITK`)

## Data

Paired MR/CT scans from the [SynthRAD2025 Challenge](https://synthrad2025.grand-challenge.org/), Task 1 (Head & Neck). Set the data path in `config.m` (`cfg.data_dir`).

The raw data arrives as `.mha` files. The preprocessing pipeline (below) converts them to NIfTI and produces the `pp_mr.nii.gz` and `pp_ct.nii.gz` files expected by the main pipeline.

## Setup

1. Edit `config.m` to set local paths (`cfg.data_dir`, `cfg.spm_preproc_dir`, `cfg.predict_pronto_dir`).
2. Run preprocessing (see below).
3. Review subjects and update `cfg.exclude` in `config.m` with any that should be excluded (e.g. incomplete brain coverage).

## Preprocessing

Preprocessing must be run before the main pipeline. The scripts are in `preprocessing/`:

| Order | Script | Description | Language |
|-------|--------|-------------|----------|
| 1 | `mha_to_nifti.m` | Convert raw `.mha` files to NIfTI | MATLAB |
| 2 | `preproc.m` | Rigid align to MNI, co-register MR/CT, resample to SPM bounding box | MATLAB |
| 3 | `deform_ct_to_mr.py` | Deformable B-spline registration of CT to MR (replicates SynthRAD2025 preprocessing) | Python |

After preprocessing, each subject folder contains:
- `pp_mr.nii.gz` — preprocessed MRI
- `pp_ct.nii.gz` — preprocessed CT (rigidly aligned)
- `pp_ct_def.nii` — preprocessed CT (deformably registered to MR)

Subjects can be visually inspected with `review/review_imgs.m`.

## Main Pipeline

```matlab
cd experiments
run_all
```

To re-run a specific step:
```matlab
run_all('step4')    % re-compute metrics
run_all('figures')  % regenerate figures
run_all('tables')   % regenerate LaTeX tables
```

To test on a single subject first (recommended):
```matlab
run_all('test')     % all steps, one subject
run_all('test3')    % only step 3, one subject
```

| Step | Script | Description |
|------|--------|-------------|
| 1 | `steps/step1_segment_mr.m` | SPM unified segmentation on MR (silver standard) |
| 2 | `steps/step2_segment_ct_spm.m` | SPM unified segmentation on CT + normalised CT |
| 3 | `steps/step3_segment_ct_ctseg.m` | CTseg on CT + normalised CT |
| 4 | `steps/step4_compute_metrics.m` | Dice, HD95, ASSD between CT and MR segmentations |
| 4b | `steps/step4b_normalisation_metrics.m` | Normalised-space segmentation metrics |
| 5 | `steps/step5_average_normalised.m` | Group-average normalised CTs (sharpness comparison) |
| 6 | `steps/step6_volumetrics.m` | TBV/TIV estimation and agreement |
| 7 | `steps/step7_prediction.m` | Sex classification and age prediction from tissue maps |
| — | `analysis/run_stats.m` | Wilcoxon tests, ICC, summary statistics |
| — | `analysis/make_figures.m` | Publication figures |
| — | `analysis/make_tables.m` | LaTeX tables |

Segmentation results can be inspected with `review/review_segs.m`.

## Atlas Creation

The CTseg atlas must be aligned to MNI space before running the pipeline. Pre-computed atlases are distributed with CTseg; to recreate them:

```matlab
% SPM-aligned (single T1 channel)
create_atlas('mni', 1.5, 'mu_CTseg_spm15')
create_atlas('mni', 1.0, 'mu_CTseg_spm10')

% ICBM 2009c-aligned (T1 + T2 + PD multi-channel)
create_atlas('icbm_asym', 1.5, 'mu_CTseg_icbm15asym')
create_atlas('icbm_asym', 1.0, 'mu_CTseg_icbm10asym')
create_atlas('icbm_sym',  1.5, 'mu_CTseg_icbm15sym')
create_atlas('icbm_sym',  1.0, 'mu_CTseg_icbm10sym')
```

The ICBM variants require the [ICBM 152 NonLinear 2009c](https://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152NLin2009) atlas files (T1, T2, PD) to be placed in `experiments/`.

Grid search scripts for tuning registration parameters are in `atlas/`.

## Output

Intermediate data is saved to `experiments/results/`:
- `metrics.mat` — Dice, HD95, ASSD per subject and tissue
- `norm_metrics.mat` — normalised-space metrics
- `volumetrics.mat` — TBV, TIV per subject
- `prediction.mat` — sex/age prediction results
- `stats.mat` — statistical test results
- `avg_normalised_ct_*.nii` — group-average normalised CTs

Publication figures and tables are saved to `manuscript/figures/`.

## File Structure

```
experiments/
  run_all.m              Top-level orchestrator
  config.m               Paths and settings (edit this)
  README.md              This file
  steps/                 Main pipeline (steps 1-7)
  analysis/              Statistical tests, figures, tables
  utils/                 Shared helper functions
  atlas/                 Atlas creation and parameter grid searches
  preprocessing/         Data conversion and preprocessing
  review/                Interactive image and segmentation viewers
```
