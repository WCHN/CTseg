# CTseg Validation Experiments

Reproduces all results for the CTseg validation paper.

## Prerequisites

- MATLAB (R2020a+)
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/download/) on MATLAB path
  - Shoot toolbox (`spm12/toolbox/Shoot`)
  - Longitudinal toolbox (`spm12/toolbox/Longitudinal`)
  - [Multi-Brain toolbox](https://github.com/WTCN-computational-anatomy-group/mb) in `spm12/toolbox/mb`
- [CTseg](https://github.com/WCHN/CTseg) on MATLAB path
- [spm-hospital-preproc](https://github.com/WTCN-computational-anatomy-group/spm-hospital-preproc) — set path in `config.m` (`cfg.spm_preproc_dir`)
- [PredictPRoNTo](https://github.com/brudfors/PredictPRoNTo) — set path in `config.m` (`cfg.predict_pronto_dir`); includes bundled PRoNTo v2
- Image Processing Toolbox (for `bwperim`, `bwdist`)
- Statistics Toolbox (for `signrank`, `prctile`)
- Python 3 with SimpleElastix (for `deform_ct_to_mr.py` preprocessing step only)

## Data

Paired MR/CT scans from the [SynthRAD2025 Challenge](https://synthrad2025.grand-challenge.org/), Task 1 (Head & Neck). The data path is set in `config.m`.

Each subject folder should contain preprocessed files:
- `pp_mr.nii.gz` — preprocessed MRI
- `pp_ct.nii.gz` — preprocessed CT

## Quick Start

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

To test the full pipeline on a single subject (recommended before a full run):
```matlab
run_all('test')     % all steps, one subject
run_all('test3')    % only step 3, one subject
```

To choose which subject to test on, edit `config.m`:
```matlab
cfg.test_subject = '1HNA001';  % or leave empty for the first available
```

## Workflow

| Step | Script | Description | Speed |
|------|--------|-------------|-------|
| 1 | `steps/step1_segment_mr.m` | SPM unified segmentation on MR (silver standard) | Slow |
| 2 | `steps/step2_segment_ct_spm.m` | SPM unified segmentation on CT + normalised CT | Slow |
| 3 | `steps/step3_segment_ct_ctseg.m` | CTseg on CT + normalised CT | Slow |
| 4 | `steps/step4_compute_metrics.m` | Dice, HD95, ASSD | Moderate |
| 5 | `steps/step5_average_normalised.m` | Average normalised CTs (sharpness comparison) | Moderate |
| 6 | `steps/step6_volumetrics.m` | TBV/TIV from MR-SPM and CTseg | Moderate |
| — | `analysis/run_stats.m` | Wilcoxon tests, ICC, summary statistics | Fast |
| — | `analysis/make_figures.m` | Boxplots, Bland-Altman plots | Fast |
| — | `analysis/make_tables.m` | LaTeX tables | Fast |

## Output

Intermediate data is saved to `experiments/results/`:
- `metrics.mat` — Dice, HD95, ASSD per subject/tissue
- `volumetrics.mat` — TBV, TIV per subject
- `stats.mat` — statistical test results
- `avg_normalised_ct_spm.nii` — group average normalised CT (SPM)
- `avg_normalised_ct_ctseg.nii` — group average normalised CT (CTseg)

Publication figures and tables are saved to `manuscript/figures/`:
- `fig_*.png`, `fig_*.pdf` — publication figures
- `table_*.tex` — LaTeX tables

## File Structure

```
experiments/
  run_all.m              Top-level orchestrator
  config.m               Paths and settings (edit this)
  README.md              This file
  steps/                 Processing and metric computation
  analysis/              Statistical tests, figures, tables
  utils/                 Shared helper functions
  atlas/                 Atlas creation and grid searches
  preprocessing/         Data conversion and preprocessing
  review/                Interactive image/segmentation review tools
```
