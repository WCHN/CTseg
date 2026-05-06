#!/usr/bin/env python3
"""Deformable registration of CT to MR using Elastix.

Replicates the SynthRAD2025 challenge preprocessing (Section 2.4.6 of
Thummerer et al. 2025) using the same Elastix B-spline parameter file.

For each subject, registers pp_ct to pp_mr (fixed) and saves pp_ct_def.nii.

Usage:
    python deform_ct_to_mr.py              # all subjects
    python deform_ct_to_mr.py --test       # test mode (first N subjects)
    python deform_ct_to_mr.py 1HNA004      # single subject
    python deform_ct_to_mr.py --force       # overwrite existing pp_ct_def.nii

Requires: SimpleITK (pip install SimpleITK) in the ctseg conda environment.
"""

import argparse
import os
import sys
import tempfile
import shutil

import SimpleITK as sitk

# ---------------------------------------------------------------------------
# Configuration (matching experiments/config.m)
# ---------------------------------------------------------------------------
DATA_DIR = r'C:\Users\mbrudfors\Data\HN'
# Convert Windows path for WSL if needed
if not os.path.exists(DATA_DIR):
    DATA_DIR = '/mnt/c/Users/mbrudfors/Data/HN'

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARAM_FILE = os.path.join(SCRIPT_DIR, 'configs', 'param_def_mr_HN.txt')

TEST_N_SUBJECTS = 1

EXCLUDE = {
    '1HNA098','1HNA099','1HNA100','1HNA102','1HNA103','1HNA104',
    '1HNA105','1HNA106','1HNA107','1HNA108','1HNA109','1HNA110',
    '1HNA113','1HNA115','1HNA116','1HNA117','1HNA119','1HNA120',
    '1HNA121','1HNA124','1HNA126','1HNA129','1HNA130','1HNA132',
    '1HNA133','1HNA135','1HNA136','1HNA138','1HNA139','1HNA141',
    '1HNA142','1HNA143',
    '1HNC001','1HNC002','1HNC003','1HNC004','1HNC005','1HNC007',
    '1HNC008','1HNC012','1HNC014','1HNC017','1HNC019','1HNC020',
    '1HNC021','1HNC022','1HNC023','1HNC025','1HNC029','1HNC031',
    '1HNC035','1HNC036','1HNC037','1HNC038','1HNC039','1HNC040',
    '1HNC043','1HNC045','1HNC046','1HNC050','1HNC061','1HNC066',
    '1HNC067','1HNC068','1HNC071','1HNC072','1HNC073','1HNC076',
    '1HNC082','1HNC083','1HNC084','1HNC085','1HNC087','1HNC088',
    '1HNC094','1HNC098','1HNC099','1HNC101','1HNC102','1HNC103',
    '1HNC104','1HNC105','1HNC107','1HNC109','1HNC110','1HNC111',
    '1HNC112','1HNC114','1HNC117','1HNC118','1HNC120','1HNC121',
    '1HNC124','1HNC125','1HNC127','1HNC128','1HNC130',
}


def get_subjects(test_mode=False):
    """Return sorted list of subject directory names, excluding flagged ones."""
    dirs = sorted([
        d for d in os.listdir(DATA_DIR)
        if os.path.isdir(os.path.join(DATA_DIR, d))
        and not d.startswith('.')
        and d not in ('overviews', 'results')
        and d not in EXCLUDE
    ])
    if test_mode:
        dirs = dirs[:TEST_N_SUBJECTS]
    return dirs


def find_nii(subj_dir, name):
    """Find a NIfTI file, handling both .nii and .nii.gz."""
    for ext in ('.nii', '.nii.gz'):
        path = os.path.join(subj_dir, name + ext)
        if os.path.exists(path):
            return path
    return None


def deformable_registration(fixed, moving, parameter_file):
    """Run Elastix deformable registration (replicating SynthRAD2025 Stage 2).

    Args:
        fixed: SimpleITK image (MR — reference anatomy)
        moving: SimpleITK image (CT — to be deformed)
        parameter_file: path to Elastix parameter file

    Returns:
        SimpleITK image of the deformed moving image
    """
    temp_dir = tempfile.mkdtemp()
    orig_dir = os.getcwd()
    os.chdir(temp_dir)

    try:
        param = sitk.ReadParameterFile(parameter_file)
        elastix = sitk.ElastixImageFilter()
        elastix.SetParameterMap(param)
        elastix.SetFixedImage(fixed)
        elastix.SetMovingImage(moving)
        elastix.LogToConsoleOn()
        elastix.LogToFileOff()
        elastix.Execute()
        result = elastix.GetResultImage()
    finally:
        os.chdir(orig_dir)
        shutil.rmtree(temp_dir, ignore_errors=True)

    return result


def process_subject(subj_name, force=False):
    """Register pp_ct to pp_mr for one subject."""
    subj_dir = os.path.join(DATA_DIR, subj_name)
    out_file = os.path.join(subj_dir, 'pp_ct_def.nii')

    if os.path.exists(out_file) and not force:
        print(f'  {subj_name} — SKIP (pp_ct_def.nii exists)')
        return True

    mr_file = find_nii(subj_dir, 'pp_mr')
    ct_file = find_nii(subj_dir, 'pp_ct')
    if mr_file is None or ct_file is None:
        print(f'  {subj_name} — SKIP (missing pp_mr or pp_ct)')
        return False

    print(f'  {subj_name} — registering...')
    fixed = sitk.ReadImage(mr_file, sitk.sitkFloat32)
    moving = sitk.ReadImage(ct_file, sitk.sitkFloat32)

    result = deformable_registration(fixed, moving, PARAM_FILE)

    # Cast to int16 (matching SynthRAD2025 convention) and save
    result = sitk.Cast(result, sitk.sitkInt16)
    sitk.WriteImage(result, out_file)
    print(f'  {subj_name} — saved pp_ct_def.nii')
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Deformable registration of CT to MR (Elastix)')
    parser.add_argument('subject', nargs='?', default=None,
                        help='Single subject ID (e.g. 1HNA004)')
    parser.add_argument('--test', action='store_true',
                        help=f'Test mode: process first {TEST_N_SUBJECTS} subject(s)')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite existing pp_ct_def.nii')
    args = parser.parse_args()

    if not os.path.exists(PARAM_FILE):
        sys.exit(f'Parameter file not found: {PARAM_FILE}')

    if args.subject:
        subjects = [args.subject]
    else:
        subjects = get_subjects(test_mode=args.test)

    n = len(subjects)
    print(f'Deformable registration: {n} subject(s)')
    if args.test:
        print(f'  [TEST MODE]')

    n_ok = 0
    for i, subj in enumerate(subjects):
        print(f'[{i+1}/{n}]', end='')
        if process_subject(subj, force=args.force):
            n_ok += 1

    print(f'\nDone. {n_ok}/{n} subjects processed.')


if __name__ == '__main__':
    main()
