"""
Run manifest and label-chemistry QC for the procainamide LC-IM-MS runs.

For every ion in a run's survey hits (20-40 min), the neutral mass is compared with the base
mass of common fetuin and transferrin N-glycans. The difference identifies the reducing-end
chemistry that produced the ion:

    procainamide           +219.1735 Da   (the intended label)
    procainamide + 1.020   +220.194 Da    (unassigned series)
    procainamide + C2H3N   +260.200 Da    (acetonitrile adduct series)
    2-AB                   +120.0687 Da   (control for carry-over of the other label)

A run passes when its procainamide-labelled signal is at least MIN_SIGNAL and the acetonitrile
series is less than MAX_ACN_SHARE of the labelled signal. Both thresholds were set from the
runs in this manifest and are reported as such; they are not general acceptance limits.

Usage:
    python cluster/procA_qc.py --hits-dir <folder with *_survey_hits.csv> --out <prefix>
"""

import argparse
import glob
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from glycan_targets import glycan_mass  # noqa: E402

PROTON = 1.007276
SHIFTS = {'procainamide': 219.1735, 'procainamide+1.020': 220.1939,
          'procainamide+ACN': 260.2000, '2-AB': 120.0687}
COMPOSITIONS = ['N4H5', 'N4H5S1', 'N4H5S2', 'N5H6', 'N5H6S1', 'N5H6S2', 'N5H6S3', 'N5H6S4',
                'N4H4', 'N5H5', 'N5H5S1', 'N3H4', 'N3H5', 'N3H3', 'N4H3', 'N4H4S1',
                'N4H5F1', 'N4H5F1S1', 'N4H5F1S2']
MIN_SIGNAL = 1e7
MAX_ACN_SHARE = 0.25

# survey file stem, sample name, raw file, analyte, role, note
MANIFEST = [
    ('Fet1_052124', 'FET-ProA-S1', '2024-05-21-15.22.10-Fet1_052124[result].d', 'Bovine fetuin', 'included',
     'Fetuin sample 1 (Fet1); same sequence as S2, injected 8 h earlier'),
    ('Fet2_052124', 'FET-ProA-S2', '2024-05-21-23.22.19-Fet2_052124[result].d', 'Bovine fetuin', 'included',
     'Fetuin sample 2 (Fet2); separate sample from S1'),
    ('HTferin052124', 'TFN-ProA-P1', '2024-05-22-00.30.37-HTferin052124[result].d', 'Human transferrin',
     'included', 'Second glycoprotein, same sequence'),
    ('Dex5_052124', 'DEX-ProA-P1', '2024-05-21-14.13.43-Dex5_052124[result].d', 'Dextran ladder', 'included',
     'Glucose unit calibrant, same sequence'),
    ('Fetuin051724', 'FET-ProA-ACN-I1', '2024-05-17-19.32.29-Fetuin051724[result].d', 'Bovine fetuin',
     'excluded', 'Reconstituted in 50:50 acetonitrile:water and frozen (Ajayi and Aduloju, 2026)'),
    ('Fet1_051724', 'FET-ProA-ACN-I2', '2024-05-30-22.49.43-Fet1_051724[result].d', 'Bovine fetuin',
     'excluded', 'Same aliquot as FET-ProA-ACN-I1, re-injected'),
    ('Dex5_051724', 'DEX-ProA-ACN', '2024-05-17-17.15.39-Dex5_051724[result].d', 'Dextran ladder', 'excluded',
     'Only GU2-GU4 detected'),
    ('Fetuin_12142023', 'FET-ProA-X1', 'OA_12142023-02.10.24-Fetuin[result].d', 'Bovine fetuin', 'excluded',
     'Low labelled signal'),
    ('Fet_2024_01_24', 'FET-ProA-X2', 'OA_2024_01_24_15.49.40-Fet[result].d', 'Bovine fetuin', 'excluded',
     'Low labelled signal; elution order inconsistent with the method'),
    ('BFet052824', 'FET-ProA-X3', '2024-06-11-21.26.01-BFet052824[result].d', 'Bovine fetuin', 'excluded',
     'Low labelled signal'),
    ('Fet-ca', 'FET-ProA-X4', '2024-04-11-20.22.16-Fet-ca[result].d', 'Bovine fetuin', 'excluded',
     'Low labelled signal; acetonitrile series dominant'),
    ('Fet-cb', 'FET-ProA-X5', '2024-04-11-23.46.26-Fet-cb[result].d', 'Bovine fetuin', 'excluded',
     'Low labelled signal; acetonitrile series dominant'),
    ('Fet-cc', 'FET-ProA-X6', '2024-04-12-03.10.33-Fet-cc[result].d', 'Bovine fetuin', 'excluded',
     'Low labelled signal'),
    ('BTferin052124', 'TFN-ProA-X1', '2024-05-22-01.39.06-BTferin052124[result].d', 'Transferrin (BTferin)',
     'excluded', 'Low labelled signal'),
]


def label_shift_totals(hits, rt_lo=20.0, rt_hi=40.0, tol=0.01):
    """Summed intensity of survey ions assigned to each reducing-end chemistry."""
    hits = hits[(hits.rt_start >= rt_lo) & (hits.rt_end <= rt_hi)]
    base = {c: glycan_mass(c, 'none') for c in COMPOSITIONS}
    totals = dict.fromkeys(SHIFTS, 0.0)
    for r in hits.itertuples():
        neutral = (r.mz - PROTON) * r.charge
        for b in base.values():
            for name, shift in SHIFTS.items():
                if abs(neutral - b - shift) < tol:
                    totals[name] += r.intensity
    return totals


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--hits-dir', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args(argv)

    rows = []
    for stem, sample, raw, analyte, role, note in MANIFEST:
        row = {'sample': sample, 'analyte': analyte, 'label': 'procainamide', 'raw_file': raw,
               'role': role, 'note': note}
        files = glob.glob(f"{a.hits_dir}/*{stem}_survey_hits.csv")
        if files and not stem.startswith('Dex'):
            t = label_shift_totals(pd.read_csv(files[0]))
            labelled = t['procainamide']
            acn_share = t['procainamide+ACN'] / labelled if labelled else float('inf')
            qc = 'pass' if labelled >= MIN_SIGNAL and acn_share < MAX_ACN_SHARE else 'fail'
            row.update({'procainamide_signal': labelled,
                        'plus1020_signal': t['procainamide+1.020'],
                        'acn_series_signal': t['procainamide+ACN'],
                        '2AB_signal': t['2-AB'],
                        'acn_to_procainamide': round(acn_share, 3), 'qc': qc})
        rows.append(row)
    rows.append({'sample': 'FET-2AB-CRL', 'analyte': 'Bovine fetuin', 'label': '2-AB',
                 'raw_file': 'MASOP00220_AB_CRL_07JAN2021_Fetuin_001 (public)', 'role': 'included',
                 'note': 'Public 2-AB data with fluorescence trace'})
    m = pd.DataFrame(rows)
    m.to_csv(f"{a.out}_qc_summary.csv", index=False)
    m[['sample', 'analyte', 'label', 'raw_file', 'role', 'note']].to_csv(f"{a.out}_manifest.csv", index=False)
    cols = ['sample', 'procainamide_signal', 'acn_to_procainamide', 'qc', 'role']
    print(m[[c for c in cols if c in m]].to_string(index=False))
    bad = m[(m.get('qc') == 'pass') & (m.role == 'excluded') | (m.get('qc') == 'fail') & (m.role == 'included')]
    if len(bad):
        print('\nWARNING: QC result disagrees with manifest role for:', ', '.join(bad['sample']))


if __name__ == '__main__':
    main()
