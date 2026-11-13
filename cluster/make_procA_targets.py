"""
Write the procainamide target lists for fetuin and human transferrin.

Retention windows were read from the untargeted survey of the same runs (FET-ProA-S1/I2 and
TFN-ProA-P1) and checked against the retention times recorded for this method
(A3G3 26.8, A3G3S1 29.9, A3G3S2 32.2, A3G3S3 34.8 min). Each window covers the apex of the
dominant charge state in every run plus the neighbouring isomer peaks, and stops short of the
apex of the next larger glycan, so that in-source fragments of a larger glycan (for example
A2G2 ions at the A3G3 retention time) are not counted as the smaller glycan.

Usage:
    python cluster/make_procA_targets.py
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from glycan_targets import glycan_mz  # noqa: E402

# name: (composition, rt_min, rt_max, charges)
FETUIN = {
    'A2G2':   ('N4H5',   24.50, 26.25, (2, 3)),
    'A3G3':   ('N5H6',   26.50, 27.90, (2, 3)),
    'A2G2S1': ('N4H5S1', 27.90, 29.10, (2, 3)),
    'A3G3S1': ('N5H6S1', 29.10, 31.00, (2, 3)),
    'A2G2S2': ('N4H5S2', 31.00, 31.90, (2, 3)),
    'A3G3S2': ('N5H6S2', 31.40, 33.90, (2, 3)),
    'A3G3S3': ('N5H6S3', 33.90, 36.50, (2, 3, 4)),
    'A3G3S4': ('N5H6S4', 36.50, 39.50, (3, 4)),
}

TRANSFERRIN = {
    'A2G1':    ('N4H4',     23.25, 24.50, (2, 3)),
    'A2G2':    ('N4H5',     24.50, 25.60, (2, 3)),
    'FA2G2':   ('N4H5F1',   25.60, 26.60, (2, 3)),
    'A3G3':    ('N5H6',     26.60, 27.90, (2, 3)),
    'A2G1S1':  ('N4H4S1',   26.75, 27.90, (2, 3)),
    'A2G2S1':  ('N4H5S1',   27.90, 29.10, (2, 3)),
    'FA2G2S1': ('N4H5F1S1', 28.60, 29.80, (2, 3)),
    'A3G3S1':  ('N5H6S1',   29.10, 31.00, (2, 3)),
    'A2G2S2':  ('N4H5S2',   31.00, 32.40, (2, 3)),
    'FA2G2S2': ('N4H5F1S2', 31.50, 32.75, (2, 3)),
    'A3G3S2':  ('N5H6S2',   32.40, 33.90, (2, 3)),
}


def table(spec, label='procainamide'):
    rows = []
    for glycan, (comp, lo, hi, charges) in spec.items():
        for z in charges:
            rows.append({'name': f'{glycan}_{z}+', 'glycan': glycan, 'composition': comp, 'charge': z,
                         'mz': round(glycan_mz(comp, z, label), 5), 'rt_min': lo, 'rt_max': hi})
    return pd.DataFrame(rows)


if __name__ == '__main__':
    here = Path(__file__).resolve().parent
    table(FETUIN).to_csv(here / 'procA_fetuin_targets.csv', index=False)
    table(TRANSFERRIN).to_csv(here / 'procA_transferrin_targets.csv', index=False)
    print('Wrote procA_fetuin_targets.csv and procA_transferrin_targets.csv')
