"""
Summarise the procainamide quantitation: the fetuin sample-to-sample comparison, the transferrin
profile, mass errors and charge-state distribution.

Usage:
    python cluster/procA_summary.py --results <folder> --out <prefix>
"""

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd

NAMES = {'Fet1_052124': 'FET-ProA-S1', 'Fet2_052124': 'FET-ProA-S2', 'HTferin052124': 'TFN-ProA'}


def sample_name(file):
    for stem, name in NAMES.items():
        if stem in file or name in file:
            return name
    return Path(file).stem


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--results', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args(argv)
    r = a.results

    # Relative abundance per glycan, every ppm setting in the run.
    g = pd.concat([pd.read_csv(f) for f in glob.glob(f"{r}/*_glycans.csv")], ignore_index=True)
    g['sample'] = g['file'].map(sample_name)

    fet = g[g['sample'].str.startswith('FET')]
    wide = fet.pivot_table(index=['ppm', 'glycan'], columns='sample', values='relative_abundance_pct').reset_index()
    i1, i2 = 'FET-ProA-S1', 'FET-ProA-S2'
    wide['difference_pct_points'] = wide[i2] - wide[i1]
    wide.round(3).to_csv(f"{a.out}_fetuin_comparison.csv", index=False)
    print('Fetuin, sample-to-sample comparison (relative abundance, %):')
    print(wide.round(2).to_string(index=False))

    tfn = g[g['sample'] == 'TFN-ProA'][['ppm', 'glycan', 'peak_area', 'relative_abundance_pct']]
    tfn.round(3).to_csv(f"{a.out}_transferrin_profile.csv", index=False)
    print('\nTransferrin profile (relative abundance, %):')
    print(tfn.round(2).to_string(index=False))

    # Charge-state distribution of each glycan.
    t = pd.concat([pd.read_csv(f) for f in glob.glob(f"{r}/*_targets.csv")
                   if not f.endswith('_target_list.csv')], ignore_index=True)
    t['sample'] = t['file'].map(sample_name)
    t = t[t['ppm'] == t['ppm'].max()]  # widest window: least affected by mass-calibration offset
    cs = t.pivot_table(index=['sample', 'glycan'], columns='charge', values='peak_area', aggfunc='sum').fillna(0)
    cs = cs.div(cs.sum(axis=1).replace(0, np.nan), axis=0).mul(100).round(1)
    cs.columns = [f'pct_{int(c)}+' for c in cs.columns]
    cs.reset_index().to_csv(f"{a.out}_charge_states.csv", index=False)

    # Mass errors.
    me = []
    for f in glob.glob(f"{r}/*_mass_errors.csv"):
        d = pd.read_csv(f)
        d['sample'] = sample_name(Path(f).name)
        me.append(d)
    me = pd.concat(me, ignore_index=True)
    me.to_csv(f"{a.out}_mass_errors_all.csv", index=False)
    ok = me.dropna(subset=['error_ppm'])
    s = ok.groupby('sample')['error_ppm'].agg(n='count', median='median', min='min', max='max',
                                              mean_abs=lambda x: x.abs().mean()).round(2)
    s.reset_index().to_csv(f"{a.out}_mass_errors_summary.csv", index=False)
    print('\nMass error (ppm) at the apex of each target ion:')
    print(s.to_string())


if __name__ == '__main__':
    main()
