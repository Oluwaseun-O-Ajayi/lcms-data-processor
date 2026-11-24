"""
Converter cross-check: the same raw file converted by ThermoRawFileParser and by MSConvert,
quantified with the same targets and settings, compared glycan by glycan.

Inputs are the `_glycans.csv` files written by lcms-batch for each conversion and the
`_fld_per_glycan.csv` file of the same run (fluorescence reference).

Usage:
    python converter_check.py --trfp A_glycans.csv --msconvert B_glycans.csv \
        --fld run_fld_per_glycan.csv --out data/converter_check
"""

import argparse
from pathlib import Path

import pandas as pd

KEYS = ['ppm', 'width', 'baseline_window']
PRIMARY = (10.0, 5, 0.1)


def load(path):
    d = pd.read_csv(path)
    d = d[KEYS + ['glycan', 'peak_area']].copy()
    d['pct'] = 100 * d['peak_area'] / d.groupby(KEYS)['peak_area'].transform('sum')
    return d


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument('--trfp', required=True)
    p.add_argument('--msconvert', required=True)
    p.add_argument('--fld', required=True)
    p.add_argument('--out', required=True, help='output prefix')
    a = p.parse_args()

    t, m = load(a.trfp), load(a.msconvert)
    both = t.merge(m, on=KEYS + ['glycan'], suffixes=('_trfp', '_msconvert'), how='outer')
    if both.isna().any().any():
        raise SystemExit('Glycans or settings differ between the two conversions:\n'
                         + both[both.isna().any(axis=1)].to_string(index=False))
    both['diff_pct'] = both['pct_msconvert'] - both['pct_trfp']
    both['area_ratio'] = both['peak_area_msconvert'] / both['peak_area_trfp']

    per_setting = both.groupby(KEYS).agg(
        n_glycans=('glycan', 'size'),
        mad_pct=('diff_pct', lambda x: x.abs().mean()),
        max_abs_diff_pct=('diff_pct', lambda x: x.abs().max()),
        area_ratio_min=('area_ratio', 'min'),
        area_ratio_max=('area_ratio', 'max')).reset_index()

    fld = pd.read_csv(a.fld)[['glycan', 'fld_area']]
    prim = both[(both['ppm'] == PRIMARY[0]) & (both['width'] == PRIMARY[1])
                & (both['baseline_window'] == PRIMARY[2])].merge(fld, on='glycan', how='left')
    if prim.empty or prim['fld_area'].isna().any():
        raise SystemExit('Primary setting or fluorescence values missing for some glycans')
    prim['fld_pct'] = 100 * prim['fld_area'] / prim['fld_area'].sum()
    prim['trfp_minus_fld'] = prim['pct_trfp'] - prim['fld_pct']
    prim['msconvert_minus_fld'] = prim['pct_msconvert'] - prim['fld_pct']

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    both.to_csv(f'{out}_all_settings.csv', index=False)
    per_setting.to_csv(f'{out}_per_setting.csv', index=False)
    cols = ['glycan', 'pct_trfp', 'pct_msconvert', 'diff_pct', 'area_ratio', 'fld_pct',
            'trfp_minus_fld', 'msconvert_minus_fld']
    prim[cols].to_csv(f'{out}_primary.csv', index=False)

    pd.set_option('display.width', 200)
    print('== Per setting (MSConvert minus ThermoRawFileParser, percentage points)')
    print(per_setting.round(4).to_string(index=False))
    print(f'\n== Primary setting ppm {PRIMARY[0]:g}, width {PRIMARY[1]}, baseline window {PRIMARY[2]:g}')
    print(prim[cols].round(3).to_string(index=False))
    print(f"\nMS vs fluorescence MAD: ThermoRawFileParser {prim['trfp_minus_fld'].abs().mean():.3f}, "
          f"MSConvert {prim['msconvert_minus_fld'].abs().mean():.3f}")
    print(f"Converter difference MAD (all {len(per_setting)} settings): "
          f"{per_setting['mad_pct'].min():.4f} to {per_setting['mad_pct'].max():.4f}")


if __name__ == '__main__':
    main()
