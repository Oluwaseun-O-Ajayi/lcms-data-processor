"""
Suggest glycan compositions for an unassigned chromatographic peak.

Sums MS1 spectra over a retention-time window, lists the most intense m/z values, and
matches each against every composition in a search space (Hex, HexNAc, dHex, NeuAc
counts) with the chosen reducing-end label and charge states. Matches are candidates
only: confirm them with retention behaviour, isotope pattern and MS/MS.

Usage:
    lcms-glycan-search --mzml run.mzML --rt 26.6 26.95 --label 2-AB --charges 1 2 3 --ppm 10
"""

import argparse
import itertools
import sys

import numpy as np
import pandas as pd

import glycan_targets as gt
from lcms_data_processor import MzMLReader


def summed_spectrum(reader, rt_lo, rt_hi, ms_level=1, bin_ppm=5.0):
    """Sum spectra in [rt_lo, rt_hi] minutes; centroids are merged within bin_ppm."""
    mzs, ints = [], []
    for s in reader.experiment:
        if s.getMSLevel() != ms_level:
            continue
        t = s.getRT() / 60.0
        if rt_lo <= t <= rt_hi:
            m, i = s.get_peaks()
            mzs.append(m)
            ints.append(i)
    if not mzs:
        raise ValueError(f"No MS{ms_level} spectra between {rt_lo} and {rt_hi} min")
    m = np.concatenate(mzs).astype(float)
    i = np.concatenate(ints).astype(float)
    order = np.argsort(m)
    m, i = m[order], i[order]
    # Merge neighbouring centroids closer than bin_ppm (intensity-weighted m/z).
    groups = np.concatenate(([0], np.cumsum(np.diff(m) / m[1:] * 1e6 > bin_ppm)))
    df = pd.DataFrame({'g': groups, 'mz': m, 'i': i, 'mi': m * i})
    agg = df.groupby('g').agg(intensity=('i', 'sum'), mzi=('mi', 'sum'))
    return pd.DataFrame({'mz': agg['mzi'] / agg['intensity'], 'intensity': agg['intensity']})


def composition_space(hex_range=(3, 9), hexnac_range=(2, 7), dhex_range=(0, 2), neuac_range=(0, 5)):
    for h, n, f, s in itertools.product(range(hex_range[0], hex_range[1] + 1),
                                        range(hexnac_range[0], hexnac_range[1] + 1),
                                        range(dhex_range[0], dhex_range[1] + 1),
                                        range(neuac_range[0], neuac_range[1] + 1)):
        # N-glycan core needs HexNAc >= 2 and Hex >= 3; sialic acids need antennae.
        if s > max(n - 2, 0) + (1 if h > 3 else 0) or (n - 2) > h:
            continue
        yield f"N{n}H{h}" + (f"F{f}" if f else "") + (f"S{s}" if s else "")


def search(spectrum, label='2-AB', charges=(1, 2, 3), ppm=10.0, top=15, adduct='H'):
    table = []
    comps = list(composition_space())
    theo = [(c, z, gt.glycan_mz(c, z, label, adduct)) for c in comps for z in charges]
    theo_mz = np.array([t[2] for t in theo])
    for _, pk in spectrum.nlargest(top, 'intensity').iterrows():
        err = (theo_mz - pk['mz']) / theo_mz * 1e6
        hits = np.where(np.abs(err) <= ppm)[0]
        if len(hits) == 0:
            table.append({'mz': pk['mz'], 'intensity': pk['intensity'], 'composition': '', 'charge': None,
                          'theoretical_mz': None, 'error_ppm': None})
        for k in hits[np.argsort(np.abs(err[hits]))]:
            c, z, mz = theo[k]
            table.append({'mz': pk['mz'], 'intensity': pk['intensity'], 'composition': c, 'charge': z,
                          'theoretical_mz': mz, 'error_ppm': err[k]})
    return pd.DataFrame(table)


def main(argv=None):
    p = argparse.ArgumentParser(prog='lcms-glycan-search', description=__doc__.split('\n\n')[0])
    p.add_argument('--mzml', required=True)
    p.add_argument('--rt', type=float, nargs=2, required=True, metavar=('START', 'END'), help='minutes')
    p.add_argument('--label', default='2-AB')
    p.add_argument('--charges', type=int, nargs='+', default=[1, 2, 3])
    p.add_argument('--ppm', type=float, default=10.0)
    p.add_argument('--top', type=int, default=15, help='number of most intense m/z values to search')
    p.add_argument('--out', help='optional CSV output')
    a = p.parse_args(argv)
    spec = summed_spectrum(MzMLReader(a.mzml), *a.rt)
    res = search(spec, a.label, a.charges, a.ppm, a.top)
    if a.out:
        res.to_csv(a.out, index=False)
    with pd.option_context('display.width', 140):
        print(res.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == '__main__':
    sys.exit(main())
