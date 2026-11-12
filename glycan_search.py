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


ISOTOPE_SPACING = 1.003355  # 13C - 12C


def survey(reader, label='2-AB', charges=(2, 3, 4), ppm=10.0, step=0.25, min_rel=0.002,
           isotope_ratio=(0.2, 3.0), max_prev_ratio=0.3, space=None):
    """
    Scan the run in `step`-minute slices. In each summed spectrum, a composition at charge z
    is a hit when its monoisotopic m/z and its M+1 isotope (+1.003355/z) are both found
    within `ppm`, the M+1/M intensity ratio lies within `isotope_ratio` (glycans of this size
    have M+1/M near 1), the signal one isotope spacing below is under `max_prev_ratio` of the
    monoisotopic intensity, and the monoisotopic intensity is at least `min_rel` of the most
    intense peak in that slice. Returns one row per composition, charge and slice.
    """
    comps = list(space or composition_space())
    theo = [(c, z, gt.glycan_mz(c, z, label)) for c in comps for z in charges]
    mz_mono = np.array([t[2] for t in theo])
    spectra = [s for s in reader.experiment if s.getMSLevel() == 1]
    times = np.array([s.getRT() / 60.0 for s in spectra])
    rows = []
    for lo in np.arange(times.min(), times.max(), step):
        idx = np.where((times >= lo) & (times < lo + step))[0]
        if len(idx) == 0:
            continue
        m = np.concatenate([spectra[i].get_peaks()[0] for i in idx]).astype(float)
        inten = np.concatenate([spectra[i].get_peaks()[1] for i in idx]).astype(float)
        if len(m) == 0:
            continue
        order = np.argsort(m)
        m, csum = m[order], np.concatenate(([0.0], np.cumsum(inten[order])))
        top = inten.max() * len(idx)  # scale so the threshold refers to summed intensity

        def window_sum(targets):
            tol = targets * ppm * 1e-6
            return csum[np.searchsorted(m, targets + tol, 'right')] - csum[np.searchsorted(m, targets - tol, 'left')]

        z_arr = np.array([t[1] for t in theo])
        mono = window_sum(mz_mono)
        iso = window_sum(mz_mono + ISOTOPE_SPACING / z_arr)
        prev = window_sum(mz_mono - ISOTOPE_SPACING / z_arr)
        ratio = np.divide(iso, mono, out=np.zeros_like(iso), where=mono > 0)
        # A monoisotopic peak has no strong peak one isotope spacing below it. This rejects
        # compositions that land on another glycan's isotope peak (for example two dHex vs
        # one NeuAc, which differ by 1.0204 Da, close to the 13C spacing).
        ok = ((mono >= min_rel * top) & (ratio >= isotope_ratio[0]) & (ratio <= isotope_ratio[1])
              & (prev < max_prev_ratio * mono))
        for k in np.where(ok)[0]:
            c, z, mz = theo[k]
            rows.append({'composition': c, 'charge': z, 'mz': mz, 'rt_start': lo, 'rt_end': lo + step,
                         'intensity': mono[k], 'isotope_ratio': ratio[k]})
    return pd.DataFrame(rows)


def summarize_survey(hits):
    """One row per composition: charges seen, apex slice, total and maximum intensity."""
    if hits.empty:
        return hits
    out = []
    for c, g in hits.groupby('composition'):
        best = g.loc[g['intensity'].idxmax()]
        out.append({'composition': c, 'charges_seen': ','.join(str(z) for z in sorted(g['charge'].unique())),
                    'n_slices': g[['rt_start']].drop_duplicates().shape[0],
                    'apex_rt': 0.5 * (best['rt_start'] + best['rt_end']), 'max_intensity': best['intensity'],
                    'total_intensity': g['intensity'].sum()})
    return pd.DataFrame(out).sort_values('total_intensity', ascending=False)


def main(argv=None):
    p = argparse.ArgumentParser(prog='lcms-glycan-search', description=__doc__.split('\n\n')[0])
    p.add_argument('--mzml', required=True)
    p.add_argument('--rt', type=float, nargs=2, metavar=('START', 'END'), help='minutes')
    p.add_argument('--survey', action='store_true', help='scan the whole run (see module docstring)')
    p.add_argument('--step', type=float, default=0.25, help='survey slice width (min)')
    p.add_argument('--min-rel', type=float, default=0.002, help='survey intensity threshold (fraction)')
    p.add_argument('--fld', help='survey: also list fluorescence peaks of this chromatogram (e.g. FLD1)')
    p.add_argument('--label', default='2-AB')
    p.add_argument('--charges', type=int, nargs='+', default=[1, 2, 3])
    p.add_argument('--ppm', type=float, default=10.0)
    p.add_argument('--top', type=int, default=15, help='number of most intense m/z values to search')
    p.add_argument('--out', help='optional CSV output')
    a = p.parse_args(argv)
    if a.survey:
        reader = MzMLReader(a.mzml)
        hits = survey(reader, a.label, tuple(a.charges), a.ppm, a.step, a.min_rel)
        summary = summarize_survey(hits)
        if a.out:
            hits.to_csv(a.out.replace('.csv', '_hits.csv'), index=False)
            summary.to_csv(a.out, index=False)
        with pd.option_context('display.width', 160, 'display.max_rows', 200):
            print(summary.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
        if a.fld:
            from lcms_data_processor import PeakDetector
            det = PeakDetector(reader.chromatogram(a.fld))
            det.detect_peaks(prominence=0.005, width=5)
            det.integrate_peaks()
            fp = det.peaks[['retention_time', 'peak_height', 'peak_area', 'fwhm']]
            if a.out:
                fp.to_csv(a.out.replace('.csv', '_fld_peaks.csv'), index=False)
            print("\nFluorescence peaks over the whole run:")
            print(fp.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
        return
    if not a.rt:
        p.error('--rt is required unless --survey is given')
    spec = summed_spectrum(MzMLReader(a.mzml), *a.rt)
    res = search(spec, a.label, a.charges, a.ppm, a.top)
    if a.out:
        res.to_csv(a.out, index=False)
    with pd.option_context('display.width', 140):
        print(res.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == '__main__':
    sys.exit(main())
