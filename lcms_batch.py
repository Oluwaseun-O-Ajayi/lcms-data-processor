"""
Batch processing of mzML files over a grid of processing settings.

For every file, extraction window (ppm) and peak-detection setting, the script extracts
one chromatogram per target m/z, finds the largest peak (inside the target's retention
time window when one is given), integrates it, and writes one row per
file x setting x target. A second table sums charge states per glycan and gives
relative abundances (percent of the summed area within each file and setting).

Usage:
    lcms-batch --mzml data/*.mzML --targets targets.csv --ppm 5 10 20 50 \\
               --baseline-window 0.05 0.1 0.2 --width 3 5 --out results/sweep

    # or build targets from glycan compositions:
    lcms-batch --mzml data/*.mzML --glycans fetuin --label 2-AB --charges 1 2 3 ...

targets.csv columns: name, mz, and optionally glycan, rt_min, rt_max (minutes).
"""

import argparse
import itertools
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from lcms_data_processor import MzMLReader, PeakDetector
import glycan_targets as gt

PRESETS = {'igg': gt.IGG_CANDIDATES, 'fetuin': gt.FETUIN_CANDIDATES,
           'all': {**gt.IGG_CANDIDATES, **gt.FETUIN_CANDIDATES}}


def load_targets(args):
    if args.targets:
        t = pd.read_csv(args.targets)
        if not {'name', 'mz'} <= set(t.columns):
            raise SystemExit("targets file needs columns: name, mz (optional: glycan, rt_min, rt_max)")
        if 'glycan' not in t:
            t['glycan'] = t['name']
    else:
        comps = PRESETS[args.glycans]
        t = gt.target_table(comps, label=args.label, charges=args.charges)
        t['name'] = t['glycan'] + '_' + t['charge'].astype(str) + '+'
    for col in ('rt_min', 'rt_max'):
        if col not in t:
            t[col] = np.nan
    return t.reset_index(drop=True)


def best_peak(times, intensity, rt_min, rt_max, width, baseline_window, prominence, mode='largest'):
    """
    Integrated signal for one target.

    mode='largest' returns the largest peak (inside the RT window if given).
    mode='sum' sums every peak inside the RT window, for example partly resolved
    linkage isomers of one composition; retention time, height, FWHM and S/N are
    those of the largest peak, and 'peak_rts' lists the apex of every peak summed.
    """
    chrom = pd.DataFrame({'time': times, 'intensity': intensity})
    if np.nanmax(intensity) <= 0:
        return None, 0
    det = PeakDetector(chrom)
    peaks = det.detect_peaks(prominence=prominence, width=width, baseline_window=baseline_window)
    if peaks.empty:
        return None, 0
    det.integrate_peaks()
    det.calculate_snr()
    peaks = det.peaks
    if np.isfinite(rt_min):
        peaks = peaks[peaks['retention_time'] >= rt_min]
    if np.isfinite(rt_max):
        peaks = peaks[peaks['retention_time'] <= rt_max]
    if peaks.empty:
        return None, len(det.peaks)
    best = peaks.loc[peaks['peak_area'].idxmax()].copy()
    if mode == 'sum':
        best['peak_area'] = peaks['peak_area'].sum()
        best['n_summed'] = len(peaks)
        best['peak_rts'] = ';'.join(f"{t:.4f}" for t in sorted(peaks['retention_time']))
    else:
        best['n_summed'] = 1
        best['peak_rts'] = f"{best['retention_time']:.4f}"
    return best, len(det.peaks)


def run(args):
    targets = load_targets(args)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    grid = list(itertools.product(args.ppm, args.width, args.baseline_window))
    for f in args.mzml:
        t0 = time.time()
        reader = MzMLReader(f)
        print(f"{reader.metadata['filename']}: {reader.metadata['num_spectra']} spectra; "
              f"{len(targets)} targets x {len(grid)} settings", flush=True)
        for ppm in args.ppm:
            times, mat = reader.xic_matrix(targets['mz'].to_numpy(), ppm=ppm, ms_level=args.ms_level)
            for width, bw in itertools.product(args.width, args.baseline_window):
                for j, tg in targets.iterrows():
                    peak, n_peaks = best_peak(times, mat[:, j], tg['rt_min'], tg['rt_max'],
                                              width, bw, args.prominence, args.peak_mode)
                    rows.append({
                        'file': Path(f).name, 'ppm': ppm, 'width': width, 'baseline_window': bw,
                        'target': tg['name'], 'glycan': tg['glycan'], 'mz': tg['mz'],
                        'charge': tg['charge'] if 'charge' in tg else np.nan,
                        'n_peaks': n_peaks,
                        'retention_time': None if peak is None else peak['retention_time'],
                        'peak_area': 0.0 if peak is None else peak['peak_area'],
                        'peak_height': None if peak is None else peak['peak_height'],
                        'fwhm': None if peak is None else peak['fwhm'],
                        'snr': None if peak is None else peak['snr'],
                        'n_summed': 0 if peak is None else int(peak['n_summed']),
                        'peak_rts': '' if peak is None else peak['peak_rts'],
                        'rt_min': tg['rt_min'], 'rt_max': tg['rt_max'],
                        'peak_mode': args.peak_mode,
                    })
        print(f"  done in {time.time() - t0:.1f} s", flush=True)

    long = pd.DataFrame(rows)
    keys = ['file', 'ppm', 'width', 'baseline_window']
    per_glycan = long.groupby(keys + ['glycan'], as_index=False)['peak_area'].sum()
    total = per_glycan.groupby(keys)['peak_area'].transform('sum')
    per_glycan['relative_abundance_pct'] = np.where(total > 0, 100 * per_glycan['peak_area'] / total, np.nan)

    long.to_csv(f"{out}_targets.csv", index=False)
    per_glycan.to_csv(f"{out}_glycans.csv", index=False)
    targets.to_csv(f"{out}_target_list.csv", index=False)
    print(f"Wrote {out}_targets.csv, {out}_glycans.csv and {out}_target_list.csv")
    return long, per_glycan


def build_parser():
    p = argparse.ArgumentParser(prog='lcms-batch', description=__doc__.split('\n\n')[0])
    p.add_argument('--mzml', nargs='+', required=True, help='mzML files')
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument('--targets', help='CSV with name, mz (optional glycan, rt_min, rt_max)')
    src.add_argument('--glycans', choices=sorted(PRESETS), help='built-in candidate glycan compositions')
    p.add_argument('--label', default='2-AB', help='reducing-end label for --glycans (default 2-AB)')
    p.add_argument('--charges', type=int, nargs='+', default=[1, 2], help='charge states for --glycans')
    p.add_argument('--ppm', type=float, nargs='+', default=[10.0], help='extraction half-widths (ppm)')
    p.add_argument('--width', type=int, nargs='+', default=[5], help='minimum peak width (points)')
    p.add_argument('--baseline-window', type=float, nargs='+', default=[0.1],
                   help='baseline window as a fraction of the run')
    p.add_argument('--prominence', type=float, default=0.1, help='minimum prominence (fraction of max)')
    p.add_argument('--peak-mode', choices=['largest', 'sum'], default='largest',
                   help="'largest' peak, or 'sum' of all peaks inside each target's RT window")
    p.add_argument('--ms-level', type=int, default=1)
    p.add_argument('--out', default='results/sweep', help='output prefix')
    return p


def main(argv=None):
    run(build_parser().parse_args(argv))


if __name__ == '__main__':
    sys.exit(main())
