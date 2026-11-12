"""
Measured mass error of each target ion.

For each target (name, mz, rt_min, rt_max), finds the scan with the most signal within
±search_ppm inside the retention window, takes the spectra within ±apex_window minutes of
that scan, and reports the intensity-weighted mean m/z of the centroids within ±search_ppm
as a ppm error from the theoretical value.

Usage:
    python examples/mass_errors.py --mzml run.mzML --targets targets.csv --out mass_errors.csv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lcms_data_processor import MzMLReader  # noqa: E402


def mass_errors(reader, targets, search_ppm=20.0, apex_window=0.1):
    spectra = [s for s in reader.experiment if s.getMSLevel() == 1]
    times = np.array([s.getRT() / 60.0 for s in spectra])
    rows = []
    for _, t in targets.iterrows():
        tol = t['mz'] * search_ppm * 1e-6
        idx = np.where((times >= t['rt_min']) & (times <= t['rt_max']))[0]
        signal = []
        for i in idx:
            m, inten = spectra[i].get_peaks()
            signal.append(inten[np.abs(m - t['mz']) <= tol].sum())
        if not len(idx) or max(signal) <= 0:
            rows.append({**t.to_dict(), 'apex_rt': np.nan, 'observed_mz': np.nan, 'error_ppm': np.nan,
                         'n_centroids': 0, 'summed_intensity': 0.0})
            continue
        apex = times[idx[int(np.argmax(signal))]]
        mz_all, w_all = [], []
        for i in np.where(np.abs(times - apex) <= apex_window)[0]:
            m, inten = spectra[i].get_peaks()
            keep = np.abs(m - t['mz']) <= tol
            mz_all.append(m[keep]); w_all.append(inten[keep])
        mz_all, w_all = np.concatenate(mz_all), np.concatenate(w_all).astype(float)
        obs = float(np.sum(mz_all * w_all) / w_all.sum())
        rows.append({**t.to_dict(), 'apex_rt': apex, 'observed_mz': obs,
                     'error_ppm': (obs - t['mz']) / t['mz'] * 1e6,
                     'n_centroids': len(mz_all), 'summed_intensity': float(w_all.sum())})
    return pd.DataFrame(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--mzml', required=True)
    p.add_argument('--targets', required=True)
    p.add_argument('--search-ppm', type=float, default=20.0)
    p.add_argument('--apex-window', type=float, default=0.1, help='minutes either side of the apex')
    p.add_argument('--out', required=True)
    a = p.parse_args(argv)
    res = mass_errors(MzMLReader(a.mzml), pd.read_csv(a.targets), a.search_ppm, a.apex_window)
    res.to_csv(a.out, index=False)
    cols = [c for c in ('name', 'mz', 'apex_rt', 'observed_mz', 'error_ppm', 'n_centroids') if c in res]
    print(res[cols].to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == '__main__':
    main()
