"""Check whether a non-MS trace in an mzML file behaves like a fluorescence (FLD) detector.

For each mzML file: lists every chromatogram, then, for the chosen detector trace, finds its
peaks, estimates the detector-to-MS delay by cross-correlation with the MS base-peak trace,
and saves an overlay figure (detector above, MS base peak below, same time axis).

A fluorescence channel from a released-glycan HILIC run shows a series of sharp peaks in the
glycan elution window whose apexes line up with MS base-peak apexes after one constant delay.

Usage (compute node, inside the lcms venv):
    python cluster/fld_overlay.py /scratch/$USER/glycan_data/pride/mzml_all/*.mzML \
        --channel Emission --out /scratch/$USER/glycan_data/pride/fld_check
"""

import argparse
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyopenms as oms
from scipy.signal import find_peaks


def load_chromatograms(path):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(path), exp)
    out = {}
    for c in exp.getChromatograms():
        t, y = c.get_peaks()
        out[c.getNativeID()] = (np.asarray(t, float) / 60.0, np.asarray(y, float))  # seconds -> min
    return out


def on_grid(t, y, grid):
    order = np.argsort(t)
    return np.interp(grid, t[order], y[order], left=0.0, right=0.0)


def lag_minutes(t_ref, y_ref, t_mov, y_mov, step=0.005, max_lag=1.0):
    """Delay of the MS trace after the detector trace (positive = MS later)."""
    lo, hi = max(t_ref.min(), t_mov.min()), min(t_ref.max(), t_mov.max())
    grid = np.arange(lo, hi, step)
    a = on_grid(t_ref, y_ref, grid)
    b = on_grid(t_mov, y_mov, grid)
    n = int(max_lag / step)
    best, best_r = 0, -np.inf
    for k in range(-n, n + 1):
        x = a[max(0, -k):len(a) - max(0, k)]
        z = b[max(0, k):len(b) - max(0, -k)]
        r = np.corrcoef(x, z)[0, 1] if x.std() and z.std() else -np.inf
        if r > best_r:
            best, best_r = k, r
    return best * step, best_r


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mzml", nargs="+", type=Path)
    p.add_argument("--channel", default="Emission", help="text that identifies the detector trace ID")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--min-prominence", type=float, default=0.02, help="peak prominence as a fraction of the trace maximum")
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    rows = []
    for f in a.mzml:
        print(f"########## {f.name}")
        chroms = load_chromatograms(f)
        for cid, (t, y) in chroms.items():
            print(f"  {cid}: {len(t)} points, {t.min():.2f} to {t.max():.2f} min, max {y.max():.4g}")
        det = [k for k in chroms if a.channel.lower() in k.lower()]
        bpc = [k for k in chroms if "basepeak" in k.lower()]
        if not det:
            print(f"  no trace containing '{a.channel}'")
            rows.append({"file": f.name, "detector": "none"})
            continue
        if not bpc:
            print("  no base-peak chromatogram; cannot align")
            continue
        td, yd = chroms[det[0]]
        tb, yb = chroms[bpc[0]]
        peaks, _ = find_peaks(yd, prominence=a.min_prominence * yd.max())
        lag, r = lag_minutes(td, yd, tb, yb)
        print(f"  detector '{det[0]}': {len(peaks)} peaks (prominence >= {a.min_prominence:.0%} of max); "
              f"MS minus detector delay {lag * 60:+.1f} s, correlation {r:.3f}")
        rows.append({"file": f.name, "detector": det[0], "n_peaks": len(peaks),
                     "delay_s": round(lag * 60, 1), "correlation": round(r, 3)})

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 4.5), sharex=True, constrained_layout=True)
        ax1.plot(td, yd, lw=0.8, color="#1f4e79")
        ax1.plot(td[peaks], yd[peaks], "v", ms=3, color="#c00000")
        ax1.set_ylabel(det[0].replace("_0", ""))
        ax2.plot(tb - lag, yb, lw=0.8, color="#404040")
        ax2.set_ylabel("MS base peak")
        ax2.set_xlabel(f"Retention time (min); MS shifted by {lag * 60:+.1f} s")
        fig.suptitle(f.stem, fontsize=9)
        png = a.out / f"{f.stem}_detector_vs_bpc.png"
        fig.savefig(png, dpi=200)
        plt.close(fig)
        print(f"  figure: {png}")

    summary = a.out / "detector_summary.csv"
    keys = ["file", "detector", "n_peaks", "delay_s", "correlation"]
    with summary.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"summary: {summary}")


if __name__ == "__main__":
    main()
