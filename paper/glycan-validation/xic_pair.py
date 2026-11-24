"""MS extracted-ion chromatogram and fluorescence trace of an isomer pair, with the retention window
boundary marked, to see whether the larger isomer's tail falls inside the smaller isomer's window.

Prints, for the MS XIC and for the fluorescence trace, the area on each side of --split (min) and the
ratio second/first, then saves an overlay figure.

Usage:
    python paper/glycan-validation/xic_pair.py --mzml run.mzML --mz 1936.7791 968.8932 646.2646 --z 1 2 3 \
        --start 51.5 --end 54.5 --split 52.96 --out afli_G1F_pair.png
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyopenms as oms
from scipy import integrate


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mzml", type=Path, required=True)
    p.add_argument("--mz", type=float, nargs="+", required=True, help="monoisotopic m/z of each charge state")
    p.add_argument("--z", type=int, nargs="+", required=True, help="charge of each --mz value")
    p.add_argument("--start", type=float, required=True)
    p.add_argument("--end", type=float, required=True)
    p.add_argument("--split", type=float, required=True, help="window boundary between the two isomers (min)")
    p.add_argument("--channel", default="Emission")
    p.add_argument("--ppm", type=float, default=10.0)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()

    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(a.mzml), exp)
    if len(a.z) != len(a.mz):
        raise SystemExit("--z needs one charge per --mz value")
    targets = np.array([m + k * 1.003355 / z for m, z in zip(a.mz, a.z) for k in range(3)])   # first three isotopes
    rt, xic = [], []
    for s in exp:
        if s.getMSLevel() != 1:
            continue
        t = s.getRT() / 60.0
        if not (a.start <= t <= a.end):
            continue
        mz, it = s.get_peaks()
        hit = (np.abs(mz[:, None] - targets[None, :]) / targets[None, :] * 1e6 <= a.ppm).any(axis=1) if len(mz) else []
        rt.append(t)
        xic.append(float(it[hit].sum()) if len(mz) else 0.0)
    rt, xic = np.array(rt), np.array(xic)
    fld = [c for c in exp.getChromatograms() if a.channel.lower() in c.getNativeID().lower()][0]
    ft, fy = fld.get_peaks()
    ft = np.asarray(ft) / 60.0
    sel = (ft >= a.start) & (ft <= a.end)
    ft, fy = ft[sel], np.asarray(fy)[sel]
    fy = fy - np.percentile(fy, 5)

    def split_areas(t, y):
        left, right = t <= a.split, t > a.split
        return integrate.trapezoid(y[left], t[left]), integrate.trapezoid(y[right], t[right])

    for name, (t, y) in (("MS XIC", (rt, xic)), ("FLD", (ft, fy))):
        l, r = split_areas(t, y)
        print(f"{name:7s} area before {a.split} min: {l:.4g}   after: {r:.4g}   after/before: {r / l:.3f}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4), sharex=True, constrained_layout=True)
    ax1.plot(ft, fy, lw=0.8, color="#1f4e79"); ax1.set_ylabel("Fluorescence")
    ax2.plot(rt, xic, lw=0.8, color="#404040"); ax2.set_ylabel("MS XIC")
    for ax in (ax1, ax2):
        ax.axvline(a.split, color="#c00000", lw=0.8, ls="--")
        ax.set_yscale("log")
    ax2.set_xlabel("Retention time (min); dashed line = window boundary")
    fig.suptitle(a.mzml.stem, fontsize=8)
    fig.savefig(a.out, dpi=200)
    print(f"figure: {a.out}")


if __name__ == "__main__":
    main()
