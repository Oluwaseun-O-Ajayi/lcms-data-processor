"""Diagnostics for each target retention window: how much of the MS1 signal belongs to the target glycan,
and which other ions dominate the window.

For every target glycan (all charge states in the targets file, first three isotopes, +/- --ppm) and every
mzML file, sums over the MS1 scans inside the glycan's retention window:
    tic             total ion current
    glycan          intensity of the glycan's own ions
    glycan_frac     glycan / tic
and lists the five most intense other ions in the window (m/z rounded to 0.01).

A low glycan fraction with intense foreign ions points to co-elution, which can suppress or distort the
glycan's MS signal; a foreign ion that is itself a glycan (for example a sialylated form that loses its
sialic acid in the source) can add signal at the target m/z.

Usage:
    python paper/glycan-validation/diag_windows.py --targets cluster/ms119_cetux_RFMS_targets.csv \
        --mzml /scratch/$USER/glycan_data/pride/mzml_all/MS119_Joined_Cetux*.mzML --out diag_cetux.csv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pyopenms as oms

ISO = 1.003355


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--targets", type=Path, required=True)
    p.add_argument("--mzml", type=Path, nargs="+", required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--ppm", type=float, default=10.0)
    a = p.parse_args()

    t = pd.read_csv(a.targets)
    windows = t.groupby("glycan").agg(rt_min=("rt_min", "first"), rt_max=("rt_max", "first")).reset_index()
    ions = {g: np.array([m + k * ISO / z for m, z in zip(d["mz"], d["charge"]) for k in range(3)])
            for g, d in t.groupby("glycan")}

    rows = []
    for f in a.mzml:
        exp = oms.MSExperiment()
        oms.MzMLFile().load(str(f), exp)
        ms1 = [s for s in exp if s.getMSLevel() == 1]
        rts = np.array([s.getRT() / 60.0 for s in ms1])
        for _, w in windows.iterrows():
            target = ions[w["glycan"]]
            tic = glyc = 0.0
            other = {}
            for i in np.where((rts >= w["rt_min"]) & (rts <= w["rt_max"]))[0]:
                mz, it = ms1[i].get_peaks()
                if len(mz) == 0:
                    continue
                tic += it.sum()
                d = np.abs(mz[:, None] - target[None, :]) / target[None, :] * 1e6
                hit = (d <= a.ppm).any(axis=1)
                glyc += it[hit].sum()
                for m, v in zip(np.round(mz[~hit], 2), it[~hit]):
                    other[m] = other.get(m, 0.0) + v
            top = sorted(other.items(), key=lambda kv: -kv[1])[:5]
            rows.append({"file": f.name, "glycan": w["glycan"], "rt_min": w["rt_min"], "rt_max": w["rt_max"],
                         "tic": tic, "glycan_intensity": glyc, "glycan_frac": glyc / tic if tic else np.nan,
                         "top_other_ions": "; ".join(f"{m:.2f} ({v / tic:.1%})" for m, v in top) if tic else ""})
        print(f"{f.name}: {len(windows)} windows", flush=True)

    out = pd.DataFrame(rows)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    summary = out.groupby("glycan")["glycan_frac"].agg(["mean", "min", "max"]).sort_values("mean")
    print("\nGlycan share of the MS1 signal in its own window (mean, min, max over files):")
    print(summary.to_string(float_format=lambda v: f"{v:.3f}"))
    print(f"\nWrote {a.out}")


if __name__ == "__main__":
    main()
