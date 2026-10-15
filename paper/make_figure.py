"""Regenerate paper/figure.png.

(a) The bundled example chromatogram (data/sample1.csv) after processing: baseline,
detected apexes and integration windows.
(b) Recovered peak area on simulated chromatograms with a known Gaussian peak on a
drifting baseline, as a function of signal-to-noise ratio (30 replicates per level).
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lcms_data_processor import ChromatogramParser, PeakDetector  # noqa: E402

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.6), dpi=300, gridspec_kw={"width_ratios": [1.5, 1]})

data = ChromatogramParser(ROOT / "data" / "sample1.csv").get_chromatogram()
det = PeakDetector(data)
det.detect_peaks()
det.integrate_peaks()
t, y = data["time"].to_numpy(), data["intensity"].to_numpy()
ax1.plot(t, y, color="0.3", lw=0.8, label="Chromatogram")
ax1.plot(t, det.baseline, color="#c0504d", lw=1.2, label="Estimated baseline")
for p in det.peaks.itertuples():
    sl = slice(int(p.left_base), int(p.right_base) + 1)
    ax1.fill_between(t[sl], np.asarray(det.baseline)[sl], y[sl], color="#1f6fb2", alpha=0.3, lw=0)
    ax1.annotate(str(int(p.peak_id)), (p.retention_time, y[int(p.apex_index)]), xytext=(0, 4),
                 textcoords="offset points", ha="center", fontsize=8)
ax1.set_xlabel("Retention time (min)")
ax1.set_ylabel("Intensity")
ax1.set_title("(a) Processed example chromatogram", fontsize=9, loc="left")
ax1.legend(fontsize=7, frameon=False)

rng = np.random.default_rng(0)
T = np.linspace(0, 20, 4000)
area, sigma, noise = 1000.0, 0.2, 30.0
levels = [5, 10, 25, 50, 100, 400]
means, sds = [], []
for snr in levels:
    height = snr * noise
    a_true = height * sigma * np.sqrt(2 * np.pi)
    ratios = []
    for _ in range(30):
        yy = 1000 + 2 * T + height * np.exp(-0.5 * ((T - 8) / sigma) ** 2) + rng.normal(0, noise, T.size)
        d = PeakDetector(pd.DataFrame({"time": T, "intensity": yy}))
        d.detect_peaks()
        d.integrate_peaks()
        pk = d.peaks.iloc[(d.peaks["retention_time"] - 8).abs().argmin()]
        ratios.append(pk["peak_area"] / a_true)
    means.append(np.mean(ratios))
    sds.append(np.std(ratios))
ax2.errorbar(levels, means, yerr=sds, fmt="o-", color="#1f6fb2", capsize=3, ms=4)
ax2.axhline(1.0, color="k", ls=":", lw=1)
ax2.set_xscale("log")
ax2.set_xlabel("Signal-to-noise ratio (height / noise SD)")
ax2.set_ylabel("Recovered area / true area")
ax2.set_title("(b) Area recovery on simulated peaks", fontsize=9, loc="left")
for a in (ax1, ax2):
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(ROOT / "paper" / "figure.png")
print("S/N levels", levels)
print("mean ratio", np.round(means, 3).tolist())
print("SD", np.round(sds, 3).tolist())
