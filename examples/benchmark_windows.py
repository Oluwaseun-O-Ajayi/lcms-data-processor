"""
Synthetic benchmark: how the extraction window trades off lost signal against interference.

Simulated LC-MS runs with a known analyte peak are written as mzML and processed with
MzMLReader.xic_matrix and PeakDetector, the same code path used for real data.

Scenario A (mass error): every centroid of the analyte is shifted by a systematic error
(0, 2 or 4 ppm) plus random scatter (SD 1 ppm per centroid). Windows narrower than the
error lose signal.

Scenario B (interference): an unrelated ion with the same elution profile and half the
analyte's intensity sits 10, 35 or 70 ppm away. Windows wider than that offset include it.

Usage:
    python examples/benchmark_windows.py --out results/benchmark
"""

import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lcms_data_processor import MzMLReader, PeakDetector  # noqa: E402

MZ = 1000.7004          # A3G3S3 [M+3H]3+ with 2-AB, as a realistic m/z
RT, SIGMA, AREA = 10.0, 0.08, 1.0e5
WINDOWS = [1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 100]
SCATTER_PPM = 1.0
SEED = 11


def simulate(path, shift_ppm=0.0, interferent_ppm=None, interferent_frac=0.5, seed=SEED):
    import pyopenms as oms
    rng = np.random.default_rng(seed)
    exp = oms.MSExperiment()
    for rt_s in np.arange(0, 20 * 60, 0.5):
        t = rt_s / 60
        prof = AREA / (SIGMA * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((t - RT) / SIGMA) ** 2)
        mz = list(rng.uniform(400, 2000, 40))
        inten = list(rng.exponential(200, 40))
        if prof > 1:
            err = shift_ppm + rng.normal(0, SCATTER_PPM)
            mz.append(MZ * (1 + err * 1e-6)); inten.append(prof)
            if interferent_ppm is not None:
                mz.append(MZ * (1 + (interferent_ppm + rng.normal(0, SCATTER_PPM)) * 1e-6))
                inten.append(interferent_frac * prof)
        o = np.argsort(mz)
        s = oms.MSSpectrum(); s.setRT(float(rt_s)); s.setMSLevel(1)
        s.set_peaks((np.asarray(mz)[o], np.asarray(inten)[o]))
        exp.addSpectrum(s)
    oms.MzMLFile().store(str(path), exp)


def measure(path):
    reader = MzMLReader(path)
    rows = []
    for w in WINDOWS:
        times, mat = reader.xic_matrix([MZ], ppm=w)
        det = PeakDetector(pd.DataFrame({'time': times, 'intensity': mat[:, 0]}))
        peaks = det.detect_peaks(width=5, baseline_window=0.1, threshold_range=(RT - 0.6, RT + 0.6))
        area = 0.0
        if not peaks.empty:
            det.integrate_peaks()
            near = det.peaks[(det.peaks['retention_time'] - RT).abs() < 0.3]
            area = near['peak_area'].sum()
        rows.append({'window_ppm': w, 'area': area, 'recovery_pct': 100 * area / AREA})
    return pd.DataFrame(rows)


def run(out):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    frames = []
    with tempfile.TemporaryDirectory() as d:
        for shift in (0.0, 2.0, 4.0):
            p = Path(d) / f"shift{shift}.mzML"
            simulate(p, shift_ppm=shift)
            frames.append(measure(p).assign(scenario='mass error', condition=f'{shift:g} ppm shift'))
        for off in (10.0, 35.0, 70.0):
            p = Path(d) / f"interf{off}.mzML"
            simulate(p, interferent_ppm=off)
            frames.append(measure(p).assign(scenario='interference', condition=f'ion at {off:g} ppm'))
    res = pd.concat(frames, ignore_index=True)
    res.to_csv(f"{out}.csv", index=False)
    return res


def plot(res, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                         "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "savefig.dpi": 600, "pdf.fonttype": 42})
    colors = ["#2a78d6", "#eb6834", "#1baf7a"]
    marks = ["o", "s", "^"]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), sharey=True)
    for ax, scen, title in [(axes[0], 'mass error', 'Mass error (1 ppm scatter)'),
                            (axes[1], 'interference', 'Co-eluting ion, 50% intensity')]:
        sub = res[res.scenario == scen]
        for k, (cond, g) in enumerate(sub.groupby('condition', sort=False)):
            ax.plot(g.window_ppm, g.recovery_pct, color=colors[k], marker=marks[k], ms=4, mec="white",
                    mew=0.6, lw=1.4, label=cond)
        ax.axhline(100, color="#52514e", lw=0.8, ls="--", zorder=0)
        ax.set_xscale("log"); ax.set_xticks(WINDOWS); ax.set_xticklabels([str(w) for w in WINDOWS]); ax.minorticks_off()
        ax.set_xlabel("Extraction window (± ppm)"); ax.set_title(title, fontsize=8)
        ax.yaxis.grid(True, color="#d9d8d4", lw=0.5); ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=7, loc="lower right" if scen == 'mass error' else "upper left")
    axes[0].set_ylabel("Recovered area (% of true)")
    axes[0].set_ylim(0, 165)
    for ax, letter in zip(axes, "AB"):
        ax.text(-0.14, 1.04, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
    fig.tight_layout(w_pad=2)
    fig.savefig(f"{out}.png"); fig.savefig(f"{out}.pdf")


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--out', default='results/benchmark_windows')
    a = p.parse_args()
    res = run(a.out)
    plot(res, a.out)
    print(res.pivot_table(index='window_ppm', columns='condition', values='recovery_pct').round(1).to_string())
