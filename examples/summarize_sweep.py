"""
Summarize a processing-settings sweep from lcms-batch.

Reads one or more <run>_glycans.csv files and reports, for every file and glycan, how much
the relative abundance moves across settings: the median, the minimum and maximum, and
the range in percentage points. It also writes a figure of relative abundance against the
extraction window (median over the other settings, with the min-max span as error bars).

Usage:
    python examples/summarize_sweep.py results/*_glycans.csv --out results/summary
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402


def summarize(frames):
    df = pd.concat(frames, ignore_index=True)
    g = df.groupby(["file", "glycan"])["relative_abundance_pct"]
    table = g.agg(median="median", min="min", max="max").reset_index()
    table["range_pct_points"] = table["max"] - table["min"]
    by_ppm = (df.groupby(["file", "glycan", "ppm"])["relative_abundance_pct"]
                .agg(median="median", min="min", max="max").reset_index())
    return table.sort_values(["file", "median"], ascending=[True, False]), by_ppm


def plot(by_ppm, path):
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                         "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8})
    files = sorted(by_ppm["file"].unique())
    fig, axes = plt.subplots(1, len(files), figsize=(3.4 * len(files), 2.8), squeeze=False)
    for ax, f in zip(axes[0], files):
        sub = by_ppm[by_ppm["file"] == f]
        order = sub.groupby("glycan")["median"].max().sort_values(ascending=False).index
        for gname in order:
            s = sub[sub["glycan"] == gname].sort_values("ppm")
            ax.errorbar(s["ppm"], s["median"], yerr=[s["median"] - s["min"], s["max"] - s["median"]],
                        marker="o", ms=3.5, lw=1, capsize=2, label=gname)
        ax.set_xscale("log")
        ticks = sorted(sub["ppm"].unique())
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t:g}" for t in ticks])
        ax.minorticks_off()
        ax.set_xlabel("Extraction window (± ppm)")
        ax.set_ylabel("Relative abundance (%)")
        ax.set_title(Path(f).stem, fontsize=9)
        ax.legend(frameon=False, fontsize=7, ncol=1, loc="best")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=600)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("glycan_tables", nargs="+")
    p.add_argument("--out", default="results/summary")
    a = p.parse_args(argv)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    table, by_ppm = summarize([pd.read_csv(f) for f in a.glycan_tables])
    table.to_csv(f"{out}_spread.csv", index=False)
    by_ppm.to_csv(f"{out}_by_ppm.csv", index=False)
    plot(by_ppm, Path(f"{out}_abundance_vs_ppm"))
    print(table.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    print(f"\nWrote {out}_spread.csv, {out}_by_ppm.csv and {out}_abundance_vs_ppm.png/.pdf")


if __name__ == "__main__":
    main()
