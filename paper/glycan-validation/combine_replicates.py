"""Combine replicate runs: MS-vs-fluorescence agreement, repeatability and settings sensitivity.

Reads, for every run R in --results:
    R_fld_per_glycan.csv   (lcms-compare fld)   MS area and FLD area per glycan
    R_glycans.csv          (lcms-batch)         MS area per glycan for every processing setting
Run names must contain a day tag and a replicate tag, matched by --day-regex and --rep-regex.

Writes to --out:
    per_run_per_glycan.csv      MS % and FLD % per glycan per run, on the common glycan set
    agreement_per_run.csv       MAD, limits of agreement, r per run (default setting)
    agreement_per_glycan.csv    mean MS %, mean FLD %, bias and 95% limits per glycan (all runs)
    repeatability.csv           within-day and between-day CV per glycan, MS and FLD separately
    settings_sensitivity.csv    MAD vs FLD for every processing setting, per run
    calibration_cross_day.csv   response factors fitted on one day, tested on the other days
    Figure_replicates.png/.pdf  (a) MS % vs FLD % per glycan, mean +/- SD; (b) Bland-Altman, all runs

Percentages are recomputed over the glycans matched to a fluorescence peak in every run, so
every run is normalised over the same set and runs are comparable.

Usage:
    python paper/glycan-validation/combine_replicates.py --results /scratch/$USER/glycan_data/pride/results/trast \
        --out paper/glycan-validation/data/trastuzumab_RFMS
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sst
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def stats(ms, fld):
    d = np.asarray(ms) - np.asarray(fld)
    sd = d.std(ddof=1) if len(d) > 1 else np.nan
    return {"n": len(d), "mean_abs_diff": np.abs(d).mean(), "bias": d.mean(),
            "loa_low": d.mean() - 1.96 * sd, "loa_high": d.mean() + 1.96 * sd,
            "pearson_r": np.corrcoef(ms, fld)[0, 1] if len(d) > 2 else np.nan}


def cv(x):
    x = np.asarray(x, float)
    return 100 * x.std(ddof=1) / x.mean() if len(x) > 1 and x.mean() else np.nan


def variance_components(df, value):
    """One-way ANOVA (day as group): within-day (repeatability) and between-day CV, in % of the mean."""
    groups = [g[value].to_numpy() for _, g in df.groupby("day")]
    k = len(groups)
    n = np.array([len(g) for g in groups])
    grand = df[value].mean()
    if k < 2 or (n < 2).any():
        return np.nan, np.nan, np.nan
    ss_w = sum(((g - g.mean()) ** 2).sum() for g in groups)
    ms_w = ss_w / (n.sum() - k)
    ss_b = sum(len(g) * (g.mean() - grand) ** 2 for g in groups)
    ms_b = ss_b / (k - 1)
    n0 = (n.sum() - (n ** 2).sum() / n.sum()) / (k - 1)
    var_b = max((ms_b - ms_w) / n0, 0.0)
    within = 100 * np.sqrt(ms_w) / grand
    between = 100 * np.sqrt(var_b) / grand
    intermediate = 100 * np.sqrt(ms_w + var_b) / grand
    return within, between, intermediate


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--day-regex", default=r"(\d+)(?:st|nd|rd|th)Day")
    p.add_argument("--rep-regex", default=r"_(\d+)$")
    p.add_argument("--ppm", type=float, default=10.0)
    p.add_argument("--width", type=int, default=5)
    p.add_argument("--baseline-window", type=float, default=0.1)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    files = sorted(a.results.glob("*_fld_per_glycan.csv"))
    if not files:
        raise SystemExit(f"No *_fld_per_glycan.csv files in {a.results}")
    rows = []
    for f in files:
        run = f.name[: -len("_fld_per_glycan.csv")]
        day = re.search(a.day_regex, run)
        rep = re.search(a.rep_regex, run)
        if not (day and rep):
            raise SystemExit(f"Cannot read day/replicate from run name '{run}'; adjust --day-regex/--rep-regex")
        t = pd.read_csv(f)
        t["run"], t["day"], t["rep"] = run, int(day.group(1)), int(rep.group(1))
        rows.append(t)
    allg = pd.concat(rows, ignore_index=True)
    runs = allg["run"].unique()
    print(f"{len(runs)} runs: " + ", ".join(runs))

    status = allg.pivot_table(index="glycan", columns="run", values="status", aggfunc="first")
    print("\nFluorescence match status per glycan (runs matched / total):")
    for g, srow in status.iterrows():
        print(f"  {g:8s} {int((srow == 'matched').sum())}/{len(runs)}")
    common = [g for g, srow in status.iterrows() if (srow == "matched").all()]
    print(f"Common glycan set (matched in every run): {', '.join(common)}")

    c = allg[allg["glycan"].isin(common)].copy()
    c["ms_pct"] = 100 * c["ms_area"] / c.groupby("run")["ms_area"].transform("sum")
    c["fld_pct"] = 100 * c["fld_area"] / c.groupby("run")["fld_area"].transform("sum")
    c["diff"] = c["ms_pct"] - c["fld_pct"]
    c[["run", "day", "rep", "glycan", "ms_rt", "fld_rts", "ms_area", "fld_area", "ms_pct", "fld_pct", "diff"]] \
        .sort_values(["day", "rep", "glycan"]).to_csv(a.out / "per_run_per_glycan.csv", index=False)

    def coverage(run, glycans):
        """Share of the whole baseline-corrected fluorescence area of the glycan region (detected or not)
        taken by the given glycans; read from lcms-compare's summary (fld_region_area)."""
        f = a.results / f"{run}_fld_summary.csv"
        if not f.exists():
            return np.nan
        sm = pd.read_csv(f)
        if "fld_region_area" not in sm:
            return np.nan
        total = float(sm["fld_region_area"].iat[0])
        sel = allg[(allg["run"] == run) & (allg["glycan"].isin(glycans))]["fld_area"].sum()
        return 100 * sel / total if total else np.nan

    per_run = pd.DataFrame([{"run": r, "day": g["day"].iat[0], "rep": g["rep"].iat[0],
                             **stats(g["ms_pct"], g["fld_pct"]),
                             "fld_coverage_pct": coverage(r, common)} for r, g in c.groupby("run")])
    per_run = per_run.sort_values(["day", "rep"])
    per_run.to_csv(a.out / "agreement_per_run.csv", index=False)

    per_gly = []
    for g, d in c.groupby("glycan"):
        s = stats(d["ms_pct"], d["fld_pct"])
        half = sst.t.ppf(0.975, len(d) - 1) * d["diff"].std(ddof=1) / np.sqrt(len(d))
        per_gly.append({"glycan": g, "n_runs": len(d), "ms_pct_mean": d["ms_pct"].mean(), "ms_pct_sd": d["ms_pct"].std(ddof=1),
                        "fld_pct_mean": d["fld_pct"].mean(), "fld_pct_sd": d["fld_pct"].std(ddof=1),
                        "bias": s["bias"], "bias_ci95_low": s["bias"] - half, "bias_ci95_high": s["bias"] + half,
                        "loa_low": s["loa_low"], "loa_high": s["loa_high"]})
    per_gly = pd.DataFrame(per_gly).sort_values("fld_pct_mean", ascending=False)
    per_gly.to_csv(a.out / "agreement_per_glycan.csv", index=False)

    rep_rows = []
    for g, d in c.groupby("glycan"):
        r = {"glycan": g}
        for col, lab in (("ms_pct", "ms"), ("fld_pct", "fld")):
            w, b, i = variance_components(d, col)
            r.update({f"{lab}_mean": d[col].mean(), f"{lab}_cv_within_day": w,
                      f"{lab}_cv_between_day": b, f"{lab}_cv_intermediate": i})
        rep_rows.append(r)
    rep = pd.DataFrame(rep_rows).sort_values("fld_mean", ascending=False)
    rep.to_csv(a.out / "repeatability.csv", index=False)

    pooled = stats(c["ms_pct"], c["fld_pct"])
    cov = per_run["fld_coverage_pct"]
    print(f"\nFluorescence coverage of the common glycan set (share of the whole baseline-corrected FLD area of the glycan region): "
          f"mean {cov.mean():.1f}%, range {cov.min():.1f} to {cov.max():.1f}%")
    print("\nPooled over all runs (default setting): " +
          ", ".join(f"{k}={v:.4g}" for k, v in pooled.items()))
    print("\nPer run:\n" + per_run.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("\nPer glycan:\n" + per_gly.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    print("\nRepeatability (CV %, of the mean relative abundance):\n" +
          rep.to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    # cross-day calibration: per-glycan response factor RF = mean(MS % / FLD %) on the training day;
    # on every other day the MS areas are divided by RF, renormalised, and compared with FLD again.
    cal_rows = []
    for train_day in sorted(c["day"].unique()):
        tr = c[c["day"] == train_day]
        rf = (tr["ms_pct"] / tr["fld_pct"]).groupby(tr["glycan"]).mean()
        te = c[c["day"] != train_day].copy()
        te["ms_area_cal"] = te["ms_area"] / te["glycan"].map(rf)
        te["ms_pct_cal"] = 100 * te["ms_area_cal"] / te.groupby("run")["ms_area_cal"].transform("sum")
        for run, d in te.groupby("run"):
            before = stats(d["ms_pct"], d["fld_pct"])
            after = stats(d["ms_pct_cal"], d["fld_pct"])
            cal_rows.append({"train_day": train_day, "test_run": run, "test_day": d["day"].iat[0],
                             "mad_before": before["mean_abs_diff"], "mad_after": after["mean_abs_diff"],
                             "loa_halfwidth_before": before["loa_high"] - before["bias"],
                             "loa_halfwidth_after": after["loa_high"] - after["bias"]})
        print(f"\nResponse factors fitted on day {train_day} (MS % / FLD %): " +
              ", ".join(f"{g} {v:.3f}" for g, v in rf.items()))
    cal = pd.DataFrame(cal_rows)
    cal.to_csv(a.out / "calibration_cross_day.csv", index=False)
    print("\nCross-day calibration, test runs only (points):\n" +
          cal.groupby("train_day")[["mad_before", "mad_after", "loa_halfwidth_before", "loa_halfwidth_after"]]
          .mean().to_string(float_format=lambda v: f"{v:.3f}"))
    print(f"All train/test pairs: MAD {cal['mad_before'].mean():.3f} -> {cal['mad_after'].mean():.3f} points")

    # settings sensitivity: MS % per setting (lcms-batch) vs the run's FLD % (setting-independent)
    sens = []
    fld_ref = c.set_index(["run", "glycan"])["fld_pct"]
    for f in sorted(a.results.glob("*_glycans.csv")):
        run = f.name[: -len("_glycans.csv")]
        if run not in runs:
            continue
        b = pd.read_csv(f)
        b = b[b["glycan"].isin(common)]
        for key, d in b.groupby(["ppm", "width", "baseline_window"]):
            pct = 100 * d["peak_area"] / d["peak_area"].sum()
            ref = fld_ref.loc[run].reindex(d["glycan"]).to_numpy()
            s = stats(pct.to_numpy(), ref)
            sens.append({"run": run, "ppm": key[0], "width": key[1], "baseline_window": key[2],
                         "mean_abs_diff": s["mean_abs_diff"], "pearson_r": s["pearson_r"]})
    if sens:
        sens = pd.DataFrame(sens)
        sens.to_csv(a.out / "settings_sensitivity.csv", index=False)
        by = sens.groupby(["ppm", "width", "baseline_window"])["mean_abs_diff"].agg(["mean", "min", "max"])
        print("\nSettings sensitivity, MAD vs FLD across runs (points):\n" + by.to_string(float_format=lambda v: f"{v:.3f}"))

    # figure
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.8,
                         "xtick.major.width": 0.8, "ytick.major.width": 0.8})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.3), constrained_layout=True)
    order = per_gly["glycan"].tolist()
    x = np.arange(len(order))
    pg = per_gly.set_index("glycan").loc[order]
    w = 0.38
    ax1.bar(x - w / 2, pg["fld_pct_mean"], w, yerr=pg["fld_pct_sd"], capsize=2, color="#9e9e9e",
            edgecolor="black", linewidth=0.6, label="Fluorescence", error_kw={"elinewidth": 0.6})
    ax1.bar(x + w / 2, pg["ms_pct_mean"], w, yerr=pg["ms_pct_sd"], capsize=2, color="#1f4e79",
            edgecolor="black", linewidth=0.6, label="MS (this work)", error_kw={"elinewidth": 0.6})
    ax1.set_xticks(x, [o.replace("_a", " (a)").replace("_b", " (b)") for o in order])
    ax1.set_ylabel("Relative abundance (%)")
    ax1.legend(frameon=False)
    ax1.set_title("a", loc="left", fontweight="bold")
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)

    markers = "osD^v<>ph"
    for k, g in enumerate(order):
        d = c[c["glycan"] == g]
        ax2.scatter((d["ms_pct"] + d["fld_pct"]) / 2, d["diff"], s=14, marker=markers[k % len(markers)],
                    facecolors="none", edgecolors="black", linewidths=0.7, label=g.replace("_a", " (a)").replace("_b", " (b)"))
    for yv, ls in ((pooled["bias"], "-"), (pooled["loa_low"], "--"), (pooled["loa_high"], "--")):
        ax2.axhline(yv, color="#1f4e79", lw=0.8, ls=ls)
    ax2.text(0.02, pooled["loa_high"], f"+1.96 SD ({pooled['loa_high']:.2f})", va="bottom", ha="left",
             transform=ax2.get_yaxis_transform(), fontsize=7)
    ax2.text(0.02, pooled["loa_low"], f"\u22121.96 SD ({pooled['loa_low']:.2f})", va="top", ha="left",
             transform=ax2.get_yaxis_transform(), fontsize=7)
    ax2.set_xlabel("Mean of MS and fluorescence (%)")
    ax2.set_ylabel("MS minus fluorescence (points)")
    ax2.legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=len(order),
               handletextpad=0.2, columnspacing=0.8)
    lim = max(abs(c["diff"]).max(), abs(pooled["loa_low"]), abs(pooled["loa_high"])) * 1.25
    ax2.set_ylim(-lim, lim)
    ax2.set_title("b", loc="left", fontweight="bold")
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    for ext in ("png", "pdf"):
        fig.savefig(a.out / f"Figure_replicates.{ext}", dpi=600)
    print(f"\nWrote tables and Figure_replicates.png/.pdf to {a.out}")


if __name__ == "__main__":
    main()
