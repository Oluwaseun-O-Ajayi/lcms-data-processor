"""Do response factors fitted on one sequence transfer to another sequence and to other proteins?

Reference: per_run_per_glycan.csv written by combine_replicates.py for the MS154 trastuzumab runs.
For each glycan g, RF_g = mean over the reference runs of (MS % / FLD %).

Test runs: every R_fld_per_glycan.csv in --results (lcms-compare fld output). The protein is read
from the run name (Trast, Afli, Cetux). Three analyses:

  1. full profile      MS vs FLD over every glycan matched to a fluorescence peak (no correction)
  2. transfer          glycans that have a reference RF; MS and FLD renormalised over that subset;
                       agreement before and after dividing MS areas by RF_g
  3. within protein    RFs fitted on one replicate of a protein, applied to its other replicates,
                       over the glycans matched in all replicates of that protein

Writes full_profile.csv, transfer.csv, within_protein.csv and a printed summary to --out.

Usage:
    python paper/glycan-validation/transfer_test.py \
        --ref paper/glycan-validation/data/trastuzumab_RFMS/per_run_per_glycan.csv \
        --results /scratch/$USER/glycan_data/pride/results/ms119 \
        --out paper/glycan-validation/data/ms119_transfer
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def stats(ms, fld):
    ms, fld = np.asarray(ms, float), np.asarray(fld, float)
    d = ms - fld
    sd = d.std(ddof=1) if len(d) > 1 else np.nan
    return {"n": len(d), "mean_abs_diff": np.abs(d).mean(), "loa_halfwidth": 1.96 * sd,
            "pearson_r": np.corrcoef(ms, fld)[0, 1] if len(d) > 2 else np.nan}


def pct(x):
    x = np.asarray(x, float)
    return 100 * x / x.sum()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", type=Path, required=True)
    p.add_argument("--results", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    ref = pd.read_csv(a.ref)
    rf = (ref["ms_pct"] / ref["fld_pct"]).groupby(ref["glycan"]).mean()
    print("Reference response factors (MS % / FLD %, mean of reference runs): " +
          ", ".join(f"{g} {v:.3f}" for g, v in rf.items()))

    rows = []
    for f in sorted(a.results.glob("*_fld_per_glycan.csv")):
        run = f.name[: -len("_fld_per_glycan.csv")]
        m = re.search(r"(Trast|Afli|Cetux)", run)
        t = pd.read_csv(f)
        t["run"], t["protein"] = run, (m.group(1) if m else "unknown")
        rows.append(t)
    if not rows:
        raise SystemExit(f"No *_fld_per_glycan.csv files in {a.results}")
    allg = pd.concat(rows, ignore_index=True)
    ok = allg[allg["status"] == "matched"].copy()

    print("\nMatched glycans per run (of targets):")
    for run, d in allg.groupby("run"):
        print(f"  {run}: {int((d['status'] == 'matched').sum())}/{len(d)}"
              + ("" if (d["status"] == "matched").all() else
                 "   not matched: " + ", ".join(f"{g} ({s})" for g, s in zip(d["glycan"], d["status"]) if s != "matched")))

    # 1. full profile
    full = []
    for run, d in ok.groupby("run"):
        full.append({"run": run, "protein": d["protein"].iat[0], **stats(pct(d["ms_area"]), pct(d["fld_area"]))})
    full = pd.DataFrame(full)
    full.to_csv(a.out / "full_profile.csv", index=False)

    # 2. transfer of reference RFs
    tr = []
    for run, d in ok.groupby("run"):
        d = d[d["glycan"].isin(rf.index)]
        if len(d) < 2:
            continue
        before = stats(pct(d["ms_area"]), pct(d["fld_area"]))
        after = stats(pct(d["ms_area"] / d["glycan"].map(rf)), pct(d["fld_area"]))
        tr.append({"run": run, "protein": d["protein"].iat[0], "glycans": ";".join(d["glycan"]), "n": len(d),
                   "mad_before": before["mean_abs_diff"], "mad_after": after["mean_abs_diff"],
                   "loa_halfwidth_before": before["loa_halfwidth"], "loa_halfwidth_after": after["loa_halfwidth"]})
    tr = pd.DataFrame(tr)
    tr.to_csv(a.out / "transfer.csv", index=False)

    # 3. within-protein calibration across replicates
    wp = []
    for prot, d in ok.groupby("protein"):
        runs = sorted(d["run"].unique())
        if len(runs) < 2:
            continue
        common = set.intersection(*(set(d[d["run"] == r]["glycan"]) for r in runs))
        d = d[d["glycan"].isin(common)].copy()
        d["ms_pct"] = d.groupby("run")["ms_area"].transform(lambda x: 100 * x / x.sum())
        d["fld_pct"] = d.groupby("run")["fld_area"].transform(lambda x: 100 * x / x.sum())
        for train in runs:
            tdat = d[d["run"] == train]
            rf_p = (tdat["ms_pct"] / tdat["fld_pct"]).groupby(tdat["glycan"]).mean()
            for test in runs:
                if test == train:
                    continue
                e = d[d["run"] == test]
                before = stats(e["ms_pct"], e["fld_pct"])
                after = stats(pct(e["ms_area"] / e["glycan"].map(rf_p)), e["fld_pct"])
                wp.append({"protein": prot, "train": train, "test": test, "n": len(e),
                           "mad_before": before["mean_abs_diff"], "mad_after": after["mean_abs_diff"],
                           "loa_halfwidth_before": before["loa_halfwidth"],
                           "loa_halfwidth_after": after["loa_halfwidth"]})
    wp = pd.DataFrame(wp)
    wp.to_csv(a.out / "within_protein.csv", index=False)

    fmt = lambda v: f"{v:.3f}"
    print("\n1. Full profile, no correction (points):\n" +
          full.groupby("protein")[["n", "mean_abs_diff", "loa_halfwidth", "pearson_r"]].mean().to_string(float_format=fmt))
    if len(tr):
        print("\n2. Reference RFs applied to shared glycans (points):\n" +
              tr.groupby("protein")[["n", "mad_before", "mad_after", "loa_halfwidth_before", "loa_halfwidth_after"]]
              .mean().to_string(float_format=fmt))
        print("   glycans used: " + "; ".join(f"{p}: {g}" for p, g in tr.groupby("protein")["glycans"].first().items()))
    if len(wp):
        print("\n3. Within-protein calibration, one replicate -> the others (points):\n" +
              wp.groupby("protein")[["n", "mad_before", "mad_after", "loa_halfwidth_before", "loa_halfwidth_after"]]
              .mean().to_string(float_format=fmt))
    print(f"\nWrote full_profile.csv, transfer.csv, within_protein.csv to {a.out}")


if __name__ == "__main__":
    main()
