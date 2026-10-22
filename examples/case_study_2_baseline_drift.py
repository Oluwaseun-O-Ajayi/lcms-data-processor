"""Case study 2: the same peak on flat, drifting and curved baselines.

Question: does baseline correction recover the true area when the baseline is
not flat?
Input:  data/case_studies/case2_baselines.csv (true area 3000, five replicates each)
Output: examples/case_study_results/case2_baseline_drift.csv
"""

import pandas as pd
from _case_study_common import DATA, RESULTS, analyze, nearest

data = pd.read_csv(DATA / "case2_baselines.csv")
TRUE_AREA, RT = 3000.0, 6.0

rows = []
for (baseline, rep), trace in data.groupby(["baseline", "replicate"], sort=False):
    p = nearest(analyze(trace).peaks, RT)
    rows.append({"baseline": baseline, "replicate": rep, "recovery_pct": round(100 * p["peak_area"] / TRUE_AREA, 1)})

table = pd.DataFrame(rows)
table.to_csv(RESULTS / "case2_baseline_drift.csv", index=False)
summary = table.groupby("baseline", sort=False)["recovery_pct"].agg(["mean", "std", "min", "max"]).round(1)
print(summary.to_string())
