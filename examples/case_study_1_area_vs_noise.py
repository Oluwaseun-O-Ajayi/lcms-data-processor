"""Case study 1: peak areas at three noise levels.

Question: how accurately are known peak areas recovered as noise increases, and
how does the error relate to signal-to-noise?
Input:  data/case_studies/case1_noise_levels.csv (four peaks, noise SD 20, 100, 300)
Output: examples/case_study_results/case1_area_vs_noise.csv
"""

import pandas as pd
from _case_study_common import DATA, RESULTS, TRUTH, analyze, nearest

data = pd.read_csv(DATA / "case1_noise_levels.csv")
truth = TRUTH[TRUTH["case"] == 1]

rows = []
for noise, trace in data.groupby("noise_sd"):
    peaks = analyze(trace).peaks
    for _, t in truth[truth["scenario"] == f"noise_sd={noise}"].iterrows():
        p = nearest(peaks, t["retention_time"])
        rows.append({"noise_sd": noise, "retention_time": t["retention_time"], "true_area": t["area"],
                     "measured_area": None if p is None else round(p["peak_area"], 1),
                     "recovery_pct": None if p is None else round(100 * p["peak_area"] / t["area"], 1),
                     "snr": None if p is None else round(p["snr"], 1)})

table = pd.DataFrame(rows)
table.to_csv(RESULTS / "case1_area_vs_noise.csv", index=False)
print(table.to_string(index=False))
