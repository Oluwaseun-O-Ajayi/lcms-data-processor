"""Case study 4: resolution and peak shape against theory.

Question: does the reported resolution match the theoretical value for two
Gaussian peaks, and does the asymmetry factor rise with tailing?
Input:  data/case_studies/case4_separation_and_shape.csv
        (peak pairs with sigma 0.1 min at four separations, theoretical
        Rs = separation / (4 sigma); one peak at three tailing levels)
Output: examples/case_study_results/case4_resolution.csv
        examples/case_study_results/case4_peak_shape.csv
"""

import pandas as pd
from _case_study_common import DATA, RESULTS, TRUTH, analyze
from lcms_data_processor import QualityControl

data = pd.read_csv(DATA / "case4_separation_and_shape.csv")
truth = TRUTH[TRUTH["case"] == 4].set_index("scenario")
qc = QualityControl()

res_rows, shape_rows = [], []
for scenario, trace in data.groupby("scenario", sort=False):
    det = analyze(trace)
    if scenario.startswith("separation"):
        r = qc.check_resolution(det.peaks)
        res_rows.append({"separation_min": float(scenario.split("_")[1]),
                         "true_resolution": truth.loc[scenario, "true_resolution"],
                         "peaks_detected": len(det.peaks),
                         "reported_resolution": round(r["resolution"].iloc[0], 2) if len(r) else None,
                         "quality": r["quality"].iloc[0] if len(r) else None})
    else:
        s = qc.check_peak_shape(det.peaks, det.data, det.baseline)
        shape_rows.append({"tailing_K": truth.loc[scenario, "tailing_K"],
                           "asymmetry_factor": round(s["asymmetry_factor"].iloc[0], 2),
                           "quality": s["quality"].iloc[0]})

res, shape = pd.DataFrame(res_rows), pd.DataFrame(shape_rows)
res.to_csv(RESULTS / "case4_resolution.csv", index=False)
shape.to_csv(RESULTS / "case4_peak_shape.csv", index=False)
print(res.to_string(index=False))
print()
print(shape.to_string(index=False))
