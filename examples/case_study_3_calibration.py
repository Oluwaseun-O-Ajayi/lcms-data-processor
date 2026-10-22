"""Case study 3: calibration and quantification of two unknowns.

Question: how accurately are unknowns quantified from a five-point calibration,
and how do the standards behave near the detection limit?
Input:  data/case_studies/case3_calibration.csv (true area = 30 x concentration,
        standards 1 to 100 ng/mL, unknowns at 20 and 75 ng/mL)
Output: examples/case_study_results/case3_calibration.csv
"""

import pandas as pd
from _case_study_common import DATA, RESULTS, analyze, nearest
from lcms_data_processor import CalibrationCurve

data = pd.read_csv(DATA / "case3_calibration.csv")
RT = 5.2

rows = []
for (kind, conc), trace in data.groupby(["sample_type", "concentration"], sort=False):
    p = nearest(analyze(trace).peaks, RT)
    rows.append({"sample_type": kind, "true_conc": conc, "true_area": 30 * conc,
                 "measured_area": None if p is None else round(p["peak_area"], 1),
                 "snr": None if p is None else round(p["snr"], 1)})
table = pd.DataFrame(rows)

std = table[table["sample_type"] == "standard"].dropna()
curve = CalibrationCurve()
curve.fit_calibration(std["true_conc"], std["measured_area"])
table["back_calculated"] = curve.quantify_samples(table["measured_area"]).round(2)
table["accuracy_pct"] = (100 * table["back_calculated"] / table["true_conc"]).round(1)

table.to_csv(RESULTS / "case3_calibration.csv", index=False)
print(f"Calibration: area = {curve.slope:.3f} x conc + {curve.intercept:.2f}, R^2 = {curve.r_squared:.5f}")
print("(true relationship: area = 30 x conc)\n")
print(table.to_string(index=False))
