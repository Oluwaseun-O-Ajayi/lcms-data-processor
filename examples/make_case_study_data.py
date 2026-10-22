"""Create the bundled case-study chromatograms (fixed seed, known answers).

Every trace is baseline + Gaussian (or exponentially modified Gaussian) peaks +
normal noise on a 0 to 15 min grid of 1501 points. Peak areas are specified
directly, so the true area of every peak is known. The true values are written
to data/case_studies/truth.csv.

Run: python examples/make_case_study_data.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import exponnorm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "case_studies"
OUT.mkdir(parents=True, exist_ok=True)
T = np.linspace(0, 15, 1501)
rng = np.random.default_rng(42)
truth = []


def gaussian(area, center, sigma):
    return area / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((T - center) / sigma) ** 2)


def trace(y, **labels):
    return pd.DataFrame({**labels, "time": np.round(T, 3), "intensity": np.round(y, 2)})


# Case 1: four peaks at three noise levels
PEAKS = [(3.5, 800, 0.15), (5.2, 3000, 0.20), (7.8, 2000, 0.25), (10.5, 500, 0.18)]
frames = []
for noise in (20, 100, 300):
    y = 1000 + sum(gaussian(a, c, s) for c, a, s in PEAKS) + rng.normal(0, noise, T.size)
    frames.append(trace(y, noise_sd=noise))
    truth += [{"case": 1, "scenario": f"noise_sd={noise}", "retention_time": c, "area": a, "sigma": s}
              for c, a, s in PEAKS]
pd.concat(frames).to_csv(OUT / "case1_noise_levels.csv", index=False)

# Case 2: one peak on three baselines, five replicates each
BASELINES = {"flat": 1000 + 0 * T, "linear_drift": 1000 + 40 * T, "curved": 1000 + 300 * np.sin(T / 5)}
frames = []
for name, base in BASELINES.items():
    for rep in range(1, 6):
        frames.append(trace(base + gaussian(3000, 6.0, 0.2) + rng.normal(0, 50, T.size),
                            baseline=name, replicate=rep))
    truth.append({"case": 2, "scenario": name, "retention_time": 6.0, "area": 3000, "sigma": 0.2})
pd.concat(frames).to_csv(OUT / "case2_baselines.csv", index=False)

# Case 3: calibration standards and two unknowns (true area = 30 x concentration)
frames = []
for kind, conc in [("standard", c) for c in (1, 5, 10, 50, 100)] + [("unknown", 20), ("unknown", 75)]:
    frames.append(trace(1000 + gaussian(30 * conc, 5.2, 0.2) + rng.normal(0, 30, T.size),
                        sample_type=kind, concentration=conc))
    truth.append({"case": 3, "scenario": f"{kind}_{conc}", "retention_time": 5.2, "area": 30 * conc, "sigma": 0.2})
pd.concat(frames).to_csv(OUT / "case3_calibration.csv", index=False)

# Case 4: peak pairs at four separations, and one peak at three tailing levels
frames = []
for dt in (0.3, 0.5, 0.8, 1.2):
    y = 1000 + gaussian(1000, 6.0, 0.1) + gaussian(1000, 6.0 + dt, 0.1) + rng.normal(0, 10, T.size)
    frames.append(trace(y, scenario=f"separation_{dt}"))
    truth.append({"case": 4, "scenario": f"separation_{dt}", "retention_time": 6.0, "area": 1000, "sigma": 0.1,
                  "true_resolution": round(dt / (4 * 0.1), 3)})
for k in (0.01, 1.0, 2.0):
    y = 1000 + 1000 * exponnorm.pdf(T, k, loc=6.0, scale=0.1) + rng.normal(0, 5, T.size)
    frames.append(trace(y, scenario=f"tailing_K{k}"))
    truth.append({"case": 4, "scenario": f"tailing_K{k}", "retention_time": 6.0, "area": 1000, "sigma": 0.1,
                  "tailing_K": k})
pd.concat(frames).to_csv(OUT / "case4_separation_and_shape.csv", index=False)

pd.DataFrame(truth).to_csv(OUT / "truth.csv", index=False)
print(f"Wrote case-study data to {OUT}")
