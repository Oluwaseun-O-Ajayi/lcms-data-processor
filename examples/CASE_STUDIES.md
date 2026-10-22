# Case studies

Four worked examples that check the processing against known answers. Each runs
on synthetic chromatograms bundled in `data/case_studies/`, whose true peak areas,
retention times and widths are listed in `data/case_studies/truth.csv`. Results are
written to `examples/case_study_results/`. Reproduce them with:

```bash
python examples/case_study_1_area_vs_noise.py
python examples/case_study_2_baseline_drift.py
python examples/case_study_3_calibration.py
python examples/case_study_4_resolution_and_shape.py
```

The bundled data were created by `examples/make_case_study_data.py` (fixed seed).
The numbers below were produced with version 1.1.1.

| Case | Application | Input |
|---|---|---|
| 1 | Peak area accuracy as noise increases | Four peaks at noise SD 20, 100 and 300 |
| 2 | Area accuracy on non-flat baselines | One peak on flat, drifting and curved baselines, five replicates each |
| 3 | Calibration and quantification of unknowns | Standards at 1 to 100 ng/mL, unknowns at 20 and 75 ng/mL |
| 4 | Resolution and peak shape against theory | Peak pairs at four separations; one peak at three tailing levels |

## Case 1: Peak areas at three noise levels

**Question.** How accurately are known peak areas recovered as noise increases?

**Result.** At noise SD 20, all four areas are within 1% of the truth (99.5 to
100.7%). At noise SD 100, recoveries range from 96.1 to 103.4%, and at noise SD 300
from 98.6 to 104.6%. The smallest peak (true area 500) is still detected at a
signal-to-noise ratio of 4.0.

**What it shows.** Integration is unbiased when the signal is strong. As noise
grows, individual areas scatter by a few percent, which sets the precision a user
can expect from a single injection at that signal-to-noise ratio.

## Case 2: The same peak on different baselines

**Question.** Does baseline correction recover the true area when the baseline
drifts or curves?

**Result.** Mean recovery over five replicates is 99.5% on a flat baseline, 99.5%
on a linear drift and 100.0% on a curved baseline, with standard deviations of 0.9
to 1.3%.

**What it shows.** A drifting or curved baseline does not bias the area, because
the baseline is re-estimated from the points outside each peak window. The
remaining spread comes from noise, not from the baseline shape.

## Case 3: Calibration and quantification

**Question.** How accurately are unknowns quantified, and how do standards behave
near the detection limit?

**Result.** The 1 ng/mL standard (signal-to-noise about 2) is not detected, so the
calibration uses the four standards from 5 to 100 ng/mL: area = 29.68 x
concentration + 10.85, R^2 = 1.00000, against a true slope of 30. The unknowns
are back-calculated as 19.83 ng/mL (true 20, 99.1%) and 75.28 ng/mL (true 75,
100.4%). The 5 ng/mL standard back-calculates to 4.93 ng/mL (98.6%) at a
signal-to-noise ratio of 9.9.

**What it shows.** Quantification within 1% is achievable inside the calibrated
range. A standard below the detection limit is reported as missing rather than
given an invented area, which is how it should be handled in a real calibration.

## Case 4: Resolution and peak shape against theory

**Question.** Does the reported resolution match theory, and does the asymmetry
factor rise with tailing?

**Result.**

| Separation (min) | Theoretical Rs | Reported Rs | Classification |
|---|---|---|---|
| 0.3 | 0.75 | 0.53 | Poor |
| 0.5 | 1.25 | 1.29 | Partial |
| 0.8 | 2.00 | 2.00 | Baseline |
| 1.2 | 3.00 | 3.00 | Baseline |

| Tailing (exponential modifier K) | Asymmetry factor | Classification |
|---|---|---|
| 0.01 (symmetric) | 1.01 | Good |
| 1.0 | 1.34 | Acceptable |
| 2.0 | 1.91 | Poor |

**What it shows.** Resolution matches theory once peaks are reasonably separated.
For strongly overlapping peaks it is underestimated, because the half-height
widths are measured on the merged profile, so a low reported value should be read
as "poorly resolved" rather than as an exact number. The asymmetry factor rises
with tailing as expected.
