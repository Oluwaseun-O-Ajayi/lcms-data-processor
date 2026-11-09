# Changelog

## 1.2.0 (2026-11)

### Added
- `MzMLReader`: reads mzML with pyOpenMS and returns total-ion and extracted-ion chromatograms (m/z ± ppm, any MS level) as time-intensity tables that `PeakDetector` accepts directly. Retention times are converted from seconds to minutes.
- `examples/make_mzml_example.py` and `examples/mzml_xic_workflow.py`: a synthetic mzML run with known areas and a co-eluting ion 35 ppm from one analyte, processed at 10 and 50 ppm windows.
- `tests/test_mzml.py`: accuracy tests for mzML input (skipped when pyOpenMS is not installed).
- Optional dependency group `mzml` (`pip install ".[mzml]"`).

## 1.1.1 (2026-10)

### Changed
- Paper: added an AI usage disclosure.

## 1.1.0 (2026-10)

### Fixed
- Baseline: the rolling-minimum baseline followed the lower edge of the noise and inflated peak areas (about +5% at a noise standard deviation of 500 on the test chromatograms). The baseline now uses a morphological opening with a noise-offset correction and a second pass that excludes detected peaks.
- Asymmetry factor: version 1.0 placed the apex at the midpoint of the integration limits, so every peak scored about 1.0. Asymmetry is now measured at 10% of height from the true apex.
- Resolution: version 1.0 divided by an arbitrary factor and mixed units (points and minutes). Resolution now uses half-height widths in time units, Rs = 1.18 (t2 - t1) / (w1 + w2).
- Signal-to-noise: noise was the standard deviation of the first 100 raw points; it is now a robust estimate from the whole baseline-corrected trace.
- `examples/calibration_workflow.py` crashed when printing the slope; fixed, and both examples now run from any working directory.
- Plotting functions no longer call `plt.show()` (which blocked batch runs) and the module no longer changes global warning filters or plot styles on import.

### Changed
- Peak detection uses a Savitzky-Golay smoothed trace with noise-aware height and prominence thresholds; integration limits follow each peak's own half-widths and stop at the valley between neighbouring peaks.
- `CalibrationCurve` requires at least three standards and exposes `slope` and `intercept`.
- The example data generator uses a fixed seed and writes to the repository's `data/` folder.
- Installable with `pip install .`; seaborn is no longer a dependency.

### Added
- Accuracy tests on simulated chromatograms with known answers, and continuous integration on Linux, macOS and Windows.
- `CONTRIBUTING.md`, this changelog and `paper/make_figure.py`.

## 1.0.0 (2026-07)
- Initial archived release.
