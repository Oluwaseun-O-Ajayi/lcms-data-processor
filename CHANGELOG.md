# Changelog

## Unreleased

### Added
- `lcms-batch --peak-mode sum`: integrate every peak inside a target's retention window (for partly resolved isomers).
- `PeakDetector.detect_peaks(threshold_range=...)`: relative detection thresholds scaled to the target's retention window; `lcms-batch` uses it whenever a window is given, so a small peak next to a much larger peak at the same m/z is still detected (found by cross-checking against Skyline).
- `glycan_targets.parse_adducts` and `glycan_mz_adducts`: mixed adduct sets such as `H+Na` and `2H+Na`.
- `lcms-glycan-search`: composition search for unassigned peaks, and `--survey` for a whole-run composition survey with isotope and monoisotopic checks.
- `examples/mass_errors.py`: measured ppm error of each target ion.
- Fluorescence comparison matches every MS apex of a glycan, sums the matched fluorescence peaks, restricts fluorescence peak detection to the glycan elution region, and writes the detected fluorescence peaks.
- Skyline comparison reads the charge from Skyline's `Precursor Adduct` column and accepts the `Molecule` header.
- Cluster job files for the fetuin case study (`cluster/sweep.slurm`, `cluster/extras.slurm`) and target lists.

## 1.3.0 (2026-11)

### Added
- `glycan_targets`: monoisotopic m/z for released N-glycans from compositions (long form `HexNAc4Hex5dHex1NeuAc2` or short form `N4H5F1S2`), with reducing-end labels (none, reduced, 2-AB, 2-AA, procainamide), H+, Na+ or NH4+ adducts and any charge state. Masses are computed from element formulas and checked against tabulated residue masses.
- `MzMLReader.xic_matrix`: extracts many ion chromatograms in one pass over the spectra.
- `lcms-batch` command: processes mzML files over a grid of extraction windows, peak widths and baseline windows, picks the largest peak in an optional retention-time window, and writes per-target results and per-glycan relative abundances.
- `reference_compare` and the `lcms-compare` command: compare MS relative abundances with fluorescence (FLD trace read from the mzML, transfer delay estimated by cross-correlation, co-eluting glycans flagged) and with Skyline (writes a transition list, reads a Skyline report); reports mean absolute difference, Bland-Altman bias and limits of agreement, and correlations.
- `MzMLReader.list_chromatograms` and `MzMLReader.chromatogram`: read detector traces and other stored chromatograms.
- `examples/summarize_sweep.py`: spread of relative abundance across settings, with a figure of abundance against extraction window.
- `cluster/`: download script with checksum verification for a public released N-glycan dataset (doi:10.5281/zenodo.4705546), conversion notes, and a SLURM array job for parameter sweeps.

### Changed
- `ChromatogramParser` raises `ValueError` (chained to the original error) for unreadable or malformed CSV files.

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
