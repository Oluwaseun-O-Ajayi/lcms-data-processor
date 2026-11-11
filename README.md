# LC-MS Data Processor

**Transparent processing of LC-MS chromatograms, from exported CSV traces or mzML files: baseline, peak detection and integration, signal-to-noise, peak shape, resolution and linear calibration.**

[![tests](https://github.com/Oluwaseun-O-Ajayi/lcms-data-processor/actions/workflows/tests.yml/badge.svg)](https://github.com/Oluwaseun-O-Ajayi/lcms-data-processor/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21652569.svg)](https://doi.org/10.5281/zenodo.21652569)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

LC-MS Data Processor reads a chromatogram saved as a CSV file with `time` and `intensity` columns (for example, an extracted-ion or total-ion chromatogram exported from vendor software) and runs a short, readable pipeline: baseline estimation, peak detection, integration, signal-to-noise, asymmetry and resolution checks, and linear calibration for quantification. Every step is a few lines of documented Python, so the processing can be inspected, taught and changed.

From version 1.2.0 it also reads mzML through [pyOpenMS](https://pyopenms.readthedocs.io/) and extracts total-ion and extracted-ion chromatograms (XICs) with a ppm window, so open-format instrument data can go straight into the same pipeline.

## Scope and limitations

- Works on time-intensity traces, either exported as CSV or extracted from mzML (MS1 or any MS level) by summing intensity within a ppm window. It does not read vendor raw files (convert them to mzML first, for example with ProteoWizard msconvert), and it does not perform feature finding, alignment, deconvolution or identification.
- Intended for teaching, method-development prototyping and transparent reprocessing. It is not a validated replacement for vendor or regulated quantitative software; users must validate settings and results for their own methods.
- The bundled chromatograms in `data/` and the mzML example are synthetic (Gaussian peaks, flat baseline, random noise, fixed seed) so that every number has a known answer. They are not experimental measurements.

## Method

| Step | Implementation |
|---|---|
| mzML input | `MzMLReader` loads the file with pyOpenMS, converts retention times from seconds to minutes and sums intensity within m/z ± ppm in each spectrum |
| Baseline | Morphological opening (rolling minimum, then rolling maximum) with an offset correction for noise, refined in a second pass by a rolling median of the points outside detected peaks, interpolated across the peaks |
| Noise | Robust standard deviation from first differences, 1.4826 x MAD / sqrt(2) |
| Detection | `scipy.signal.find_peaks` on a Savitzky-Golay smoothed trace; height and prominence thresholds never below three noise standard deviations |
| Integration | Trapezoidal rule between limits set 3.4 half-widths at half maximum from the apex on each side (four standard deviations for a Gaussian), never crossing the valley between neighbouring peaks |
| Signal-to-noise | Corrected peak height divided by the noise standard deviation (pharmacopeial definitions based on peak-to-peak noise give different values) |
| Asymmetry | As = b / a at 10% of peak height |
| Resolution | Rs = 1.18 (t2 - t1) / (w0.5,1 + w0.5,2) from half-height widths |
| Calibration | Ordinary least-squares line of area against concentration; inverse prediction for unknowns |

On simulated Gaussian peaks on a drifting baseline, the mean recovered area is within 1.1% of the true area for signal-to-noise ratios from 5 to 400; the scatter grows as signal-to-noise falls (`paper/make_figure.py` reproduces this check).

## Installation

Requires Python 3.9 or later.

```bash
git clone https://github.com/Oluwaseun-O-Ajayi/lcms-data-processor.git
cd lcms-data-processor
pip install .            # or: pip install ".[test]" to include the test dependencies
python generate_sample_data.py
python examples/basic_analysis.py
python examples/calibration_workflow.py

pip install ".[mzml]"    # optional: adds pyOpenMS for mzML input
python examples/mzml_xic_workflow.py
```

Four reproducible case studies that check areas, baseline correction, calibration and resolution against known answers are described in [`examples/CASE_STUDIES.md`](examples/CASE_STUDIES.md).

## Quick start

```python
from lcms_data_processor import ChromatogramParser, PeakDetector, QualityControl

chrom = ChromatogramParser("data/sample1.csv").get_chromatogram()
det = PeakDetector(chrom)
det.detect_peaks()
det.integrate_peaks()
det.calculate_snr()
print(det.peaks[["peak_id", "retention_time", "peak_area", "fwhm", "snr"]])

qc = QualityControl()
print(qc.check_peak_shape(det.peaks, chrom, det.baseline))
print(qc.check_resolution(det.peaks))
```

Released N-glycans, over a grid of processing settings:

```bash
pip install ".[mzml]"
lcms-batch --mzml run.mzML --glycans fetuin --label 2-AB --charges 1 2 3 \
           --ppm 5 10 20 50 --width 3 5 --baseline-window 0.05 0.1 --out results/run
python examples/summarize_sweep.py results/run_glycans.csv --out results/summary
```

`lcms-compare` checks MS-based abundances against fluorescence (for 2-AB and other fluorescent labels) and against Skyline.

`glycan_targets` computes candidate m/z values from compositions (for example `N4H5F1S2`) with reducing-end labels and charge states; assignments still need retention time and MS/MS evidence. `cluster/README.md` shows how to run the sweep on a SLURM cluster with a public fetuin and IgG dataset.

From mzML:

```python
from lcms_data_processor import MzMLReader, PeakDetector

reader = MzMLReader("run.mzML")
xic = reader.xic(1111.387, ppm=10)          # time (min), intensity
det = PeakDetector(xic)
det.detect_peaks()
det.integrate_peaks()
```

`examples/mzml_xic_workflow.py` writes a synthetic run with known peak areas and a co-eluting ion 35 ppm from one analyte. A 10 ppm window recovers both true areas to within 0.1%; a 50 ppm window includes the interfering ion and overstates that analyte's area by 25%, which shows why the extraction window has to be reported with the result.

`LCMSAnalysisPipeline(output_dir="results").process_sample(path, name)` runs the same steps and writes plots, a peak table and a QC report.

## Testing

```bash
pip install ".[test]"
python -m pytest
```

The tests use simulated chromatograms with known answers: recovered peak areas on a drifting baseline at three noise levels, unbiased areas at moderate signal-to-noise, the noise and signal-to-noise estimates, full width at half maximum against Gaussian theory, the resolution of two Gaussians with Rs = 1.5, splitting of overlapping peaks, the asymmetry of symmetric and tailing peaks, calibration slope and inverse prediction, input validation, and that all scripts run. When pyOpenMS is installed, further tests check mzML reading, time units, recovered XIC areas and retention times, and the effect of the ppm window on a nearby co-eluting ion. Tests run on Linux, macOS and Windows for every push.

## Contributing and support

Bug reports, questions and pull requests are welcome; see [CONTRIBUTING.md](CONTRIBUTING.md). Please use [GitHub Issues](https://github.com/Oluwaseun-O-Ajayi/lcms-data-processor/issues); the author aims to respond within two weeks.

## Citation

Ajayi, O. O. (2026). *LC-MS Data Processor*. Zenodo. https://doi.org/10.5281/zenodo.21652569 (see `CITATION.cff`).

## License

MIT. See [LICENSE](LICENSE).
