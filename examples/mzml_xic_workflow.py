"""
From mzML to quantified peaks: extract ion chromatograms with pyOpenMS and
integrate them with LC-MS Data Processor.

Runs on the synthetic file written by examples/make_mzml_example.py, whose true
peak areas are known, and prints recovered areas for a 10 ppm and a 50 ppm
extraction window. The 50 ppm window also captures a co-eluting ion 35 ppm from
analyte A and overstates its area; the 10 ppm window excludes it.

Requires: pip install pyopenms
"""

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lcms_data_processor import MzMLReader, PeakDetector  # noqa: E402
from make_mzml_example import ANALYTES, write_example  # noqa: E402


def main():
    with tempfile.TemporaryDirectory() as tmp:
        path = write_example(Path(tmp) / 'example_run.mzML')
        reader = MzMLReader(path)
        print(f"{reader.metadata['filename']}: {reader.metadata['num_spectra']} spectra, "
              f"MS levels {reader.metadata['ms_levels']}")

        for ppm in (10, 50):
            print(f"\nExtraction window +/- {ppm} ppm")
            print(f"{'analyte':<8}{'m/z':>11}{'RT (min)':>10}{'area':>12}{'true':>12}{'recovery':>10}")
            for name, xic in reader.xics({k: v['mz'] for k, v in ANALYTES.items()}, ppm=ppm).items():
                det = PeakDetector(xic)
                det.detect_peaks()
                det.integrate_peaks()
                peak = det.peaks.loc[det.peaks['peak_area'].idxmax()]
                true = ANALYTES[name]['area']
                print(f"{name:<8}{ANALYTES[name]['mz']:>11.4f}{peak['retention_time']:>10.2f}"
                      f"{peak['peak_area']:>12.0f}{true:>12.0f}{100 * peak['peak_area'] / true:>9.1f}%")


if __name__ == '__main__':
    main()
