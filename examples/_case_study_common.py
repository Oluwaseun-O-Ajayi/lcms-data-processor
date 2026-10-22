"""Shared helpers for the case-study scripts."""

from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from lcms_data_processor import PeakDetector  # noqa: E402

DATA = ROOT / "data" / "case_studies"
RESULTS = Path(__file__).resolve().parent / "case_study_results"
RESULTS.mkdir(exist_ok=True)
TRUTH = pd.read_csv(DATA / "truth.csv")


def analyze(trace):
    """Detect, integrate and compute S/N for one chromatogram (DataFrame with time, intensity)."""
    det = PeakDetector(trace[["time", "intensity"]].reset_index(drop=True))
    det.detect_peaks()
    det.integrate_peaks()
    det.calculate_snr()
    return det


def nearest(peaks, retention_time, tolerance=0.1):
    """The detected peak closest to retention_time, or None if none lies within tolerance."""
    if peaks is None or len(peaks) == 0:
        return None
    row = peaks.iloc[(peaks["retention_time"] - retention_time).abs().argmin()]
    return row if abs(row["retention_time"] - retention_time) <= tolerance else None
