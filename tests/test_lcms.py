"""Accuracy tests for LC-MS Data Processor on synthetic chromatograms with known answers."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import exponnorm

from lcms_data_processor import CalibrationCurve, ChromatogramParser, PeakDetector, QualityControl

ROOT = Path(__file__).resolve().parents[1]
T = np.linspace(0, 20, 4000)


def gaussian(area, center, sigma):
    return area / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((T - center) / sigma) ** 2)


def detect(y):
    det = PeakDetector(pd.DataFrame({"time": T, "intensity": y}))
    det.detect_peaks()
    det.integrate_peaks()
    det.calculate_snr()
    return det


@pytest.mark.parametrize("noise", [5, 50, 500])
def test_area_of_a_known_peak_is_recovered_on_a_drifting_baseline(noise):
    rng = np.random.default_rng(0)
    ratios = []
    for _ in range(10):
        y = 1000 + 2 * T + gaussian(1e5, 8.0, 0.1) + rng.normal(0, noise, T.size)
        ratios.append(detect(y).peaks["peak_area"].iloc[0] / 1e5)
    assert np.mean(ratios) == pytest.approx(1.0, abs=0.005)


def test_area_is_unbiased_at_moderate_signal_to_noise():
    rng = np.random.default_rng(1)
    ratios = []
    for _ in range(30):
        y = 1000 + gaussian(376, 8.0, 0.2) + rng.normal(0, 30, T.size)  # height/noise = 25
        peaks = detect(y).peaks
        ratios.append(peaks["peak_area"].iloc[(peaks["retention_time"] - 8).abs().argmin()] / 376)
    assert np.mean(ratios) == pytest.approx(1.0, abs=0.03)


def test_noise_estimate_and_signal_to_noise():
    rng = np.random.default_rng(2)
    y = 1000 + gaussian(1e4, 8.0, 0.1) + rng.normal(0, 20, T.size)
    det = detect(y)
    assert det.noise_sd == pytest.approx(20, rel=0.15)
    height = 1e4 / (0.1 * np.sqrt(2 * np.pi))
    assert det.peaks["snr"].iloc[0] == pytest.approx(height / 20, rel=0.15)


def test_fwhm_matches_gaussian_theory():
    y = 1000 + gaussian(1e4, 8.0, 0.2) + np.random.default_rng(3).normal(0, 1, T.size)
    assert detect(y).peaks["fwhm"].iloc[0] == pytest.approx(2.3548 * 0.2, rel=0.02)


def test_resolution_of_two_gaussians():
    # Rs = (t2 - t1) / (4 sigma) = 0.6 / 0.4 = 1.5 for equal sigma
    y = 1000 + gaussian(1e3, 8.0, 0.1) + gaussian(1e3, 8.6, 0.1) + np.random.default_rng(4).normal(0, 1, T.size)
    det = detect(y)
    res = QualityControl().check_resolution(det.peaks)
    assert res["resolution"].iloc[0] == pytest.approx(1.5, rel=0.05)


def test_overlapping_peaks_are_split_at_the_valley():
    y = 1000 + gaussian(1e3, 8.0, 0.1) + gaussian(5e2, 8.5, 0.1) + np.random.default_rng(5).normal(0, 1, T.size)
    areas = detect(y).peaks["peak_area"].to_numpy()
    np.testing.assert_allclose(areas, [1e3, 5e2], rtol=0.03)


def test_asymmetry_factor_detects_tailing():
    qc = QualityControl()
    values = []
    for k in (1e-4, 2.0):
        shape = exponnorm.pdf(T, k, loc=8.0, scale=0.1)
        y = 1000 + 1e4 * shape / shape.max() + np.random.default_rng(6).normal(0, 2, T.size)
        det = detect(y)
        values.append(qc.check_peak_shape(det.peaks, det.data, det.baseline)["asymmetry_factor"].iloc[0])
    assert values[0] == pytest.approx(1.0, abs=0.1)
    assert values[1] > 1.5


def test_calibration_recovers_slope_and_quantifies():
    conc = np.array([1, 5, 10, 50, 100.0])
    areas = 75.2 * conc + 3.0
    cal = CalibrationCurve()
    _, r2 = cal.fit_calibration(conc, areas)
    assert cal.slope == pytest.approx(75.2) and cal.intercept == pytest.approx(3.0)
    assert r2 == pytest.approx(1.0)
    np.testing.assert_allclose(cal.quantify_samples([75.2 * 20 + 3.0]), [20.0])
    with pytest.raises(ValueError):
        CalibrationCurve().fit_calibration([1, 2], [1, 2])


def test_parser_validates_columns(tmp_path):
    bad = tmp_path / "bad.csv"
    pd.DataFrame({"t": [0, 1], "i": [1, 2]}).to_csv(bad, index=False)
    with pytest.raises(Exception):
        ChromatogramParser(bad)
    with pytest.raises(FileNotFoundError):
        ChromatogramParser(tmp_path / "missing.csv")


@pytest.mark.parametrize("script", ["generate_sample_data.py", "examples/basic_analysis.py",
                                    "examples/calibration_workflow.py"])
def test_scripts_run(script, tmp_path):
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONIOENCODING="utf-8")
    proc = subprocess.run([sys.executable, str(ROOT / script)], cwd=tmp_path, env=env,
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-2000:]
