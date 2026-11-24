"""The fluorescence comparison must detect small peaks down to the --fld-prominence setting,
not only peaks above the 5% height floor used for MS chromatograms."""
import numpy as np
import pandas as pd

from lcms_data_processor import PeakDetector


def _trace():
    t = np.arange(0, 30, 1 / 600)            # 10 Hz, minutes
    y = np.zeros_like(t)
    for rt, h in ((8.0, 100.0), (14.0, 4.0), (20.0, 1.0)):   # main peak, 4% peak, 1% peak
        y += h * np.exp(-0.5 * ((t - rt) / 0.08) ** 2)
    rng = np.random.default_rng(0)
    y += 1000 + rng.normal(0, 0.02, len(t))
    return pd.DataFrame({"time": t, "intensity": y})


def _found(prom, frac):
    d = PeakDetector(_trace())
    d.detect_peaks(prominence=prom, width=5, baseline_window=0.1, height_fraction=frac)
    return sorted(round(r) for r in d.peaks["retention_time"])


def test_default_height_floor_drops_small_peaks():
    assert _found(0.005, 0.05) == [8]          # the old behaviour: 4% and 1% peaks are lost


def test_height_floor_follows_prominence():
    assert _found(0.005, 0.005) == [8, 14, 20]
    assert _found(0.02, 0.02) == [8, 14]
