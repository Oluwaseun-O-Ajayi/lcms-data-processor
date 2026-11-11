"""Tests for mzML input: extracted-ion chromatograms from a synthetic run with known answers."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("pyopenms")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples"))

from lcms_data_processor import MzMLReader, PeakDetector  # noqa: E402
from make_mzml_example import ANALYTES, INTERFERENT, write_example  # noqa: E402


@pytest.fixture(scope="module")
def reader(tmp_path_factory):
    return MzMLReader(write_example(tmp_path_factory.mktemp("mzml") / "run.mzML"))


def integrate(xic):
    det = PeakDetector(xic)
    det.detect_peaks()
    det.integrate_peaks()
    return det.peaks.loc[det.peaks["peak_area"].idxmax()]


def test_metadata_and_time_units(reader):
    assert reader.metadata["ms_levels"] == [1]
    tic = reader.tic()
    assert len(tic) == reader.metadata["num_spectra"] == 1800
    assert tic["time"].iloc[0] == pytest.approx(0.0)
    assert tic["time"].iloc[-1] == pytest.approx(15.0 - 0.5 / 60)  # seconds converted to minutes


@pytest.mark.parametrize("name", list(ANALYTES))
def test_xic_recovers_known_area_and_retention_time(reader, name):
    truth = ANALYTES[name]
    peak = integrate(reader.xic(truth["mz"], ppm=10))
    assert peak["retention_time"] == pytest.approx(truth["rt_min"], abs=0.02)
    assert peak["peak_area"] / truth["area"] == pytest.approx(1.0, abs=0.02)


def test_narrow_window_excludes_nearby_ion_and_wide_window_includes_it(reader):
    a = ANALYTES["A"]
    narrow = integrate(reader.xic(a["mz"], ppm=10))["peak_area"]
    wide = integrate(reader.xic(a["mz"], ppm=50))["peak_area"]
    assert narrow / a["area"] == pytest.approx(1.0, abs=0.02)
    assert wide / (a["area"] + INTERFERENT["area"]) == pytest.approx(1.0, abs=0.02)


def test_xic_away_from_analytes_is_background_only(reader):
    xic = reader.xic(1500.0, ppm=10)
    assert xic["intensity"].max() < 0.01 * ANALYTES["B"]["area"]


def test_xics_returns_one_chromatogram_per_target(reader):
    out = reader.xics({k: v["mz"] for k, v in ANALYTES.items()})
    assert set(out) == set(ANALYTES)
    assert all(list(df.columns) == ["time", "intensity"] for df in out.values())


def test_invalid_input(reader, tmp_path):
    with pytest.raises(ValueError):
        reader.xic(-1)
    with pytest.raises(ValueError):
        reader.xic(500.0, ms_level=2)
    with pytest.raises(FileNotFoundError):
        MzMLReader(tmp_path / "missing.mzML")


def test_mzml_example_runs(tmp_path):
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONIOENCODING="utf-8")
    proc = subprocess.run([sys.executable, str(ROOT / "examples" / "mzml_xic_workflow.py")],
                          cwd=tmp_path, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=300)
    assert proc.returncode == 0, proc.stderr[-2000:]
