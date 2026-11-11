"""Tests for glycan m/z targets, multi-target extraction and the batch command."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import glycan_targets as gt


# Tabulated monoisotopic residue masses
@pytest.mark.parametrize("name, mass", [("Hex", 162.0528), ("HexNAc", 203.0794),
                                        ("dHex", 146.0579), ("NeuAc", 291.0954), ("NeuGc", 307.0903)])
def test_residue_masses(name, mass):
    assert gt.RESIDUE[name] == pytest.approx(mass, abs=1e-4)


def test_g0f_free_and_2ab_masses():
    # G0F = HexNAc4 Hex3 dHex1: 4(203.07937) + 3(162.05282) + 146.05791 + 18.01056
    assert gt.glycan_mass("N4H3F1") == pytest.approx(1462.5444, abs=1e-3)
    assert gt.glycan_mz("N4H3F1", 1, label="2-AB") == pytest.approx(1583.6204, abs=1e-3)
    assert gt.glycan_mz("N4H3F1", 2, label="2-AB") == pytest.approx(792.3138, abs=1e-3)


def test_label_shift_matches_elemental_formula():
    # 2-AB G0F as an elemental formula: C56H90N6O39 (glycan C49H82N4O39 + C7H8N2O - H2O + 2H)
    glycan = {"C": 8 * 4 + 6 * 3 + 6, "H": 13 * 4 + 10 * 3 + 10 + 2, "N": 4, "O": 5 * 4 + 5 * 3 + 4 + 1}
    labelled = {k: glycan.get(k, 0) + gt.LABEL_FORMULA["2-AB"].get(k, 0) for k in "CHNO"}
    labelled["O"] -= 1  # - H2O + 2H leaves the H count unchanged
    assert gt.glycan_mass("N4H3F1", label="2-AB") == pytest.approx(gt.formula_mass(labelled), abs=1e-6)


def test_long_and_short_forms_agree():
    assert gt.parse_composition("HexNAc4Hex5dHex1NeuAc2") == gt.parse_composition("N4H5F1S2")


@pytest.mark.parametrize("bad", ["N4X3", "N4H3N1", "", "hello"])
def test_bad_compositions_raise(bad):
    with pytest.raises(ValueError):
        gt.parse_composition(bad)


def test_unknown_label_raises():
    with pytest.raises(ValueError):
        gt.label_shift("RapiFluor")


def test_target_table_shape():
    t = gt.target_table(gt.FETUIN_CANDIDATES, label="2-AB", charges=(2, 3))
    assert len(t) == 2 * len(gt.FETUIN_CANDIDATES)
    assert (t.loc[t.charge == 3, "mz"].to_numpy() < t.loc[t.charge == 2, "mz"].to_numpy()).all()


# --- mzML-based tests -------------------------------------------------------

pytest.importorskip("pyopenms")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
from make_mzml_example import ANALYTES, write_example  # noqa: E402
from lcms_data_processor import MzMLReader  # noqa: E402
import lcms_batch  # noqa: E402


@pytest.fixture(scope="module")
def mzml(tmp_path_factory):
    return write_example(tmp_path_factory.mktemp("batch") / "run.mzML")


@pytest.mark.parametrize("ppm", [10, 50])
def test_xic_matrix_matches_single_xic(mzml, ppm):
    reader = MzMLReader(mzml)
    mzs = [v["mz"] for v in ANALYTES.values()] + [1500.0]
    times, mat = reader.xic_matrix(mzs, ppm=ppm)
    for j, mz in enumerate(mzs):
        single = reader.xic(mz, ppm=ppm)
        np.testing.assert_allclose(times, single["time"].to_numpy())
        np.testing.assert_allclose(mat[:, j], single["intensity"].to_numpy(), rtol=1e-6, atol=1e-3)  # mzML intensities are float32


def test_batch_sweep_recovers_areas_and_abundances(mzml, tmp_path):
    targets = pd.DataFrame({"name": list(ANALYTES), "mz": [v["mz"] for v in ANALYTES.values()],
                            "rt_min": [v["rt_min"] - 0.5 for v in ANALYTES.values()],
                            "rt_max": [v["rt_min"] + 0.5 for v in ANALYTES.values()]})
    tfile = tmp_path / "targets.csv"
    targets.to_csv(tfile, index=False)
    long, per = lcms_batch.run(lcms_batch.build_parser().parse_args(
        ["--mzml", str(mzml), "--targets", str(tfile), "--ppm", "10", "50",
         "--width", "5", "--baseline-window", "0.05", "0.1", "--out", str(tmp_path / "sweep")]))
    assert len(long) == 2 * 2 * len(ANALYTES)
    assert (tmp_path / "sweep_glycans.csv").exists()

    at10 = long[long.ppm == 10].groupby("target")["peak_area"].mean()
    for name, truth in ANALYTES.items():
        assert at10[name] / truth["area"] == pytest.approx(1.0, abs=0.02)

    # Relative abundance of A at 10 ppm is its true share; at 50 ppm the interferent inflates it.
    share = ANALYTES["A"]["area"] / sum(v["area"] for v in ANALYTES.values())
    rel = per[per.glycan == "A"].groupby("ppm")["relative_abundance_pct"].mean()
    assert rel[10] == pytest.approx(100 * share, abs=1.0)
    assert rel[50] > rel[10] + 2


def test_batch_glycan_preset_runs(mzml, tmp_path):
    long, per = lcms_batch.run(lcms_batch.build_parser().parse_args(
        ["--mzml", str(mzml), "--glycans", "fetuin", "--charges", "2", "3", "--out", str(tmp_path / "g")]))
    assert set(per["glycan"]) == set(gt.FETUIN_CANDIDATES)


def test_summarize_sweep_writes_tables_and_figure(mzml, tmp_path):
    lcms_batch.run(lcms_batch.build_parser().parse_args(
        ["--mzml", str(mzml), "--glycans", "igg", "--ppm", "5", "50", "--out", str(tmp_path / "s")]))
    import summarize_sweep
    summarize_sweep.main([str(tmp_path / "s_glycans.csv"), "--out", str(tmp_path / "sum")])
    for suffix in ("_spread.csv", "_by_ppm.csv", "_abundance_vs_ppm.png", "_abundance_vs_ppm.pdf"):
        assert (tmp_path / f"sum{suffix}").exists()
