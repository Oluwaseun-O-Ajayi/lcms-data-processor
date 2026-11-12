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


# --- reference comparisons --------------------------------------------------

import reference_compare as rc  # noqa: E402
from make_mzml_example import FLD_DELAY_MIN  # noqa: E402


@pytest.fixture(scope="module")
def fld_run(tmp_path_factory):
    d = tmp_path_factory.mktemp("fld")
    path = write_example(d / "fld.mzML", fld=True)
    targets = pd.DataFrame({"name": list(ANALYTES), "mz": [v["mz"] for v in ANALYTES.values()],
                            "rt_min": [v["rt_min"] - 0.5 for v in ANALYTES.values()],
                            "rt_max": [v["rt_min"] + 0.5 for v in ANALYTES.values()]})
    targets.to_csv(d / "t.csv", index=False)
    long, _ = lcms_batch.run(lcms_batch.build_parser().parse_args(
        ["--mzml", str(path), "--targets", str(d / "t.csv"), "--ppm", "10", "50", "--out", str(d / "s")]))
    return path, long


def test_fld_chromatogram_is_listed_and_read(fld_run):
    reader = MzMLReader(fld_run[0])
    listing = reader.list_chromatograms()
    assert listing["type"].str.contains("EMISSION").any()
    fld = reader.chromatogram("fld")
    assert fld["time"].max() == pytest.approx(15.0, abs=0.01)
    with pytest.raises(ValueError):
        reader.chromatogram("UV")


def test_fld_delay_estimate_and_agreement(fld_run):
    reader = MzMLReader(fld_run[0])
    table, s10 = rc.compare_fld(reader, fld_run[1], "FLD", ppm=10, width=5, baseline_window=0.1)
    assert s10["delay_min"] == pytest.approx(FLD_DELAY_MIN, abs=0.02)
    assert (table["status"] == "matched").all()
    assert s10["mean_abs_diff_pct_points"] < 1.0
    # The 50 ppm window lets the unlabelled interferent inflate A's MS share; FLD exposes it.
    _, s50 = rc.compare_fld(reader, fld_run[1], "FLD", ppm=50, width=5, baseline_window=0.1)
    assert s50["mean_abs_diff_pct_points"] > s10["mean_abs_diff_pct_points"] + 2


def test_coeluting_glycans_on_one_fld_peak_are_flagged_not_compared():
    t = np.arange(0, 10, 0.01)

    class FakeReader:
        def chromatogram(self, match):
            return pd.DataFrame({"time": t, "intensity": 1000 * np.exp(-0.5 * ((t - 4.9) / 0.05) ** 2)
                                 + 800 * np.exp(-0.5 * ((t - 7.9) / 0.05) ** 2)})

    ms = pd.DataFrame({"glycan": ["X", "Y", "Z"], "ppm": 10.0, "width": 5, "baseline_window": 0.1,
                       "peak_area": [2.0, 1.0, 3.0], "retention_time": [5.0, 5.02, 8.0]})
    table, s = rc.compare_fld(FakeReader(), ms, "FLD", 10, 5, 0.1, delay=0.1)
    status = table.set_index("glycan")["status"]
    assert status["X"] == status["Y"] == "shared FLD peak"
    assert status["Z"] == "matched"
    assert s["n"] == 1


def test_skyline_round_trip(tmp_path):
    tl = gt.target_table(gt.IGG_CANDIDATES, label="2-AB", charges=(1, 2))
    tl["name"] = tl["glycan"] + "_" + tl["charge"].astype(str) + "+"
    sky_list = rc.skyline_transition_list(tl)
    assert list(sky_list.columns[:4]) == ["Molecule List Name", "Molecule Name", "Precursor m/z", "Precursor Charge"]

    rng = np.random.default_rng(1)
    ours = tl.assign(ppm=10.0, width=5, baseline_window=0.1, target=tl["name"],
                     peak_area=rng.uniform(1e4, 1e6, len(tl)))
    report = pd.DataFrame({"Molecule Name": tl["glycan"], "Precursor Charge": tl["charge"],
                           "Replicate Name": "r1",
                           "Total Area MS1": ours["peak_area"] * 2.0 * rng.normal(1, 0.02, len(tl))})
    _, per, s = rc.compare_skyline(ours, report, 10, 5, 0.1)
    assert s["n_precursors"] == len(tl)
    assert s["median_area_ratio"] == pytest.approx(0.5, rel=0.05)
    assert s["spearman_rho_areas"] > 0.98
    assert s["mean_abs_diff_pct_points"] < 0.5


def test_sum_mode_adds_isomer_peaks_inside_window_only():
    t = np.arange(0, 20, 0.005)
    g = lambda mu, area, s=0.08: area / (s * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((t - mu) / s) ** 2)
    y = g(10.0, 3000) + g(10.6, 1000) + g(15.0, 5000)  # isomer pair at 10.0/10.6, other glycan at 15
    largest, _ = lcms_batch.best_peak(t, y, 9.5, 11.0, 5, 0.1, 0.05, "largest")
    summed, _ = lcms_batch.best_peak(t, y, 9.5, 11.0, 5, 0.1, 0.05, "sum")
    assert largest["peak_area"] == pytest.approx(3000, rel=0.03)
    assert summed["peak_area"] == pytest.approx(4000, rel=0.03)
    assert summed["n_summed"] == 2 and summed["retention_time"] == pytest.approx(10.0, abs=0.01)
    assert [float(x) for x in summed["peak_rts"].split(";")] == pytest.approx([10.0, 10.6], abs=0.01)


def test_fld_sums_peaks_for_each_ms_apex():
    t = np.arange(0, 20, 0.005)
    g = lambda mu, area, s=0.08: area / (s * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((t - mu) / s) ** 2)

    class FakeReader:
        def chromatogram(self, match):
            return pd.DataFrame({"time": t, "intensity": g(9.9, 30) + g(10.5, 10) + g(14.9, 60)})

    ms = pd.DataFrame({"glycan": ["P", "Q"], "ppm": 10.0, "width": 5, "baseline_window": 0.1,
                       "peak_area": [4.0, 6.0], "retention_time": [10.0, 15.0],
                       "peak_rts": ["10.0000;10.6000", "15.0000"]})
    table, s = rc.compare_fld(FakeReader(), ms, "FLD", 10, 5, 0.1, delay=0.1)
    row = table.set_index("glycan")
    assert row.loc["P", "fld_area"] == pytest.approx(40, rel=0.03)
    assert (table["status"] == "matched").all()
    assert s["mean_abs_diff_pct_points"] < 1.0  # MS 40/60 vs FLD 40/60


def test_large_free_label_peak_does_not_hide_glycan_fld_peaks():
    t = np.arange(0, 45, 0.005)
    g = lambda mu, area, s=0.08: area / (s * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((t - mu) / s) ** 2)

    class FakeReader:
        def chromatogram(self, match):  # free 2-AB at 2 min, 500x the glycan peaks
            return pd.DataFrame({"time": t, "intensity": g(2.0, 20000) + g(19.4, 30) + g(23.6, 10)})

    ms = pd.DataFrame({"glycan": ["S1", "S2"], "ppm": 10.0, "width": 5, "baseline_window": 0.1,
                       "peak_area": [3.0, 1.0], "retention_time": [19.43, 23.63]})
    table, s = rc.compare_fld(FakeReader(), ms, "FLD", 10, 5, 0.1, delay=0.03)
    assert (table["status"] == "matched").all()
    assert s["fld_region_min"] > 2.5
    assert s["mean_abs_diff_pct_points"] < 1.0


def test_glycan_search_finds_known_composition(tmp_path):
    import pyopenms as oms
    import glycan_search as gs
    target = gt.glycan_mz("N5H6S3", 3, label="2-AB")  # A3G3S3 3+
    exp = oms.MSExperiment()
    rng = np.random.default_rng(3)
    for rt_s in np.arange(0, 120, 1.0):
        s = oms.MSSpectrum()
        s.setRT(float(rt_s))
        s.setMSLevel(1)
        mz = np.sort(np.append(rng.uniform(400, 2000, 30), target * (1 + 2e-6)))
        inten = np.where(np.isclose(mz, target * (1 + 2e-6)), 1e5, rng.exponential(100, len(mz)))
        s.set_peaks((mz, inten))
        exp.addSpectrum(s)
    oms.MzMLFile().store(str(tmp_path / "s.mzML"), exp)
    spec = gs.summed_spectrum(MzMLReader(tmp_path / "s.mzML"), 0.5, 1.5)
    res = gs.search(spec, "2-AB", (1, 2, 3), ppm=10, top=3)
    best = res.iloc[0]
    assert best["composition"] == "N5H6S3" and best["charge"] == 3
    assert abs(best["error_ppm"]) < 5
