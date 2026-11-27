# Analysis package for the GlyVAS validation paper

Every figure and every number in the manuscript is regenerated from the tables in `data/` by

```bash
python paper/glycan-validation/make_figures.py
```

which writes `figures/` (Figures 1–3 and S1–S4 as PNG and PDF, and the TOC graphic as PNG, PDF and TIF), `numbers.txt`
(all quoted values), `table1_agreement.csv` and `tableS_skyline_per_glycan.csv`.
Glycan structures in Figure 2 and in Table S1 of the Supporting Information are drawn by `snfg.py` in the Symbol
Nomenclature for Glycans (SNFG).

## Data

| Folder | Content | Produced by |
|---|---|---|
| `data/fetuin_2AB/` | Settings sweep, fluorescence and Skyline comparisons, mass errors and ion-set comparison for the public 2-AB fetuin run (FET-2AB-CRL; doi:10.5281/zenodo.4705546) | `cluster/sweep.slurm`, `cluster/extras.slurm`, Skyline report, `examples/benchmark_windows.py` |
| `data/trastuzumab_RFMS/` | MS154 trastuzumab (PXD078897; 9 runs, 3 days): agreement per run and per glycan, variance components, cross-day calibration, settings sensitivity | `cluster/trast_all.slurm`, `combine_replicates.py` |
| `data/ms119_per_glycan/` | MS119 trastuzumab, aflibercept and cetuximab: MS and fluorescence areas per glycan and per run, fluorescence region summaries | `cluster/ms119_quant.slurm` (`lcms-compare fld`) |
| `data/ms119_transfer/` | MS119 full-profile agreement and coverage, calibration within protein, transfer of MS154 response factors | `transfer_test.py` |
| `data/sensitivity/` | The same analyses at fluorescence peak thresholds of 0.01 and 0.02 | `cluster/rerun_quant.sh` |
| `data/diag/` | Share of the MS1 signal taken by each glycan in its window, top other ions, G1F isomer-pair traces | `diag_windows.py`, `xic_pair.py` (`cluster/diag_windows.slurm`, `cluster/xic_pairs.slurm`) |
| `data/converter_check/` | ThermoRawFileParser vs MSConvert on one run | `converter_check.py` (`cluster/converter_check.slurm`) |
| `data/procainamide/` | Run manifest, label-chemistry check, relative abundances at ±10 and ±20 ppm, mass errors and dextran ladder for the procainamide runs | `cluster/procA_explore.slurm`, `cluster/procA_replicates.slurm` |

## Sample names

| Name | Sample |
|---|---|
| FET-2AB-CRL | Public 2-AB-labeled bovine fetuin run |
| FET-ProA-S1, FET-ProA-S2 | Procainamide-labeled bovine fetuin, samples 1 and 2 |
| TFN-ProA | Procainamide-labeled human transferrin |
| DEX-ProA | Procainamide-labeled dextran ladder, same sequence |
| FET-ProA-ACN-I1, -I2 | Fetuin stored in 50:50 acetonitrile:water and frozen (excluded) |
| DEX-ProA-ACN | Dextran ladder from the same sequence (excluded) |
| FET-ProA-X1 to X6, TFN-ProA-X1 | Runs that failed the label-chemistry check (excluded) |

Original raw file names are in `data/procainamide/manifest.csv`.
