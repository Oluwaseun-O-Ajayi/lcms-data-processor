# Analysis package for the JASMS application note

Every figure and every number in the manuscript is regenerated from the tables in `data/` by

```bash
python paper/jasms/make_figures.py
```

which writes `figures/` (Figures 1–3 and S1, PNG and PDF), `numbers.txt` (all quoted values),
`table1_agreement.csv` and `tableS_skyline_per_glycan.csv`.

## Data

| Folder | Content | Produced by |
|---|---|---|
| `data/fetuin_2AB/` | Settings sweep, fluorescence and Skyline comparisons, mass errors and ion-set comparison for the public 2-AB fetuin run (FET-2AB-CRL; doi:10.5281/zenodo.4705546) | `cluster/sweep.slurm`, `cluster/extras.slurm`, Skyline report, `examples/benchmark_windows.py` |
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
