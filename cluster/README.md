# Running the real-data study on a cluster (UGA GACRC Sapelo2)

These steps process the public released N-glycan dataset
(doi:10.5281/zenodo.4705546, CC BY 4.0) over a grid of extraction windows and
peak-detection settings. The same steps work for your own mzML files.

## 1. Download (transfer node)

```bash
ssh <MyID>@xfer.gacrc.uga.edu
cd /path/to/lcms-data-processor
bash cluster/download_glycan_data.sh /scratch/$USER/glycan_data
```

The script checks every file against the MD5 checksums published on Zenodo.

## 2. Convert Agilent .d to mzML

Agilent and Thermo raw formats need vendor libraries, which run on Windows.
The simplest route is ProteoWizard MSConvert on a Windows PC:

1. Copy the two `.d` folders to the PC (for example with WinSCP from `xfer.gacrc.uga.edu`).
2. In MSConvert, choose output format **mzML**, add the filter **Peak Picking**
   (algorithm: Vendor, MS levels 1-), and convert.
3. Copy the `.mzML` files back to `/scratch/$USER/glycan_data/mzml/`.

Record the MSConvert version and settings; they belong in the paper's methods.

## 3. Run the sweep

```bash
cd /path/to/lcms-data-processor
mkdir -p logs
sbatch cluster/sweep.slurm
squeue --me
```

Each array task writes three tables per file to `/scratch/$USER/glycan_data/results/`:

| File | Content |
|---|---|
| `<run>_targets.csv` | one row per setting and target m/z: retention time, area, height, FWHM, S/N, number of peaks |
| `<run>_glycans.csv` | areas summed over charge states and relative abundance (%) per glycan |
| `<run>_target_list.csv` | the exact m/z values used |

## 4. Check assignments before interpreting

The built-in compositions are candidates. Before using a glycan in the analysis,
confirm it from the data: a peak at a plausible HILIC retention time, consistent
charge states and isotope pattern, and MS/MS where available. Add confirmed
retention windows to a targets CSV (`name, mz, glycan, rt_min, rt_max`) and rerun
with `--targets` so every setting integrates the same peak.
