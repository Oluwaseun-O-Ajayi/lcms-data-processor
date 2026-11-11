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

## 5. Reference comparisons

### Fluorescence (2-AB and other fluorescent labels)

Fluorescence area is proportional to molar amount, so it is the accepted reference for
labelled-glycan quantitation. If the run recorded an FLD trace, MSConvert exports it as a
chromatogram in the mzML file. Check what the file contains:

```bash
lcms-compare fld --mzml run.mzML --targets results/run_targets.csv --list
```

Then compare at any setting from the sweep (repeat for several settings to show which
ones agree with fluorescence):

```bash
lcms-compare fld --mzml run.mzML --targets results/run_targets.csv --fld FLD \
                 --ppm 10 --width 5 --baseline-window 0.1 --out results/run_fld_10ppm
```

The FLD-to-MS transfer delay is estimated from the data (cross-correlation with the
MS total-ion chromatogram) unless you pass `--delay`. Glycans that land on the same
fluorescence peak cannot be told apart by fluorescence; they are marked
`shared FLD peak` and left out of the agreement statistics.

If the mzML has no FLD chromatogram, check MSConvert's options for exporting
detector traces, or export the FLD trace as CSV from the vendor software.

### Skyline

1. Write a transition list from the same m/z values:
   `lcms-compare skyline-targets --targets results/run_target_list.csv --out skyline_targets.csv`
2. In Skyline, set Full-Scan MS1 filtering to match this package: monoisotopic peak
   only, centroided data, and the same mass accuracy (ppm) as `--ppm`.
3. Import the list (Edit > Insert > Transition List), import the mzML results, and
   check the peak boundaries.
4. Export a report with the columns `Molecule Name`, `Precursor Charge`,
   `Replicate Name` and `Total Area MS1` (File > Export > Report).
5. Compare:
   `lcms-compare skyline --targets results/run_targets.csv --skyline report.csv --ppm 10 --out results/run_skyline`

Skyline menus and setting names change between versions; check them against the
version you use and record it in the methods.
