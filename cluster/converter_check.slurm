#!/bin/bash
# Converter cross-check (MS154 day 1 run 01): ThermoRawFileParser mzML vs MSConvert mzML,
# same targets and settings.  Submit from Sapelo2:  sbatch cluster/converter_check.slurm
#SBATCH --job-name=convcheck
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/%u/jobs/logs/%x_%j.out
set -eu
R=MS154_Trast_1stDay_HILIC_01
D=/scratch/$USER/glycan_data/pride
A=$D/mzml_all/$R.mzML
B=$D/msconvert_check.mzML
FLD=$D/results/trast/${R}_fld_per_glycan.csv
OUT=$D/results/convcheck; mkdir -p $OUT
for f in $A $B $FLD; do [ -s $f ] || { echo "MISSING: $f"; exit 1; }; done

echo "== File checks"
for f in $A $B; do
  echo "$f"
  echo "  size: $(du -h $f | cut -f1)"
  echo "  $(grep -m1 -o '<spectrumList count="[0-9]*"' $f)"
  echo "  centroid spectra (MS:1000127): $(grep -c 'MS:1000127' $f)"
  echo "  profile spectra (MS:1000128): $(grep -c 'MS:1000128' $f)"
done

module load Python/3.11.3-GCCcore-12.3.0
source ~/venvs/lcms/bin/activate
cd ~/lcms-data-processor
export PYTHONPATH=$HOME/lcms-data-processor${PYTHONPATH:+:$PYTHONPATH}
for pair in "trfp:$A" "msconvert:$B"; do
  lcms-batch --mzml ${pair#*:} --targets cluster/trast_RFMS_targets.csv \
    --peak-mode sum --ppm 5 10 20 --width 3 5 --baseline-window 0.05 0.1 --out $OUT/${pair%%:*}
done
python paper/glycan-validation/converter_check.py --trfp $OUT/trfp_glycans.csv \
  --msconvert $OUT/msconvert_glycans.csv --fld $FLD \
  --out paper/glycan-validation/data/converter_check/convcheck
echo "== DONE"
