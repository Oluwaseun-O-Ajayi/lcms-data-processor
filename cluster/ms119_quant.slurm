#!/bin/bash
# MS1 quantitation and fluorescence comparison for the nine MS119 "Joined" runs
# (trastuzumab, aflibercept, cetuximab; triplicate). Needs the mzML files made by cluster/ms119.slurm.
# Submit from a Sapelo2 login node:   sbatch ~/lcms-data-processor/cluster/ms119_quant.slurm
#SBATCH --job-name=ms119q
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-8
#SBATCH --output=/home/%u/jobs/logs/%x_%A_%a.out
set -eu
D=/scratch/$USER/glycan_data/pride
RUNS=($(cd $D/mzml_all && ls MS119_Joined_*DeglycPep_RFMS_Glycans_10a2i5uguL_0*.mzML | sed 's/\.mzML$//' | sort))
if [ ${#RUNS[@]} -ne 9 ]; then echo "Expected 9 MS119 mzML files, found ${#RUNS[@]}"; exit 1; fi
R=${RUNS[$SLURM_ARRAY_TASK_ID]}
case $R in
  *Trast*) T=cluster/ms119_trast_RFMS_targets.csv ;;
  *Afli*)  T=cluster/ms119_afli_RFMS_targets.csv ;;
  *Cetux*) T=cluster/ms119_cetux_RFMS_targets.csv ;;
esac
OUT=$D/results/ms119; mkdir -p $OUT
module load Python/3.11.3-GCCcore-12.3.0
source ~/venvs/lcms/bin/activate
cd ~/lcms-data-processor
echo "run: $R  targets: $T"
lcms-batch --mzml $D/mzml_all/$R.mzML --targets $T --peak-mode sum --ppm 10 --width 5 --baseline-window 0.1 --out $OUT/$R
# Fluorescence peak threshold: 0.005 (fraction of the region maximum) is the pre-set primary value for
# every dataset; 0.01 and 0.02 are run as a sensitivity check and written to sub-folders.
for P in 0.005 0.01 0.02; do
  if [ "$P" = "0.005" ]; then O=$OUT; else O=$OUT/fldp$P; mkdir -p $O; fi
  lcms-compare fld --mzml $D/mzml_all/$R.mzML --targets $OUT/${R}_targets.csv \
    --fld Emission --fld-prominence $P --ppm 10 --width 5 --baseline-window 0.1 --out $O/${R}_fld
done
echo "== DONE $R"
