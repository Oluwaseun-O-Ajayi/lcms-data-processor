#!/bin/bash
#SBATCH --job-name=procA-replicates
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#
# Targeted quantitation of the procainamide runs that passed QC:
#   FET-ProA-S1 and FET-ProA-S2 (two fetuin samples) and TFN-ProA (human transferrin).
# Also writes the run manifest and QC table for every procainamide run surveyed.
# Submit from the repository folder:   mkdir -p logs && sbatch cluster/procA_replicates.slurm
set -euo pipefail

MZML=/scratch/$USER/glycan_data/procA/mzml
EXPLORE=/scratch/$USER/glycan_data/procA/explore
OUT=/scratch/$USER/glycan_data/procA/final
PPMS="10 20 30 50"

module load Python/3.11.3-GCCcore-12.3.0
source ~/venvs/lcms/bin/activate
mkdir -p "$OUT"

pick() { ls "$MZML"/*"$1"*.mzML | head -1; }
FET1=$(pick Fet1_052124); FET2=$(pick Fet2_052124); TFN=$(pick HTferin052124)
echo "FET-ProA-S1: $FET1"; echo "FET-ProA-S2: $FET2"; echo "TFN-ProA:    $TFN"

python cluster/make_procA_targets.py

lcms-batch --mzml "$FET1" --targets cluster/procA_fetuin_targets.csv --peak-mode sum \
  --ppm $PPMS --width 5 --baseline-window 0.1 --out "$OUT/FET-ProA-S1"
lcms-batch --mzml "$FET2" --targets cluster/procA_fetuin_targets.csv --peak-mode sum \
  --ppm $PPMS --width 5 --baseline-window 0.1 --out "$OUT/FET-ProA-S2"
lcms-batch --mzml "$TFN" --targets cluster/procA_transferrin_targets.csv --peak-mode sum \
  --ppm $PPMS --width 5 --baseline-window 0.1 --out "$OUT/TFN-ProA"

python examples/mass_errors.py --mzml "$FET1" --targets cluster/procA_fetuin_targets.csv --out "$OUT/FET-ProA-S1_mass_errors.csv" > /dev/null
python examples/mass_errors.py --mzml "$FET2" --targets cluster/procA_fetuin_targets.csv --out "$OUT/FET-ProA-S2_mass_errors.csv" > /dev/null
python examples/mass_errors.py --mzml "$TFN"  --targets cluster/procA_transferrin_targets.csv --out "$OUT/TFN-ProA_mass_errors.csv" > /dev/null

python cluster/procA_summary.py --results "$OUT" --out "$OUT/procA"
python cluster/procA_qc.py --hits-dir "$EXPLORE" --out "$OUT/procA"
echo "Finished. Results in $OUT"
