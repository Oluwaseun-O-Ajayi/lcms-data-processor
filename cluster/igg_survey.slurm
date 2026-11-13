#!/bin/bash
#SBATCH --job-name=igg-survey
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#
# First look at the public 2-AB human IgG run (same Zenodo data set as the fetuin run):
# lists the stored chromatograms, surveys glycan compositions over the whole run, and lists
# the fluorescence peaks. The results are used to write the IgG target list with retention
# windows, which the follow-up job (sweep + fluorescence comparison) then uses.
# Submit from the repository folder:   mkdir -p logs && sbatch cluster/igg_survey.slurm
set -euo pipefail

MZML=/scratch/$USER/glycan_data/mzml/MASOP00011_AB_0332_16DEC2020_hIgG_NoDE_40ul_001.mzML
OUT=/scratch/$USER/glycan_data/results/igg
mkdir -p "$OUT"
[ -f "$MZML" ] || { echo "Missing $MZML: convert and upload the IgG run first"; exit 1; }

module load Python/3.11.3-GCCcore-12.3.0
source ~/venvs/lcms/bin/activate

lcms-compare fld --mzml "$MZML" --targets cluster/fetuin_2AB_targets.csv --list | tee "$OUT/igg_chromatograms.txt"

lcms-glycan-search --mzml "$MZML" --survey --label 2-AB --charges 1 2 3 --ppm 10 \
  --fld FLD1 --fld-start 5 --out "$OUT/igg_survey.csv" > "$OUT/igg_survey_log.txt"

python - "$MZML" "$OUT/igg_tic.csv" <<'PY'
import sys
from lcms_data_processor import MzMLReader
MzMLReader(sys.argv[1]).tic().to_csv(sys.argv[2], index=False)
PY
ls -l "$OUT"
echo "Finished. Results in $OUT"
