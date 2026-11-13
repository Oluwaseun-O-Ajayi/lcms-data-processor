#!/bin/bash
#SBATCH --job-name=procA-explore
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.out
#
# First look at the procainamide-labelled LC-IM-MS runs (ion mobility combined in MSConvert).
#   1. TIC of every run
#   2. dextran ladder (GU2-GU20, 1+ and 2+): apex retention time and mass error per GU
#   3. composition survey of every fetuin run (2+, 3+, 4+), to locate the glycans
# Submit from the repository folder:   sbatch cluster/procA_explore.slurm
set -euo pipefail

IN=/scratch/$USER/glycan_data/procA/mzml
OUT=/scratch/$USER/glycan_data/procA/explore

module load Python/3.11.3-GCCcore-12.3.0
source ~/venvs/lcms/bin/activate
mkdir -p "$OUT"

shopt -s nullglob
files=("$IN"/*.mzML)
echo "Found ${#files[@]} mzML files"
[ ${#files[@]} -gt 0 ] || { echo "No mzML files in $IN"; exit 1; }

for f in "${files[@]}"; do
  base=$(basename "$f" .mzML)
  short=$(echo "$base" | sed -E 's/\[result\]//; s/^OA_([0-9_]+)[-_][0-9]+\.[0-9]+\.[0-9]+-(.*)$/\2_\1/; s/^[0-9.-]+-//')
  if [ -f "$OUT/${short}_tic.csv" ]; then echo "== $short (done before, skipped)"; continue; fi
  echo "== $short"
  python - "$f" "$OUT/${short}_tic.csv" <<'PY'
import sys
from lcms_data_processor import MzMLReader
MzMLReader(sys.argv[1]).tic().to_csv(sys.argv[2], index=False)
PY
  case "$short" in
    *Dex*)
      python examples/mass_errors.py --mzml "$f" --targets cluster/procA_dextran_targets.csv \
        --out "$OUT/${short}_dextran_ladder.csv" ;;
    *Fet*|*ferin*|*Tfr*|*Trnsf*)
      lcms-glycan-search --mzml "$f" --survey --label procainamide --charges 2 3 4 --ppm 20 \
        --out "$OUT/${short}_survey.csv" > /dev/null ;;
  esac
done
echo "Finished. Results in $OUT"
ls -lh "$OUT"
