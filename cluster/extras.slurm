#!/bin/bash
#SBATCH --job-name=lcms-extras
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.out
#
# Follow-up analyses on the fetuin run:
#   1. ion coverage: 2+/3+ protonated (main study) vs adding 4+ vs adding Na adducts
#   2. measured mass error of every target ion
# Submit from the repository folder:   sbatch cluster/extras.slurm
set -euo pipefail

MZML=/scratch/$USER/glycan_data/mzml/MASOP00220_AB_CRL_07JAN2021_Fetuin_001.mzML
OUT=/scratch/$USER/glycan_data/results/extras
PPMS="5 10 20 50"

module load Python/3.11.3-GCCcore-12.3.0
source ~/venvs/lcms/bin/activate
mkdir -p "$OUT"

declare -A TARGETS=([H23]=cluster/fetuin_2AB_targets.csv
                    [H234]=cluster/fetuin_2AB_targets_2to4H.csv
                    [H234Na]=cluster/fetuin_2AB_targets_with_Na.csv)
for set in H23 H234 H234Na; do
  T=${TARGETS[$set]}
  lcms-batch --mzml "$MZML" --targets "$T" --peak-mode sum --ppm $PPMS --width 5 \
    --baseline-window 0.1 --out "$OUT/fetuin_$set"
  for p in $PPMS; do
    lcms-compare fld --mzml "$MZML" --targets "$OUT/fetuin_${set}_targets.csv" --fld FLD1 \
      --ppm "$p" --width 5 --baseline-window 0.1 --out "$OUT/fld_${set}_${p}ppm" > /dev/null
  done
done

python examples/mass_errors.py --mzml "$MZML" --targets cluster/fetuin_2AB_targets_with_Na.csv \
  --out "$OUT/mass_errors.csv"

python - "$OUT" <<'PY'
import glob, sys, pandas as pd
out = sys.argv[1]
rows = []
for f in sorted(glob.glob(f"{out}/fld_*_summary.csv")):
    s = pd.read_csv(f)
    s['ions'] = f.split('/fld_')[1].split('_')[0]
    rows.append(s)
s = pd.concat(rows)
print(s[['ions', 'ppm', 'n', 'mean_abs_diff_pct_points', 'pearson_r']].sort_values(['ppm', 'ions']).to_string(index=False))
s.to_csv(f"{out}/adduct_fld_agreement.csv", index=False)
PY
echo "Finished. Results in $OUT"
