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
#   1. sodium adducts: protonated ions only versus protonated plus [M+H+Na] and [M+2H+Na]
#   2. measured mass error of every target ion
# Submit from the repository folder:   sbatch cluster/extras.slurm
set -euo pipefail

MZML=/scratch/$USER/glycan_data/mzml/MASOP00220_AB_CRL_07JAN2021_Fetuin_001.mzML
OUT=/scratch/$USER/glycan_data/results/extras
PPMS="5 10 20 50"

module load Python/3.11.3-GCCcore-12.3.0
source ~/venvs/lcms/bin/activate
mkdir -p "$OUT"

for set in H Na; do
  T=cluster/fetuin_2AB_targets.csv
  [[ $set == Na ]] && T=cluster/fetuin_2AB_targets_with_Na.csv
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
    s['ions'] = 'protonated only' if '/fld_H_' in f else 'protonated + Na adducts'
    rows.append(s)
s = pd.concat(rows)
print(s[['ions', 'ppm', 'n', 'mean_abs_diff_pct_points', 'pearson_r']].sort_values(['ppm', 'ions']).to_string(index=False))
s.to_csv(f"{out}/adduct_fld_agreement.csv", index=False)
PY
echo "Finished. Results in $OUT"
