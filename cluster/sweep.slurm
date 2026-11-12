#!/bin/bash
#SBATCH --job-name=lcms-sweep
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#
# Settings sweep on the fetuin 2-AB run, then the fluorescence comparison at every setting.
# Submit from the repository folder:   mkdir -p logs && sbatch cluster/sweep.slurm
#
# One-time setup (on a compute node from `interact`):
#   module load Python/3.11.3-GCCcore-12.3.0
#   python3 -m venv ~/venvs/lcms && source ~/venvs/lcms/bin/activate
#   pip install -e ".[mzml]"
set -euo pipefail

MZML=/scratch/$USER/glycan_data/mzml/MASOP00220_AB_CRL_07JAN2021_Fetuin_001.mzML
TARGETS=cluster/fetuin_2AB_targets.csv
OUT=/scratch/$USER/glycan_data/results
PPMS="2 5 10 20 50 100"
WIDTHS="3 5 8"
BWS="0.02 0.05 0.1 0.2"

module load Python/3.11.3-GCCcore-12.3.0
source ~/venvs/lcms/bin/activate
mkdir -p "$OUT/fld"

lcms-batch --mzml "$MZML" --targets "$TARGETS" --peak-mode sum \
  --ppm $PPMS --width $WIDTHS --baseline-window $BWS --out "$OUT/fetuin"

python examples/summarize_sweep.py "$OUT/fetuin_glycans.csv" --out "$OUT/fetuin_summary"

for p in $PPMS; do for w in $WIDTHS; do for b in $BWS; do
  lcms-compare fld --mzml "$MZML" --targets "$OUT/fetuin_targets.csv" --fld FLD1 \
    --ppm "$p" --width "$w" --baseline-window "$b" --out "$OUT/fld/fetuin_${p}ppm_w${w}_b${b}" > /dev/null
done; done; done

python - "$OUT" <<'PY'
import glob, sys, pandas as pd
out = sys.argv[1]
s = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{out}/fld/*_summary.csv"))], ignore_index=True)
s.to_csv(f"{out}/fetuin_fld_agreement_all_settings.csv", index=False)
cols = ['ppm', 'width', 'baseline_window', 'n', 'mean_abs_diff_pct_points', 'bias_pct_points', 'pearson_r', 'delay_min']
print(s[cols].sort_values('mean_abs_diff_pct_points').to_string(index=False))
PY
echo "Finished. Results in $OUT"
