#!/bin/bash
#SBATCH --job-name=lcms-sweep
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/%x_%A_%a.out
#
# One array task per mzML file; each task sweeps the full settings grid.
# Submit from the repository folder:   mkdir -p logs && sbatch cluster/sweep.slurm
# Set the array range to 0-(number of files - 1).
#
# One-time setup (on a login node):
#   module load Python          # check available versions with: module spider Python
#   python -m venv ~/venvs/lcms && source ~/venvs/lcms/bin/activate
#   pip install -e ".[mzml]"
set -euo pipefail

DATA=/scratch/$USER/glycan_data/mzml
OUT=/scratch/$USER/glycan_data/results

module load Python
source ~/venvs/lcms/bin/activate
mkdir -p "$OUT"

FILES=("$DATA"/*.mzML)
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
NAME=$(basename "$FILE" .mzML)
PRESET=fetuin
[[ "$NAME" == *IgG* ]] && PRESET=igg

lcms-batch --mzml "$FILE" --glycans "$PRESET" --label 2-AB --charges 1 2 3 \
  --ppm 2 5 10 20 50 100 \
  --width 3 5 8 \
  --baseline-window 0.02 0.05 0.1 0.2 \
  --out "$OUT/$NAME"
