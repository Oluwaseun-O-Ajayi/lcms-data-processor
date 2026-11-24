#!/bin/bash
# Window diagnostics (glycan share of MS1 signal, dominant foreign ions) for every dataset.
#SBATCH --job-name=diag
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/home/%u/jobs/logs/%x_%j.out
set -eu
D=/scratch/$USER/glycan_data/pride/mzml_all
O=~/lcms-data-processor/paper/glycan-validation/data/diag; mkdir -p $O
module load Python/3.11.3-GCCcore-12.3.0
source ~/venvs/lcms/bin/activate
cd ~/lcms-data-processor
S=paper/glycan-validation/diag_windows.py
python $S --targets cluster/trast_RFMS_targets.csv       --mzml $D/MS154_Trast_*Day_HILIC_0*.mzML --out $O/ms154_trast.csv
python $S --targets cluster/ms119_trast_RFMS_targets.csv --mzml $D/MS119_Joined_Trast*.mzML      --out $O/ms119_trast.csv
python $S --targets cluster/ms119_afli_RFMS_targets.csv  --mzml $D/MS119_Joined_Afli*.mzML       --out $O/ms119_afli.csv
python $S --targets cluster/ms119_cetux_RFMS_targets.csv --mzml $D/MS119_Joined_Cetux*.mzML      --out $O/ms119_cetux.csv
echo "== DONE diag"
