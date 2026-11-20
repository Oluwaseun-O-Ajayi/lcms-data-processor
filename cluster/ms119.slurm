#!/bin/bash
# Convert the nine MS119 "Joined" runs (trastuzumab, aflibercept, cetuximab; triplicate) and list the
# strongest MS1 ions at each fluorescence apex in the 48-80 min glycan window.
# Submit from a Sapelo2 login node after all nine .raw files are downloaded:
#   sbatch ~/lcms-data-processor/cluster/ms119.slurm
#SBATCH --job-name=ms119
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --array=0-8
#SBATCH --output=/home/%u/jobs/logs/%x_%A_%a.out
set -eu
D=/scratch/$USER/glycan_data/pride
RUNS=($(cd $D && ls MS119_Joined_*DeglycPep_RFMS_Glycans_10a2i5uguL_0*.raw | sed 's/\.raw$//' | sort))
if [ ${#RUNS[@]} -ne 9 ]; then echo "Expected 9 MS119 raw files, found ${#RUNS[@]}"; exit 1; fi
R=${RUNS[$SLURM_ARRAY_TASK_ID]}
mkdir -p $D/mzml_all $D/fld_ms119
if [ ! -s $D/mzml_all/$R.mzML ]; then
  ldconfig -p | grep -q libicu || export DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1 DOTNET_SYSTEM_GLOBALIZATION_PREDEFINED_CULTURES_ONLY=false
  ~/tools/trfp/ThermoRawFileParser -i $D/$R.raw -a -f 2 -o $D/mzml_all -l 4
fi
module load Python/3.11.3-GCCcore-12.3.0
source ~/venvs/lcms/bin/activate
python ~/lcms-data-processor/cluster/fld_overlay.py $D/mzml_all/$R.mzML --channel Emission --window 48 80 --spectra --out $D/fld_ms119
echo "== DONE $R"
