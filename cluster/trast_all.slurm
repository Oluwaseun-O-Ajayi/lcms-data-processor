#!/bin/bash
# MS154 trastuzumab: 9 runs (3 days x triplicate). Converts .raw if needed, MS1 quantitation over the
# settings grid, and fluorescence comparison.  Submit from Sapelo2:  sbatch cluster/trast_all.slurm
#SBATCH --job-name=trast-all
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-8
#SBATCH --output=/home/%u/jobs/logs/%x_%A_%a.out
set -eu
RUNS=(); for d in 1st 2nd 3rd; do for r in 01 02 03; do RUNS+=("MS154_Trast_${d}Day_HILIC_${r}"); done; done
R=${RUNS[$SLURM_ARRAY_TASK_ID]}
D=/scratch/$USER/glycan_data/pride
OUT=$D/results/trast; mkdir -p $OUT $D/mzml_all
if [ ! -s $D/mzml_all/$R.mzML ]; then
  ldconfig -p | grep -q libicu || export DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1 DOTNET_SYSTEM_GLOBALIZATION_PREDEFINED_CULTURES_ONLY=false
  ~/tools/trfp/ThermoRawFileParser -i $D/$R.raw -a -f 2 -o $D/mzml_all -l 4
fi
module load Python/3.11.3-GCCcore-12.3.0
source ~/venvs/lcms/bin/activate
cd ~/lcms-data-processor
# use the code in this repository, not an older copy installed in the venv
export PYTHONPATH=$HOME/lcms-data-processor${PYTHONPATH:+:$PYTHONPATH}
python -c "import reference_compare, lcms_data_processor as l, inspect; print('code from:', reference_compare.__file__); assert 'height_fraction' in inspect.getsource(l.PeakDetector.detect_peaks)"
lcms-batch --mzml $D/mzml_all/$R.mzML --targets cluster/trast_RFMS_targets.csv \
  --peak-mode sum --ppm 5 10 20 --width 3 5 --baseline-window 0.05 0.1 --out $OUT/$R
# Fluorescence peak threshold: 0.005 (fraction of the region maximum) is the pre-set primary value for
# every dataset; 0.01 and 0.02 are run as a sensitivity check and written to sub-folders.
for P in 0.005 0.01 0.02; do
  if [ "$P" = "0.005" ]; then O=$OUT; else O=$OUT/fldp$P; mkdir -p $O; fi
  lcms-compare fld --mzml $D/mzml_all/$R.mzML --targets $OUT/${R}_targets.csv \
    --fld Emission --fld-prominence $P --ppm 10 --width 5 --baseline-window 0.1 --out $O/${R}_fld
done
echo "== DONE $R"
