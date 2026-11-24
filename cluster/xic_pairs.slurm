#!/bin/bash
# G1F isomer pairs: MS XIC vs fluorescence on each side of the window boundary (afli, trast MS119, trast MS154).
#SBATCH --job-name=xicpair
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/home/%u/jobs/logs/%x_%j.out
set -eu
D=/scratch/$USER/glycan_data/pride/mzml_all
O=~/lcms-data-processor/paper/glycan-validation/data/diag; mkdir -p $O
module load Python/3.11.3-GCCcore-12.3.0
source ~/venvs/lcms/bin/activate
cd ~/lcms-data-processor
X="python paper/glycan-validation/xic_pair.py --mz 1936.7791 968.8932 646.2646 --z 1 2 3"
echo "## aflibercept MS119 run 01"
$X --mzml $D/MS119_Joined_AfliDeglycPep_RFMS_Glycans_10a2i5uguL_01.mzML --start 51.5 --end 54.5 --split 52.96 --out $O/pair_afli01.png
echo "## trastuzumab MS119 run 01"
$X --mzml $D/MS119_Joined_TrastDeglycPep_RFMS_Glycans_10a2i5uguL_01.mzML --start 63.3 --end 66.7 --split 65.00 --out $O/pair_trast119_01.png
echo "## trastuzumab MS154 day 1 run 01"
$X --mzml $D/MS154_Trast_1stDay_HILIC_01.mzML --start 59.8 --end 63.1 --split 61.44 --out $O/pair_trast154_01.png
echo "== DONE xicpair"
