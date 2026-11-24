#!/bin/bash
# Re-run all quantitation and both summaries, in order, as one chain of jobs.
# Run from a Sapelo2 login node:   bash ~/lcms-data-processor/cluster/rerun_quant.sh
# Each step starts only if the step it depends on finished without error.
set -eu
cd ~/lcms-data-processor
for f in cluster/trast_RFMS_targets.csv cluster/ms119_trast_RFMS_targets.csv cluster/ms119_afli_RFMS_targets.csv \
         cluster/ms119_cetux_RFMS_targets.csv cluster/trast_all.slurm cluster/ms119_quant.slurm \
         paper/glycan-validation/combine_replicates.py paper/glycan-validation/transfer_test.py \
         paper/glycan-validation/diag_windows.py cluster/diag_windows.slurm \
         paper/glycan-validation/xic_pair.py cluster/xic_pairs.slurm tests/test_fld_detection.py; do
  [ -s "$f" ] || { echo "MISSING $f (run git pull first)"; exit 1; }
done
grep -q "^G0_1+," cluster/trast_RFMS_targets.csv || { echo "Old targets file: run git pull first"; exit 1; }
grep -q "fldp" cluster/trast_all.slurm || { echo "Old trast_all.slurm: run git pull first"; exit 1; }
grep -q "height_fraction=fld_prominence" reference_compare.py || { echo "Old reference_compare.py: run git pull first"; exit 1; }
grep -q "fld_region_area" paper/glycan-validation/combine_replicates.py || { echo "Old combine_replicates.py: run git pull first"; exit 1; }
mkdir -p ~/jobs/logs
ENV="module load Python/3.11.3-GCCcore-12.3.0 && source ~/venvs/lcms/bin/activate && cd ~/lcms-data-processor"
R=/scratch/$USER/glycan_data/pride/results
J1=$(sbatch --parsable cluster/trast_all.slurm)
J2=$(sbatch --parsable cluster/ms119_quant.slurm)
C="python paper/glycan-validation/combine_replicates.py"
T="python paper/glycan-validation/transfer_test.py"
G=paper/glycan-validation/data
J3=$(sbatch --parsable --dependency=afterok:$J1 --job-name=combine --partition=batch --ntasks=1 --mem=4G --time=00:30:00 \
     --output=$HOME/jobs/logs/combine_%j.out \
     --wrap="$ENV && $C --results $R/trast --out $G/trastuzumab_RFMS && $C --results $R/trast/fldp0.01 --out $G/sensitivity/trast_fldp0.01 && $C --results $R/trast/fldp0.02 --out $G/sensitivity/trast_fldp0.02")
J4=$(sbatch --parsable --dependency=afterok:$J2:$J3 --job-name=transfer --partition=batch --ntasks=1 --mem=4G --time=00:30:00 \
     --output=$HOME/jobs/logs/transfer_%j.out \
     --wrap="$ENV && $T --ref $G/trastuzumab_RFMS/per_run_per_glycan.csv --results $R/ms119 --out $G/ms119_transfer && $T --ref $G/trastuzumab_RFMS/per_run_per_glycan.csv --results $R/ms119/fldp0.01 --out $G/sensitivity/ms119_fldp0.01 && $T --ref $G/trastuzumab_RFMS/per_run_per_glycan.csv --results $R/ms119/fldp0.02 --out $G/sensitivity/ms119_fldp0.02")
J5=$(sbatch --parsable --dependency=afterok:$J1:$J2 cluster/diag_windows.slurm)
J6=$(sbatch --parsable --dependency=afterok:$J1:$J2 cluster/xic_pairs.slurm)
echo "Submitted: trast-all $J1, ms119q $J2, combine $J3, transfer $J4, diag $J5, xicpair $J6"
echo "When 'sq --me' is empty, collect everything into one file with:"
echo "  cat $HOME/jobs/logs/combine_$J3.out $HOME/jobs/logs/transfer_$J4.out $HOME/jobs/logs/diag_$J5.out $HOME/jobs/logs/xicpair_$J6.out > ~/results_all.txt"
