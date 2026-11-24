#!/bin/bash
# Re-run all quantitation and both summaries, in order, as one chain of jobs.
# Run from a Sapelo2 login node:   bash ~/lcms-data-processor/cluster/rerun_quant.sh
# Each step starts only if the step it depends on finished without error.
set -eu
cd ~/lcms-data-processor
for f in cluster/trast_RFMS_targets.csv cluster/ms119_trast_RFMS_targets.csv cluster/ms119_afli_RFMS_targets.csv \
         cluster/ms119_cetux_RFMS_targets.csv cluster/trast_all.slurm cluster/ms119_quant.slurm \
         paper/glycan-validation/combine_replicates.py paper/glycan-validation/transfer_test.py; do
  [ -s "$f" ] || { echo "MISSING $f (run git pull first)"; exit 1; }
done
grep -q "^G0_1+," cluster/trast_RFMS_targets.csv || { echo "Old targets file: run git pull first"; exit 1; }
mkdir -p ~/jobs/logs
ENV="module load Python/3.11.3-GCCcore-12.3.0 && source ~/venvs/lcms/bin/activate && cd ~/lcms-data-processor"
R=/scratch/$USER/glycan_data/pride/results
J1=$(sbatch --parsable cluster/trast_all.slurm)
J2=$(sbatch --parsable cluster/ms119_quant.slurm)
J3=$(sbatch --parsable --dependency=afterok:$J1 --job-name=combine --partition=batch --ntasks=1 --mem=4G --time=00:15:00 \
     --output=$HOME/jobs/logs/combine_%j.out \
     --wrap="$ENV && python paper/glycan-validation/combine_replicates.py --results $R/trast --out paper/glycan-validation/data/trastuzumab_RFMS")
J4=$(sbatch --parsable --dependency=afterok:$J2:$J3 --job-name=transfer --partition=batch --ntasks=1 --mem=4G --time=00:15:00 \
     --output=$HOME/jobs/logs/transfer_%j.out \
     --wrap="$ENV && python paper/glycan-validation/transfer_test.py --ref paper/glycan-validation/data/trastuzumab_RFMS/per_run_per_glycan.csv --results $R/ms119 --out paper/glycan-validation/data/ms119_transfer")
echo "Submitted: trast-all $J1, ms119q $J2, combine $J3 (after $J1), transfer $J4 (after $J2 and $J3)"
echo "When 'sq --me' is empty, read the results with:"
echo "  cat $HOME/jobs/logs/combine_$J3.out $HOME/jobs/logs/transfer_$J4.out"
