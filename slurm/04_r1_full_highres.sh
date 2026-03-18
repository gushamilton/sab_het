#!/bin/bash
#SBATCH --job-name=sab_r1_highres
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
mkdir -p logs

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

echo "[$(date -Is)] Starting high-resolution R1 run"

# Main + ordinal use denser sample-size grid and higher Monte Carlo reps.
export SAMPLE_SIZES=500,750,1000,1250,1500,1750,2000,2500,3000,4000,5000,7500,10000,15000,20000
export ORDINAL_SAMPLE_SIZES=500,750,1000,1250,1500,1750,2000,2500,3000,4000,5000,7500,10000,15000,20000
export N_REPS_MAIN=10000
export N_REPS_COHORT=5000
export N_REPS_ORDINAL=10000
export COHORT_FIXED_N=5000
export SAVE_COHORT_RAW=1

echo "[$(date -Is)] Running main binary simulation"
Rscript code/03_run_main_binary.R

echo "[$(date -Is)] Running enrichment supplement"
Rscript code/04_run_enrichment.R

echo "[$(date -Is)] Running cohort-variation supplement (with raw output)"
Rscript code/05_run_cohort_variation.R

echo "[$(date -Is)] Running ordinal supplement (80% accuracy)"
export ACCURACY=0.80
export RUN_ORDER=binary,ordinal
Rscript code/06_run_ordinal_comparison.R

echo "[$(date -Is)] Running ordinal supplement (100% accuracy sensitivity run)"
export ACCURACY=1.00
Rscript code/06_run_ordinal_comparison.R

echo "[$(date -Is)] Building GT summary tables"
Rscript code/07_build_gt_tables.R

echo "[$(date -Is)] Completed high-resolution R1 run"
