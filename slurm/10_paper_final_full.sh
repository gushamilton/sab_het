#!/bin/bash
#SBATCH --job-name=sab_paper_final
#SBATCH --partition=compute
#SBATCH --time=36:00:00
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

echo "[$(date -Is)] Starting final paper pipeline"

# Final manuscript defaults.
export PAPER_ACCURACY=1.00
export PAPER_BUILD_NONPO_PLOTS=1
export NONPO_RUN_TAG="paper_final_${SLURM_JOB_ID}"

# Keep all raw outputs for secondary analyses.
export SAVE_COHORT_RAW=1

# Use script defaults for full run unless overridden here:
# export N_REPS_MAIN=2000
# export N_REPS_COHORT=1000
# export N_REPS_ORDINAL=1000
# export N_REPS_NONPO=1000

Rscript code/00_run_paper_final.R

echo "[$(date -Is)] Finished final paper pipeline"
