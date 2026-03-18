#!/bin/bash
#SBATCH --job-name=sab_paper_core
#SBATCH --partition=compute
#SBATCH --time=18:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
mkdir -p logs

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

echo "[$(date -Is)] Starting core paper pipeline (excluding nonPO and final assets)"

export PAPER_ACCURACY="${PAPER_ACCURACY:-1.00}"
export PAPER_BUILD_NONPO_PLOTS=0
export PAPER_SKIP_NONPO=1
export PAPER_SKIP_ASSETS=1
export NONPO_RUN_TAG="${NONPO_BASE_TAG:-paper_final_array_${SLURM_JOB_ID}}"

export SAVE_COHORT_RAW=1

Rscript code/00_run_paper_final.R

echo "[$(date -Is)] Finished core paper pipeline"
