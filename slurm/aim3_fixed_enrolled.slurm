#!/bin/bash
#SBATCH --job-name=aim3_fixed_enrolled
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
mkdir -p logs

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

SCENARIO_SET=main \
AIM3_TARGET_ENROLLED=5000 \
N_REPS_AIM3_FIXED=500 \
Rscript code/05b_run_aim3_fixed_enrolled.R
