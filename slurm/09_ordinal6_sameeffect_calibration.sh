#!/bin/bash
#SBATCH --job-name=sab_ord6_cal
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
mkdir -p logs

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

export SAMPLE_SIZES=500
export N_REPS_ORD6=1
export RUN_TAG=calibration
export SKIP_CALIBRATION=0
export CALIBRATION_N=200000

Rscript code/14_run_ordinal6_same_effect_comparison.R
