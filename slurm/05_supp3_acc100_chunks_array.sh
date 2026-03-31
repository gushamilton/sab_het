#!/bin/bash
#SBATCH --job-name=sab_s3_a100
#SBATCH --partition=short,compute
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=sscm013902
#SBATCH --array=1-15
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

export ACCURACY=1
export ORDINAL_POINTS=6
export RUN_ORDER=binary,ordinal

sample_sizes=(500 750 1000 1250 1500 1750 2000 2500 3000 4000 5000 7500 10000 15000 20000)
idx=$((SLURM_ARRAY_TASK_ID - 1))
sample_size="${sample_sizes[$idx]}"

export RUN_TAG="chunk${SLURM_ARRAY_TASK_ID}"
export ORDINAL_SAMPLE_SIZES="${sample_size}"

Rscript code/06_run_ordinal_comparison.R
