#!/bin/bash
#SBATCH --job-name=sab_s3_a100
#SBATCH --partition=compute
#SBATCH --time=08:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH --account=sscm013902
#SBATCH --array=1-3
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

export ACCURACY=1
export RUN_ORDER=binary,ordinal

if [ "${SLURM_ARRAY_TASK_ID}" = "1" ]; then
  export RUN_TAG=chunk1
  export ORDINAL_SAMPLE_SIZES=500,750,1000,1250,1500
elif [ "${SLURM_ARRAY_TASK_ID}" = "2" ]; then
  export RUN_TAG=chunk2
  export ORDINAL_SAMPLE_SIZES=1750,2000,2500,3000,4000
else
  export RUN_TAG=chunk3
  export ORDINAL_SAMPLE_SIZES=5000,7500,10000,15000,20000
fi

Rscript code/06_run_ordinal_comparison.R
