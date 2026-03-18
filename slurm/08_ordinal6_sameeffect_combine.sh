#!/bin/bash
#SBATCH --job-name=sab_ord6_merge
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

export CHUNK_SAMPLE_SIZES=500,1000,1500,2000,3000,5000,10000,20000
export N_REPS_TOTAL=1000
export REPS_PER_CHUNK=200

Rscript code/15_combine_ordinal6_sameeffect_chunks.R
