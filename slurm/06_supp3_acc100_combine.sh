#!/bin/bash
#SBATCH --job-name=sab_s3_a100_merge
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

export ACCURACY=1
export ORDINAL_POINTS=6
export CHUNK_COUNT="${CHUNK_COUNT:-15}"

Rscript code/09_combine_supp3_chunks.R
