#!/bin/bash
#SBATCH --job-name=sab_supp
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
mkdir -p logs

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

Rscript code/04_run_enrichment.R
Rscript code/05_run_cohort_variation.R
Rscript code/06_run_ordinal_comparison.R
Rscript code/07_build_gt_tables.R
