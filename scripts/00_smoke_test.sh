#!/bin/bash
set -euo pipefail

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

export N_REPS_MAIN=50
export N_REPS_ORDINAL=50
export N_REPS_COHORT=50
export SAMPLE_SIZES=500,1000

Rscript code/03_run_main_binary.R
Rscript code/04_run_enrichment.R
Rscript code/05_run_cohort_variation.R
Rscript code/06_run_ordinal_comparison.R
Rscript code/07_build_gt_tables.R
