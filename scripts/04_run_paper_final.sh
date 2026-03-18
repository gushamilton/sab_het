#!/bin/bash
set -euo pipefail

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

# Full manuscript build:
# Rscript code/00_run_paper_final.R
#
# Quick login-node smoke:
# PAPER_QUICK=1 Rscript code/00_run_paper_final.R
Rscript code/00_run_paper_final.R
