#!/bin/bash
set -euo pipefail

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

# Quick sanity mode on login node:
# CORE_QUICK=1 Rscript code/00_run_core_paper.R
Rscript code/00_run_core_paper.R
