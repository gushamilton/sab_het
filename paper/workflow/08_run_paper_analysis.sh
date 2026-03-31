#!/usr/bin/env bash
set -euo pipefail

cd /user/work/fh6520/bc4/sab_het

module purge
module load languages/R/4.5.1

Rscript code/00_run_paper_final.R

echo "Rebuilt paper analysis assets under results/paper/final/"
