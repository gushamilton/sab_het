#!/bin/bash
set -euo pipefail

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

Rscript code/03_run_main_binary.R
