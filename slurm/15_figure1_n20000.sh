#!/bin/bash
#SBATCH --job-name=sab_fig1_n20000
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
mkdir -p logs

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

echo "[fig1] Rebuilding Figure 1 with N=20000"
FIG1_N=20000 Rscript code/11_run_figure1_accuracy_scan.R

echo "[fig1] Refreshing paper final assets"
PAPER_ACCURACY=1.00 NONPO_RUN_TAG=paper_final_independent_20260313 \
  Rscript code/28_build_paper_final_assets.R

echo "[fig1] Done"
