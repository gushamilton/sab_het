#!/bin/bash
#SBATCH --job-name=sab_ord6_final
#SBATCH --partition=short,compute
#SBATCH --time=02:00:00
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

base_tag="${ORD6_BASE_TAG:-ordinal6_figure5}"
chunk_count="${ORD6_CHUNK_COUNT:-8}"

chunk_tags=""
for i in $(seq 1 "${chunk_count}"); do
  tag="${base_tag}_chunk${i}"
  if [[ -z "${chunk_tags}" ]]; then
    chunk_tags="${tag}"
  else
    chunk_tags="${chunk_tags},${tag}"
  fi
done

export CHUNK_TAGS="${chunk_tags}"
export FINAL_RUN_TAG="${base_tag}"
Rscript code/29_combine_ordinal_nonPO_chunks.R

RUN_TAG="${base_tag}" Rscript code/22_summarise_relative_rmse_nonpo.R
RUN_TAG="${base_tag}" Rscript code/20_plot_ordinal_po_nonpo_patchwork.R
RUN_TAG="${base_tag}" Rscript code/21_plot_ordinal_po_nonpo_centered_bias_patchwork.R
RUN_TAG="${base_tag}" Rscript code/23_plot_split_bias_by_model_nonpo.R

mkdir -p results/paper/final/figures
cp "results/supp/ordinal_nonPO_comparison/plots/ordinal_nonPO_patchwork_${base_tag}.pdf" \
   "results/paper/final/figures/Figure5_ordinal_binary_summary_6point.pdf"

echo "[$(date -Is)] Finished 6-point Figure 5 finalize for ${base_tag}"
