#!/bin/bash
#SBATCH --job-name=sab_nonpo_final
#SBATCH --partition=compute
#SBATCH --time=08:00:00
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

base_tag="${NONPO_BASE_TAG:-paper_final_array}"
chunk_count="${NONPO_CHUNK_COUNT:-8}"

echo "[$(date -Is)] Combining nonPO chunks for base tag ${base_tag}"

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

echo "[$(date -Is)] Building nonPO summaries and manuscript assets for tag ${base_tag}"

RUN_TAG="${base_tag}" Rscript code/22_summarise_relative_rmse_nonpo.R
RUN_TAG="${base_tag}" Rscript code/20_plot_ordinal_po_nonpo_patchwork.R
RUN_TAG="${base_tag}" Rscript code/21_plot_ordinal_po_nonpo_centered_bias_patchwork.R
RUN_TAG="${base_tag}" Rscript code/23_plot_split_bias_by_model_nonpo.R

PAPER_ACCURACY="${PAPER_ACCURACY:-1.00}" NONPO_RUN_TAG="${base_tag}" Rscript code/28_build_paper_final_assets.R

echo "[$(date -Is)] Finished nonPO finalize + paper asset build"
