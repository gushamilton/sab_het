#!/bin/bash
#SBATCH --job-name=sab_ord6_chunk
#SBATCH --partition=short,compute
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-8
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

chunk_id="${SLURM_ARRAY_TASK_ID}"
chunk_count="${SLURM_ARRAY_TASK_COUNT}"

base_tag="${ORD6_BASE_TAG:-ordinal6_figure5}"
n_reps_total="${N_REPS_ORD6_TOTAL:-1000}"
reps_per_chunk="${N_REPS_ORD6_PER_CHUNK:-$(( n_reps_total / chunk_count ))}"
rep_offset="$(( (chunk_id - 1) * reps_per_chunk ))"

if [[ "$(( reps_per_chunk * chunk_count ))" -ne "$n_reps_total" ]]; then
  echo "N_REPS_ORD6_TOTAL (${n_reps_total}) must be divisible by array task count (${chunk_count})."
  exit 1
fi

chunk_tag="${base_tag}_chunk${chunk_id}"

echo "[$(date -Is)] Starting 6-point Figure 5 chunk ${chunk_id}/${chunk_count} tag=${chunk_tag} reps=${reps_per_chunk} offset=${rep_offset}"

export ORDINAL_POINTS=6
export RUN_TAG="${chunk_tag}"
export N_REPS_DOOR="${reps_per_chunk}"
export REP_OFFSET="${rep_offset}"
export SKIP_CALIBRATION=1

if [[ -n "${CALIBRATION_N:-}" ]]; then export CALIBRATION_N; fi
if [[ -n "${NONPO_SAMPLE_SIZES:-}" ]]; then export SAMPLE_SIZES="${NONPO_SAMPLE_SIZES}"; fi

Rscript code/15_run_ordinal_nonPO_comparison.R

echo "[$(date -Is)] Finished 6-point Figure 5 chunk ${chunk_id}/${chunk_count} tag=${chunk_tag}"
