#!/bin/bash
#SBATCH --job-name=sab_nonpo_chunk
#SBATCH --partition=compute
#SBATCH --time=12:00:00
#SBATCH --mem=16G
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

base_tag="${NONPO_BASE_TAG:-paper_final_array}"
n_reps_total="${N_REPS_NONPO_TOTAL:-1000}"
reps_per_chunk="${N_REPS_NONPO_PER_CHUNK:-$(( n_reps_total / chunk_count ))}"
rep_offset="$(( (chunk_id - 1) * reps_per_chunk ))"

if [[ "$(( reps_per_chunk * chunk_count ))" -ne "$n_reps_total" ]]; then
  echo "N_REPS_NONPO_TOTAL (${n_reps_total}) must be divisible by array task count (${chunk_count})."
  exit 1
fi

chunk_tag="${base_tag}_chunk${chunk_id}"

echo "[$(date -Is)] Starting nonPO chunk ${chunk_id}/${chunk_count} tag=${chunk_tag} reps=${reps_per_chunk} offset=${rep_offset}"

export RUN_TAG="${chunk_tag}"
export N_REPS_DOOR="${reps_per_chunk}"
export REP_OFFSET="${rep_offset}"
export SKIP_CALIBRATION=1

if [[ -n "${CALIBRATION_N:-}" ]]; then export CALIBRATION_N; fi
if [[ -n "${NONPO_SAMPLE_SIZES:-}" ]]; then export SAMPLE_SIZES="${NONPO_SAMPLE_SIZES}"; fi

Rscript code/15_run_ordinal_nonPO_comparison.R

echo "[$(date -Is)] Finished nonPO chunk ${chunk_id}/${chunk_count} tag=${chunk_tag}"
