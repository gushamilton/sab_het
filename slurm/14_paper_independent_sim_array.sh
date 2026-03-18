#!/bin/bash
#SBATCH --job-name=sab_paper_sim
#SBATCH --partition=compute
#SBATCH --time=18:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-9
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

task_id="${SLURM_ARRAY_TASK_ID}"
base_tag="${NONPO_BASE_TAG:-paper_final_independent_array}"
chunk_count=8
n_reps_total="${N_REPS_NONPO_TOTAL:-1000}"
reps_per_chunk="${N_REPS_NONPO_PER_CHUNK:-$(( n_reps_total / chunk_count ))}"

if [[ "$(( reps_per_chunk * chunk_count ))" -ne "$n_reps_total" ]]; then
  echo "N_REPS_NONPO_TOTAL (${n_reps_total}) must be divisible by ${chunk_count}."
  exit 1
fi

if [[ "${task_id}" -eq 1 ]]; then
  echo "[$(date -Is)] Task ${task_id}: running core pipeline (no nonPO, no final assets)"
  export PAPER_ACCURACY="${PAPER_ACCURACY:-1.00}"
  export PAPER_BUILD_NONPO_PLOTS=0
  export PAPER_SKIP_NONPO=1
  export PAPER_SKIP_ASSETS=1
  export NONPO_RUN_TAG="${base_tag}"
  export SAVE_COHORT_RAW=1
  Rscript code/00_run_paper_final.R
  echo "[$(date -Is)] Task ${task_id}: core pipeline complete"
  exit 0
fi

chunk_id="$(( task_id - 1 ))"
rep_offset="$(( (chunk_id - 1) * reps_per_chunk ))"
chunk_tag="${base_tag}_chunk${chunk_id}"

echo "[$(date -Is)] Task ${task_id}: nonPO chunk ${chunk_id}/${chunk_count}, tag=${chunk_tag}, reps=${reps_per_chunk}, offset=${rep_offset}"

export RUN_TAG="${chunk_tag}"
export N_REPS_DOOR="${reps_per_chunk}"
export REP_OFFSET="${rep_offset}"
export SKIP_CALIBRATION=1

if [[ -n "${CALIBRATION_N:-}" ]]; then export CALIBRATION_N; fi
if [[ -n "${NONPO_SAMPLE_SIZES:-}" ]]; then export SAMPLE_SIZES="${NONPO_SAMPLE_SIZES}"; fi

Rscript code/15_run_ordinal_nonPO_comparison.R

echo "[$(date -Is)] Task ${task_id}: nonPO chunk ${chunk_id} complete"
