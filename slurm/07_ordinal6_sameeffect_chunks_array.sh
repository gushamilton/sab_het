#!/bin/bash
#SBATCH --job-name=sab_ord6_chunks
#SBATCH --partition=compute
#SBATCH --time=06:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH --account=sscm013902
#SBATCH --array=1-40
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge
module load languages/R/4.5.1

cd /user/work/fh6520/bc4/sab_het

SAMPLE_SIZES_ALL=(500 1000 1500 2000 3000 5000 10000 20000)
REP_CHUNKS=5
REPS_PER_CHUNK=200

TASK_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE_INDEX=$((TASK_INDEX / REP_CHUNKS))
CHUNK_INDEX=$((TASK_INDEX % REP_CHUNKS + 1))

if [ "${SAMPLE_INDEX}" -lt 0 ] || [ "${SAMPLE_INDEX}" -ge "${#SAMPLE_SIZES_ALL[@]}" ]; then
  echo "Invalid SAMPLE_INDEX=${SAMPLE_INDEX} from SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

TARGET_N="${SAMPLE_SIZES_ALL[${SAMPLE_INDEX}]}"
REP_OFFSET=$(((CHUNK_INDEX - 1) * REPS_PER_CHUNK))
RUN_TAG="n${TARGET_N}_r${CHUNK_INDEX}"

export SAMPLE_SIZES="${TARGET_N}"
export N_REPS_ORD6="${REPS_PER_CHUNK}"
export REP_OFFSET="${REP_OFFSET}"
export RUN_TAG="${RUN_TAG}"
export SKIP_CALIBRATION=1
export CALIBRATION_N=200000

echo "Running SAMPLE_SIZE=${TARGET_N}, CHUNK=${CHUNK_INDEX}, REP_OFFSET=${REP_OFFSET}, RUN_TAG=${RUN_TAG}"
Rscript code/14_run_ordinal6_same_effect_comparison.R
