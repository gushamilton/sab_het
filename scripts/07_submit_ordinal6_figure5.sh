#!/usr/bin/env bash
set -euo pipefail

cd /user/work/fh6520/bc4/sab_het
mkdir -p logs

base_tag="${1:-ordinal6_figure5_$(date +%Y%m%d_%H%M%S)}"
chunks="${ORD6_CHUNK_COUNT:-8}"
total_reps="${N_REPS_ORD6_TOTAL:-1000}"

if (( total_reps % chunks != 0 )); then
  echo "N_REPS_ORD6_TOTAL (${total_reps}) must be divisible by ORD6_CHUNK_COUNT (${chunks})."
  exit 1
fi

reps_per_chunk=$(( total_reps / chunks ))

echo "Submitting 6-point Figure 5 with:"
echo "  ORD6_BASE_TAG=${base_tag}"
echo "  ORD6_CHUNK_COUNT=${chunks}"
echo "  N_REPS_ORD6_TOTAL=${total_reps}"
echo "  N_REPS_ORD6_PER_CHUNK=${reps_per_chunk}"

array_job_id=$(
  sbatch --parsable \
    --array=1-"${chunks}" \
    --export=ALL,ORD6_BASE_TAG="${base_tag}",N_REPS_ORD6_TOTAL="${total_reps}",N_REPS_ORD6_PER_CHUNK="${reps_per_chunk}",ORD6_CHUNK_COUNT="${chunks}" \
    slurm/16_ordinal6_figure5_chunks_array.sh
)
echo "6-point array job: ${array_job_id}"

final_job_id=$(
  sbatch --parsable \
    --dependency=afterok:"${array_job_id}" \
    --export=ALL,ORD6_BASE_TAG="${base_tag}",ORD6_CHUNK_COUNT="${chunks}" \
    slurm/17_ordinal6_figure5_finalize.sh
)
echo "6-point finalize job: ${final_job_id}"

echo "6-point Figure 5 submission chain complete."
