#!/usr/bin/env bash
set -euo pipefail

cd /user/work/fh6520/bc4/sab_het
mkdir -p logs

base_tag="${1:-paper_final_array_$(date +%Y%m%d_%H%M%S)}"
chunks="${NONPO_CHUNK_COUNT:-8}"
total_reps="${N_REPS_NONPO_TOTAL:-1000}"

if (( total_reps % chunks != 0 )); then
  echo "N_REPS_NONPO_TOTAL (${total_reps}) must be divisible by NONPO_CHUNK_COUNT (${chunks})."
  exit 1
fi

reps_per_chunk=$(( total_reps / chunks ))

echo "Submitting with:"
echo "  NONPO_BASE_TAG=${base_tag}"
echo "  NONPO_CHUNK_COUNT=${chunks}"
echo "  N_REPS_NONPO_TOTAL=${total_reps}"
echo "  N_REPS_NONPO_PER_CHUNK=${reps_per_chunk}"

core_job_id=$(
  sbatch --parsable \
    --export=ALL,NONPO_BASE_TAG="${base_tag}" \
    slurm/11_paper_final_core_no_nonpo.sh
)
echo "Core job: ${core_job_id}"

array_job_id=$(
  sbatch --parsable \
    --array=1-"${chunks}" \
    --export=ALL,NONPO_BASE_TAG="${base_tag}",N_REPS_NONPO_TOTAL="${total_reps}",N_REPS_NONPO_PER_CHUNK="${reps_per_chunk}",NONPO_CHUNK_COUNT="${chunks}" \
    slurm/12_paper_final_nonpo_chunks_array.sh
)
echo "NonPO array job: ${array_job_id}"

final_job_id=$(
  sbatch --parsable \
    --dependency=afterok:"${core_job_id}:${array_job_id}" \
    --export=ALL,NONPO_BASE_TAG="${base_tag}",NONPO_CHUNK_COUNT="${chunks}" \
    slurm/13_paper_final_nonpo_finalize.sh
)
echo "Finalize job: ${final_job_id}"

echo "Submission chain complete."
