#!/usr/bin/env bash
set -euo pipefail

cd /user/work/fh6520/bc4/sab_het
mkdir -p logs

base_tag="${1:-paper_final_independent_$(date +%Y%m%d_%H%M%S)}"
total_reps="${N_REPS_NONPO_TOTAL:-1000}"
chunk_count=8

if (( total_reps % chunk_count != 0 )); then
  echo "N_REPS_NONPO_TOTAL (${total_reps}) must be divisible by ${chunk_count}."
  exit 1
fi

reps_per_chunk=$(( total_reps / chunk_count ))

echo "Submitting independent simulation array:"
echo "  NONPO_BASE_TAG=${base_tag}"
echo "  N_REPS_NONPO_TOTAL=${total_reps}"
echo "  N_REPS_NONPO_PER_CHUNK=${reps_per_chunk}"

sim_array_job=$(
  sbatch --parsable \
    --export=ALL,NONPO_BASE_TAG="${base_tag}",N_REPS_NONPO_TOTAL="${total_reps}",N_REPS_NONPO_PER_CHUNK="${reps_per_chunk}" \
    slurm/14_paper_independent_sim_array.sh
)
echo "Simulation array job: ${sim_array_job}"

final_job=$(
  sbatch --parsable \
    --dependency=afterok:"${sim_array_job}" \
    --export=ALL,NONPO_BASE_TAG="${base_tag}",NONPO_CHUNK_COUNT="${chunk_count}" \
    slurm/13_paper_final_nonpo_finalize.sh
)
echo "Finalize job: ${final_job}"

echo "Submitted."
