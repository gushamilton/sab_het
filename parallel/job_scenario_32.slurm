#!/bin/bash
#SBATCH --job-name=aim3_scenario_32
#SBATCH --output=logs/aim3_scenario_32_%%j.log
#SBATCH --error=logs/aim3_scenario_32_%%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=16


# Load R module (adjust as needed for your cluster)
module load languages/R

# Set working directory to the parallel folder
cd /user/work/fh6520/sab_het/parallel

# Run the scenario
echo "Starting scenario 32 at \\$(date)"
Rscript aim3_parallel.R 32
echo "Completed scenario 32 at \\$(date)"
