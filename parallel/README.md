# Parallel AIM 3 Simulations

This directory contains scripts to run AIM 3 simulations in parallel using SLURM, which is much faster than the sequential approach.

## Files

- `aim3_parallel.R` - Runs a single scenario (called by SLURM jobs)
- `aim3_collect.R` - Collects and combines results from all scenarios
- `aim3.sh` - Bash script that launches parallel SLURM jobs
- `test_parallel.R` - Test script to verify the approach works
- `results/` - Directory where individual scenario results are saved

## Usage

### 1. Test the Approach (Recommended First)

Test that everything works with a single scenario:

```bash
cd parallel
Rscript test_parallel.R
```

This will run scenario 1 (ARREST, Perfect test, Group B) with fewer replications to verify the setup works.

### 2. Run Full Parallel Simulation

Launch all 48 scenarios in parallel:

```bash
cd parallel
./aim3.sh
```

By default, this will run up to 16 concurrent jobs. You can specify a different limit:

```bash
./aim3.sh 8    # Run max 8 concurrent jobs
./aim3.sh 32   # Run max 32 concurrent jobs
```

### 3. Monitor Progress

Check job status:
```bash
squeue -u $USER
```

Check logs:
```bash
cd parallel/logs
ls -la
```

### 4. Collect Results

The `aim3.sh` script automatically runs the collection script when all jobs complete. If you need to run it manually:

```bash
cd parallel
Rscript aim3_collect.R
```

## Output Files

### Individual Scenario Results
- `results/aim3_summary_scenario_XXX.tsv` - Summary results for each scenario
- `results/aim3_detailed_scenario_XXX.tsv` - Detailed results for each scenario

### Final Combined Results
- `../results/tables/aim3_sens_spec_summary.tsv` - Combined summary results
- `../results/tables/aim3_detailed_all_scenarios.tsv` - Combined detailed results
- `../results/plots/aim3_nns_summary.pdf` - Final plot
- `../results/objects/aim3_plot.rds` - R plot object

## Enhanced Metrics

The parallel approach now includes additional metrics in the results:

### Summary Results Include:
- `true_beta` - True log odds ratio for the target group
- `mean_beta_hat` - Mean estimated beta across simulations
- `bias` - Mean bias (empirical - true)
- `mse` - Mean squared error
- `rmse` - Root mean squared error
- `sd_beta_hat` - Standard deviation of estimated betas
- `wrong_direction` - Proportion of simulations with wrong sign
- `significant_wrong` - Proportion of significant results with wrong sign

### Detailed Results Include:
- Individual simulation results with empirical beta, p-value, bias, MSE
- Flags for significance and wrong direction
- All metrics calculated per simulation

## SLURM Configuration

The default SLURM configuration in `aim3.sh`:
- `--time=02:00:00` - 2 hours per job
- `--mem=4G` - 4GB memory per job
- `--cpus-per-task=4` - 4 CPU cores per job
- `--partition=short` - Short partition

Adjust these in the script if needed for your cluster.

## Troubleshooting

### Job Fails
Check the error logs:
```bash
cd parallel/logs
cat aim3_scenario_X_*.err
```

### Missing Results
If some scenarios fail, you can rerun just those:
```bash
cd parallel
Rscript aim3_parallel.R <scenario_id>
```

### Collection Fails
If the collection script fails, check that all summary files exist:
```bash
cd parallel/results
ls aim3_summary_scenario_*.tsv | wc -l
```

Should show 48 files if all scenarios completed successfully.

## Performance

- **Sequential approach**: ~2-4 hours for all scenarios
- **Parallel approach**: ~10-30 minutes depending on cluster resources
- **Speedup**: 4-12x faster with 16 concurrent jobs

## Notes

- Each scenario uses 200 replications for power calculation and 400 for final metrics
- Results are automatically saved to both individual files and combined files
- The approach is robust - if one scenario fails, others continue
- All scenarios use the same random seeds as the original sequential approach 