# -----------------------------------------------------------------------------
# aim3_parallel.R - Run a single AIM 3 scenario in parallel
# 
# Usage: Rscript aim3_parallel.R <scenario_id>
# 
# This script runs one specific scenario from the aim3 simulation and saves
# both summary and detailed results to the parallel/results/ directory.
# -----------------------------------------------------------------------------

# --- 1. SETUP ---
message("--- Loading Packages ---")
pacman::p_load(tidyverse, broom, furrr, scales)

theme_set(theme_minimal() + theme(legend.position = "bottom"))
source("../code/functions/functions.R")

# --- 2. GET COMMAND LINE ARGUMENTS ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript aim3_parallel.R <scenario_id>")
}
scenario_id <- as.integer(args[1])

# --- 3. DEFINE SIMULATION PARAMETERS ---
message("--- Defining Simulation Parameters ---")

# --- Repetitions and Parallel Cores ---
n_reps_aim3 <- 200
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores() - 1))
if (is.na(n_cores) || n_cores < 1) n_cores <- 1
plan(multisession, workers = n_cores)
message(paste("--- Parallel processing enabled on", n_cores, "cores ---"))

# --- Define the two main OR scenarios ---
or_arrest <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
or_conservative <- c(A = 1.0, B = 2.0, C = 0.7, D = 1.2, E = 0.8)

scenario_definitions <- tribble(
  ~scenario_name, ~or_vector,
  "ARREST",        or_arrest,
  "Conservative",  or_conservative
)

# --- Shared Parameters ---
freq_arrest        <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)
p0_arrest_raw      <- c(A = 7/60, B = 0/52, C = 11/138, D = 11/69, E = 1/69)
p0_B_corrected     <- 0.5 / (52 + 0.5)
p0_arrest_adjusted <- p0_arrest_raw
p0_arrest_adjusted["B"] <- p0_B_corrected

# --- 4. DEFINE ALL SCENARIOS ---
sens_spec_scenarios <- tribble(
  ~test_type,                  ~sensitivity, ~specificity,
  "Perfect (100%)",            1.00,         1.00,
  "Near-Perfect",              0.99,         0.99,
  "High Sens/High Spec",       0.95,         0.95,
  "High Sens/Low Spec",        0.95,         0.70,
  "Low Sens/High Spec",        0.70,         0.95,
  "Balanced/Moderate",         0.80,         0.80
)
target_groups_aim3 <- c("B", "C", "D", "E")

all_aim3_scenarios <- expand_grid(
  scenario_definitions,
  sens_spec_scenarios,
  target_group = target_groups_aim3
) %>% mutate(seed = 1:n() + 4000)

# --- 5. VALIDATE SCENARIO ID ---
if (scenario_id < 1 || scenario_id > nrow(all_aim3_scenarios)) {
  stop(paste("Invalid scenario_id:", scenario_id, ". Must be between 1 and", nrow(all_aim3_scenarios)))
}

# --- 6. RUN SINGLE SCENARIO ---
message(paste("--- Running Scenario", scenario_id, "of", nrow(all_aim3_scenarios), "---"))

scenario_row <- all_aim3_scenarios[scenario_id, ]
scenario_name <- scenario_row$scenario_name
or_vector <- unlist(scenario_row$or_vector)  # Ensure named numeric vector
target_group <- scenario_row$target_group
sensitivity <- scenario_row$sensitivity
specificity <- scenario_row$specificity
test_type <- scenario_row$test_type
seed <- scenario_row$seed

message(sprintf("Scenario: %s, Group: %s, Test: %s (Sens: %.2f, Spec: %.2f)", 
                scenario_name, target_group, test_type, sensitivity, specificity))

# Set seed for reproducibility
set.seed(seed)

# Run the scenario
start_time <- Sys.time()

# Run the main function to get summary results
summary_result <- find_nns_for_scenario_sens_spec(
  target_group = target_group, sensitivity = sensitivity, specificity = specificity, seed = seed,
  or_vector = or_vector, p0_vector = p0_arrest_adjusted, freq_vector = freq_arrest,
  n_reps_per_calc = n_reps_aim3
) %>%
  dplyr::mutate(scenario_name = scenario_name, test_type = test_type)

# After summary_result is created, select and order columns to match closed_form tsv
# Handle case where power goal wasn't reached and some columns might be missing
if (is.na(summary_result$nns_needed)) {
  # Power goal wasn't reached - add missing columns with NA values
  summary_result <- summary_result %>%
    dplyr::mutate(
      true_beta = log(or_vector[target_group]),
      mean_beta_hat = NA_real_,
      bias = NA_real_,
      mse = NA_real_,
      rmse = NA_real_,
      sd_beta_hat = NA_real_,
      wrong_direction = NA_real_,
      significant_wrong = NA_real_
    ) %>%
    dplyr::select(scenario_name, test_type, sensitivity, specificity, target_group, 
                  nns_needed, nnr_corresponding, true_beta, mean_beta_hat, bias, mse, rmse, sd_beta_hat, wrong_direction, significant_wrong)
} else {
  # Power goal was reached - all columns should exist
  summary_result <- summary_result %>%
    dplyr::select(scenario_name, test_type, sensitivity, specificity, target_group, 
                  nns_needed, nnr_corresponding, true_beta, mean_beta_hat, bias, mse, rmse, sd_beta_hat, wrong_direction, significant_wrong)
}

# Also run the final scenario to get individual simulation results
if (!is.na(summary_result$nns_needed)) {
  final_nns <- summary_result$nns_needed
  final_scenario_result <- run_enrichment_scenario_sens_spec(
    n_screened = final_nns, target_group = target_group, sensitivity = sensitivity, specificity = specificity,
    or_vector = or_vector, p0_vector = p0_arrest_adjusted, freq_vector = freq_arrest,
    n_reps = n_reps_aim3 * 2, base_seed = seed  # Use same seed as power search
  )
  
  # Create detailed results dataframe
  true_beta <- log(or_vector[target_group])
  detailed_results <- tibble(
    scenario_id = scenario_id,
    scenario_name = scenario_name,
    test_type = test_type,
    target_group = target_group,
    sensitivity = sensitivity,
    specificity = specificity,
    sim_id = 1:length(final_scenario_result$beta_hats),
    true_beta = true_beta,
    empirical_beta = final_scenario_result$beta_hats,
    p_value = final_scenario_result$p_values,
    n_enrolled = final_scenario_result$mean_n_enrolled,
    bias = empirical_beta - true_beta,
    mse = (empirical_beta - true_beta)^2,
    significant = p_value < 0.05,
    wrong_direction = sign(empirical_beta) != sign(true_beta),
    significant_wrong = significant & wrong_direction
  )
  
  # Save detailed results to file
  detailed_filename <- sprintf("results/aim3_detailed_scenario_%03d.tsv", scenario_id)
  write_tsv(detailed_results, detailed_filename)
  message(paste("Detailed results saved to:", detailed_filename))
}

# Save summary results to file
summary_filename <- sprintf("results/aim3_summary_scenario_%03d.tsv", scenario_id)
write_tsv(summary_result, summary_filename)
message(paste("Summary results saved to:", summary_filename))

# Print completion message
duration <- difftime(Sys.time(), start_time, units = "mins")
message(sprintf("--- Scenario %d COMPLETE --- (Duration: %.1f minutes) ---", scenario_id, duration)) 