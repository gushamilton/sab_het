# -----------------------------------------------------------------------------
# test_parallel.R - Test the parallel AIM 3 approach
# 
# This script tests the parallel approach by running just one scenario
# to make sure everything works correctly.
# -----------------------------------------------------------------------------

# --- 1. SETUP ---
message("--- Loading Packages ---")
pacman::p_load(tidyverse, broom, furrr, scales)

theme_set(theme_minimal() + theme(legend.position = "bottom"))
source("../code/functions/functions.R")

# --- 2. TEST PARAMETERS ---
message("--- Testing parallel approach ---")

# Test with scenario 1 (ARREST, Perfect test, Group B)
scenario_id <- 1

# --- 3. DEFINE SCENARIO PARAMETERS ---
or_arrest <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
or_conservative <- c(A = 1.0, B = 2.0, C = 0.7, D = 1.2, E = 0.8)

scenario_definitions <- tribble(
  ~scenario_name, ~or_vector,
  "ARREST",        or_arrest,
  "Conservative",  or_conservative
)

freq_arrest        <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)
p0_arrest_raw      <- c(A = 7/60, B = 0/52, C = 11/138, D = 11/69, E = 1/69)
p0_B_corrected     <- 0.5 / (52 + 0.5)
p0_arrest_adjusted <- p0_arrest_raw
p0_arrest_adjusted["B"] <- p0_B_corrected

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

# --- 4. RUN TEST SCENARIO ---
message(paste("--- Testing Scenario", scenario_id, "---"))

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

# Run with fewer replications for testing
n_reps_test <- 50

# Run the scenario
start_time <- Sys.time()

# Run the main function to get summary results
summary_result <- find_nns_for_scenario_sens_spec(
  target_group = target_group, sensitivity = sensitivity, specificity = specificity, seed = seed,
  or_vector = or_vector, p0_vector = p0_arrest_adjusted, freq_vector = freq_arrest,
  n_reps_per_calc = n_reps_test
) %>%
  dplyr::mutate(scenario_name = scenario_name, test_type = test_type)

# After summary_result is created, select and order columns to match closed_form tsv
summary_result <- summary_result %>%
  dplyr::select(scenario_name, test_type, sensitivity, specificity, target_group, nns_needed, nnr_corresponding, true_beta, mean_beta_hat, bias, mse, rmse, sd_beta_hat, wrong_direction, significant_wrong)

# Print results
message("--- Test Results ---")
print(summary_result)

# Print completion message
duration <- difftime(Sys.time(), start_time, units = "mins")
message(sprintf("--- Test COMPLETE --- (Duration: %.1f minutes) ---", duration))
message("If this works, the parallel approach is ready to use!") 