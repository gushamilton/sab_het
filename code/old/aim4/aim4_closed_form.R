# -----------------------------------------------------------------------------
# aim4_closed_form.R - AIM 4: Closed Form Analysis with Multipliers
# 
# This script uses the closed form Aim 3 approach to analyze:
# 1. Multipliers for overall event rates (1x to 3x)
# 2. Multipliers for proportion of events attributed to subgroup C (1x to 5x)
# 3. Different test types (sensitivity/specificity combinations)
# 4. Focus on subgroup C
# 
# Output: NNS, NNR, bias, and the two multiplier columns
# -----------------------------------------------------------------------------

# --- 1. SETUP ---
message("--- Loading Packages ---")
pacman::p_load(tidyverse, scales)

theme_set(theme_minimal() + theme(legend.position = "bottom"))

# Create results directories
dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)
dir.create("results/objects", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)

# --- 2. BASE PARAMETERS (from Aim 3) ---
message("--- Setting Base Parameters ---")

# Base odds ratios (ARREST scenario)
or_base <- c(A = 1.00, B = 18.8, C = 0.79, D = 1.40, E = 0.30)

# Base prevalence and baseline risks
freq_base <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)
p0_base_raw <- c(A = 7/60, B = 0/52, C = 11/138, D = 11/69, E = 1/69)

# Correctly define the adjusted baseline risks as used in the original Aim 3 analysis
p0_B_corrected <- 0.5 / (52 + 0.5)
p0_base_adjusted <- p0_base_raw
p0_base_adjusted["B"] <- p0_B_corrected

# Test types (sensitivity/specificity combinations)
test_types <- tribble(
  ~test_type,                  ~sensitivity, ~specificity,
  "Perfect (100%)",            1.00,         1.00,
  "Near-Perfect",              0.99,         0.99,
  "High Sens/High Spec",       0.95,         0.95,
  "High Sens/Low Spec",        0.95,         0.70,
  "Low Sens/High Spec",        0.70,         0.95,
  "Balanced/Moderate",         0.80,         0.80
)

# Constants for power calculation
alpha <- 0.05
z_alpha_half <- qnorm(1 - alpha/2)  # 1.96
z_power_target <- qnorm(0.80)       # 0.84

# --- 3. HELPER FUNCTIONS (from Aim 3) ---
get_or_p1 <- function(or, p0) { 
  (or * p0) / (1 - p0 + (or * p0)) 
}

calculate_mixed_cohort_properties <- function(target_group, sensitivity, specificity,
                                              or_vector, p0_vector, freq_vector) {
  
  p_vec <- freq_vector / sum(freq_vector)
  p_target <- p_vec[target_group]
  p_non_target <- 1 - p_target
  non_target_groups <- names(p_vec)[names(p_vec) != target_group]
  
  # True Target properties
  beta_target <- log(or_vector[target_group])
  p0_target <- p0_vector[target_group]
  p1_target <- get_or_p1(or_vector[target_group], p0_target)
  
  # Contaminant Mix
  if (p_non_target > 0) {
    weights_mix <- p_vec[non_target_groups] / p_non_target
    p0_mix <- sum(p0_vector[non_target_groups] * weights_mix)
    beta_mix <- sum(log(or_vector[non_target_groups]) * weights_mix)
    p1_mix <- sum(get_or_p1(or_vector[non_target_groups], p0_vector[non_target_groups]) * weights_mix)
  } else {
    p0_mix <- NA_real_
    beta_mix <- NA_real_
    p1_mix <- NA_real_
  }
  
  # Final enrolled cohort properties
  tp_rate <- sensitivity * p_target
  fp_rate <- (1 - specificity) * p_non_target
  enrol_rate <- tp_rate + fp_rate
  
  if (enrol_rate == 0) {
    return(list(beta_obs = 0, p0_obs = 0, p1_obs = 0, enrol_rate = 0,
                p_true_in_enrolled = NA_real_, p_mix_in_enrolled = NA_real_,
                beta_target = beta_target))
  }
  
  p_true_in_enrolled <- tp_rate / enrol_rate
  p_mix_in_enrolled <- fp_rate / enrol_rate
  
  # Observed risks in enrolled cohort
  if (p_non_target > 0) {
    p0_obs <- p_true_in_enrolled * p0_target + p_mix_in_enrolled * p0_mix
    p1_obs <- p_true_in_enrolled * p1_target + p_mix_in_enrolled * p1_mix
  } else {
    p0_obs <- p0_target
    p1_obs <- p1_target
  }
  
  # Observed log-OR
  if (p0_obs %in% c(0,1) || p1_obs %in% c(0,1)) {
    beta_obs <- 0
  } else {
    beta_obs <- log((p1_obs / (1 - p1_obs)) / (p0_obs / (1 - p0_obs)))
  }
  
  list(beta_obs = beta_obs,
       p0_obs = p0_obs,
       p1_obs = p1_obs,
       enrol_rate = enrol_rate,
       p_true_in_enrolled = p_true_in_enrolled,
       p_mix_in_enrolled = p_mix_in_enrolled,
       beta_target = beta_target)
}

calc_required_nns <- function(beta_obs, p0_obs, p1_obs, enrol_rate,
                              z_alpha_half, z_power_target) {
  
  if (abs(beta_obs) < 1e-6 || is.na(beta_obs) || enrol_rate == 0 ||
      p0_obs == 0 || p0_obs == 1 || p1_obs == 0 || p1_obs == 1) {
    return(list(n_total = Inf, nns = Inf))
  }
  
  var_term <- (1 / (p0_obs * (1 - p0_obs))) + (1 / (p1_obs * (1 - p1_obs)))
  if (!is.finite(var_term)) return(list(n_total = Inf, nns = Inf))
  
  n_total <- 2 * (z_alpha_half + z_power_target)^2 * var_term / (beta_obs^2)
  nns <- ceiling(n_total / enrol_rate)
  
  list(n_total = ceiling(n_total), nns = nns)
}

# --- 4. MULTIPLIER SCENARIOS ---
message("--- Defining Multiplier Scenarios ---")

# Overall event rate multipliers (1x to 2x)
event_rate_multipliers <- seq(1, 2, by = 0.1)

# Proportion of events attributed to subgroup C multipliers (1x to 2x)
proportion_multipliers <- seq(1, 2, by = 0.1)

# Target group is C
target_group <- "C"

# --- 5. MAIN CALCULATION ---
message("--- Running Aim 4 Calculations ---")

# Create all combinations
all_scenarios <- expand_grid(
  test_types,
  event_rate_multiplier = event_rate_multipliers,
  proportion_multiplier = proportion_multipliers
)

# Function to apply multipliers and calculate results
calculate_with_multipliers <- function(test_type, sensitivity, specificity, 
                                       event_rate_multiplier, proportion_multiplier) {
  
  # Apply event rate multiplier to all baseline risks
  p0_multiplied <- p0_base_adjusted * event_rate_multiplier
  # Cap at 0.99 to avoid probabilities > 1
  p0_multiplied <- pmin(p0_multiplied, 0.99)
  
  # Apply proportion multiplier to subgroup C's baseline risk
  p0_multiplied["C"] <- p0_multiplied["C"] * proportion_multiplier
  p0_multiplied["C"] <- pmin(p0_multiplied["C"], 0.99)
  
  # Calculate mixed cohort properties
  cohort_props <- calculate_mixed_cohort_properties(
    target_group = target_group,
    sensitivity = sensitivity,
    specificity = specificity,
    or_vector = or_base,
    p0_vector = p0_multiplied,
    freq_vector = freq_base
  )
  
  # Calculate required NNS/NNR
  reqs <- calc_required_nns(
    cohort_props$beta_obs,
    cohort_props$p0_obs,
    cohort_props$p1_obs,
    cohort_props$enrol_rate,
    z_alpha_half,
    z_power_target
  )
  
  # Calculate bias
  beta_target <- cohort_props$beta_target
  beta_obs <- cohort_props$beta_obs
  bias <- beta_obs - beta_target
  
  # Return results
  tibble(
    test_type = test_type,
    event_rate_multiplier = event_rate_multiplier,
    proportion_multiplier = proportion_multiplier,
    nns = reqs$nns,
    nnr = reqs$n_total,
    bias = bias
  )
}

# Run all scenarios
message(paste("Running", nrow(all_scenarios), "scenarios..."))
results <- pmap_dfr(all_scenarios, calculate_with_multipliers)

# --- 6. SAVE RESULTS ---
message("--- Saving Results ---")

write_tsv(results, "results/tables/aim4_closed_form_results.tsv")
message("Results saved to results/tables/aim4_closed_form_results.tsv")

# Summary statistics
summary_stats <- results %>%
  group_by(test_type) %>%
  summarise(
    n_scenarios = n(),
    mean_nns = mean(nns, na.rm = TRUE),
    median_nns = median(nns, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

message("--- Aim 4 Closed Form Analysis Complete ---") 