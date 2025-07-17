# -----------------------------------------------------------------------------
# Aim 4: Investigating Event Rate vs Proportion Effects
#
# This analysis investigates whether it's the overall event rates or the 
# proportion of events attributed to a specific group that matters more
# for enrichment trial design. Uses closed form solutions and focuses on
# arrest scenarios only.
# -----------------------------------------------------------------------------

# --- 1. SETUP ----------------------------------------------------------------

pacman::p_load(tidyverse, glue, readr, scales)

theme_set(theme_minimal())

# --- 2. BASE SCENARIO DATA --------------------------------------------------

# Base arrest scenario (from Aim 3)
or_arrest <- c(A = 1.00, B = 18.8, C = 0.79, D = 1.40, E = 0.30)
freq_arrest <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)
p0_arrest_raw <- c(A = 7/60, B = 0/52, C = 11/138, D = 11/69, E = 1/69)

# Apply continuity correction for group B
apply_zero_cc <- function(p_vec, freq_vec, add = 0.5) {
  out <- p_vec
  zero_idx <- which(p_vec == 0)
  if (length(zero_idx)) {
    out[zero_idx] <- add / (freq_vec[zero_idx] + add)
  }
  out
}

p0_arrest_adjusted <- apply_zero_cc(p0_arrest_raw, freq_arrest)

# Test scenarios (include worse tests)
test_scenarios <- tribble(
  ~test_type,            ~sensitivity, ~specificity,
  "Perfect (100%)",            1.00,         1.00,
  "Near-Perfect",       0.99,         0.99,
  "High Sens/High Spec",0.95,         0.95,
  "Balanced/Moderate",  0.80,         0.80,
  "Low Sens/High Spec", 0.70,         0.95,
  "High Sens/Low Spec", 0.95,         0.70
)

# Constants for power calculation
alpha          <- 0.05
z_alpha_half   <- qnorm(1 - alpha/2)
z_power_target <- qnorm(0.80)

# --- 3. HELPER FUNCTIONS -----------------------------------------------------

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

  # Contaminant Mix properties
  if (p_non_target > 0) {
    weights_mix <- p_vec[non_target_groups] / p_non_target
    p0_mix   <- sum(p0_vector[non_target_groups] * weights_mix)
    beta_mix <- sum(log(or_vector[non_target_groups]) * weights_mix)
    p1_mix   <- sum(get_or_p1(or_vector[non_target_groups], p0_vector[non_target_groups]) * weights_mix)
  } else {
    p0_mix   <- NA_real_
    beta_mix <- NA_real_
    p1_mix   <- NA_real_
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
  p_mix_in_enrolled  <- fp_rate / enrol_rate

  # Observed risks in enrolled cohort
  if (p_non_target > 0) {
    p0_obs   <- p_true_in_enrolled * p0_target + p_mix_in_enrolled * p0_mix
    p1_obs   <- p_true_in_enrolled * p1_target + p_mix_in_enrolled * p1_mix
  } else {
    p0_obs <- p0_target
    p1_obs <- p1_target
  }

  # Observed log-OR
  if (p0_obs %in% c(0,1) || p1_obs %in% c(0,1)) {
    beta_obs <- 0
  } else {
    beta_obs <- log( (p1_obs / (1 - p1_obs)) / (p0_obs / (1 - p0_obs)) )
  }

  list(beta_obs = beta_obs,
       p0_obs = p0_obs,
       p1_obs = p1_obs,
       enrol_rate = enrol_rate,
       p_true_in_enrolled = p_true_in_enrolled,
       p_mix_in_enrolled  = p_mix_in_enrolled,
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

# --- 4. SIMULATION 1: EVENT RATE MULTIPLIERS --------------------------------

# Multipliers for overall event rates (1.0x to 3x) - include baseline
event_rate_multipliers <- seq(1.0, 3.0, by = 0.1)

results_event_rate <- expand_grid(
  test_type = test_scenarios$test_type,
  sensitivity = test_scenarios$sensitivity,
  specificity = test_scenarios$specificity,
  multiplier = event_rate_multipliers
) %>%
  mutate(
    # Apply multiplier to all baseline risks
    p0_modified = map(multiplier, ~p0_arrest_adjusted * .x),
    
    # Calculate properties for each scenario
    cohort_props = pmap(list(sensitivity, specificity, p0_modified), 
                       function(sens, spec, p0_mod) {
                         calculate_mixed_cohort_properties(
                           target_group = "C", 
                           sensitivity = sens, 
                           specificity = spec,
                           or_vector = or_arrest, 
                           p0_vector = p0_mod, 
                           freq_vector = freq_arrest
                         )
                       }),
    
    # Extract results
    beta_obs = map_dbl(cohort_props, "beta_obs"),
    p0_obs = map_dbl(cohort_props, "p0_obs"),
    p1_obs = map_dbl(cohort_props, "p1_obs"),
    enrol_rate = map_dbl(cohort_props, "enrol_rate"),
    beta_target = map_dbl(cohort_props, "beta_target"),
    
    # Calculate NNS and bias
    nns_results = pmap(list(beta_obs, p0_obs, p1_obs, enrol_rate),
                      function(beta, p0, p1, enrol) {
                        calc_required_nns(beta, p0, p1, enrol, z_alpha_half, z_power_target)
                      }),
    
    nns = map_dbl(nns_results, "nns"),
    nnr = map_dbl(nns_results, "n_total"),
    bias_closed_form = beta_obs - beta_target
  ) %>%
  select(test_type, multiplier, nns, nnr, bias_closed_form) %>%
  mutate(
    scenario = "Event Rate Multiplier",
    subgroup = "C",
    multiplier_type = "Event Rate"
  ) %>%
  # Ensure we have exactly one row per test_type and multiplier combination
  distinct(test_type, multiplier, .keep_all = TRUE)

# --- 5. SIMULATION 2: PROPORTION MULTIPLIERS -------------------------------

# Multipliers for proportion of events attributed to group E (1.0x to 5x) - include baseline
proportion_multipliers <- seq(1.0, 5.0, by = 0.1)

results_proportion <- expand_grid(
  test_type = test_scenarios$test_type,
  sensitivity = test_scenarios$sensitivity,
  specificity = test_scenarios$specificity,
  multiplier = proportion_multipliers
) %>%
  mutate(
    # Modify frequency of group C while keeping others proportional
    freq_modified = map(multiplier, function(mult) {
      freq_mod <- freq_arrest
      freq_mod["C"] <- freq_arrest["C"] * mult
      # Renormalize to sum to 1
      freq_mod / sum(freq_mod)
    }),
    
    # Calculate properties for each scenario
    cohort_props = pmap(list(sensitivity, specificity, freq_modified), 
                       function(sens, spec, freq_mod) {
                         calculate_mixed_cohort_properties(
                           target_group = "C", 
                           sensitivity = sens, 
                           specificity = spec,
                           or_vector = or_arrest, 
                           p0_vector = p0_arrest_adjusted, 
                           freq_vector = freq_mod
                         )
                       }),
    
    # Extract results
    beta_obs = map_dbl(cohort_props, "beta_obs"),
    p0_obs = map_dbl(cohort_props, "p0_obs"),
    p1_obs = map_dbl(cohort_props, "p1_obs"),
    enrol_rate = map_dbl(cohort_props, "enrol_rate"),
    beta_target = map_dbl(cohort_props, "beta_target"),
    
    # Calculate NNS and bias
    nns_results = pmap(list(beta_obs, p0_obs, p1_obs, enrol_rate),
                      function(beta, p0, p1, enrol) {
                        calc_required_nns(beta, p0, p1, enrol, z_alpha_half, z_power_target)
                      }),
    
    nns = map_dbl(nns_results, "nns"),
    nnr = map_dbl(nns_results, "n_total"),
    bias_closed_form = beta_obs - beta_target
  ) %>%
  select(test_type, multiplier, nns, nnr, bias_closed_form) %>%
  mutate(
    scenario = "Proportion Multiplier",
    subgroup = "C",
    multiplier_type = "Proportion"
  ) %>%
  # Ensure we have exactly one row per test_type and multiplier combination
  distinct(test_type, multiplier, .keep_all = TRUE)

# --- 6. COMBINE RESULTS ----------------------------------------------------

results_aim4 <- bind_rows(results_event_rate, results_proportion) %>%
  arrange(test_type, scenario, multiplier)

# --- 7. OUTPUT -------------------------------------------------------------

# Ensure output directory exists
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# Save results
write_tsv(results_aim4, "results/tables/aim4_comparison_summary.tsv")

cat("Aim 4 results saved to results/tables/aim4_comparison_summary.tsv\n")

# Display summary
cat("\n=== AIM 4 SUMMARY ===\n")
cat("Event Rate Multipliers: 1.0x to 3.0x\n")
cat("Proportion Multipliers: 1.0x to 5.0x\n")
cat("Target Group: C\n")
cat("Test Types: Perfect, Near-Perfect, High Sens/High Spec, Balanced/Moderate, Low Sens/High Spec, High Sens/Low Spec\n")
cat("Total scenarios:", nrow(results_aim4), "\n\n")

# Show range of NNS values
cat("NNS Range:\n")
cat("  Event Rate Multipliers:", 
    range(results_event_rate$nns, na.rm = TRUE) %>% 
      format(big.mark = ",") %>% paste(collapse = " to "), "\n")
cat("  Proportion Multipliers:", 
    range(results_proportion$nns, na.rm = TRUE) %>% 
      format(big.mark = ",") %>% paste(collapse = " to "), "\n\n")

# Show range of bias values
cat("Bias Range:\n")
cat("  Event Rate Multipliers:", 
    range(results_event_rate$bias_closed_form, na.rm = TRUE) %>% 
      round(4) %>% paste(collapse = " to "), "\n")
cat("  Proportion Multipliers:", 
    range(results_proportion$bias_closed_form, na.rm = TRUE) %>% 
      round(4) %>% paste(collapse = " to "), "\n")

# Display first few rows
cat("\nFirst 10 rows of results:\n")
print(results_aim4 %>% head(10)) 