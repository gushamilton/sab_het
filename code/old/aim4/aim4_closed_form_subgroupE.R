# -----------------------------------------------------------------------------
# 02_closed_form_aim3_fixed.R
#
# Analytic (closed-form) approximation for Aim 3: compute the minimum number of
# patients to screen (NNS) required to achieve 80 % power for an enrichment
# trial. This version corrects two internal-logic errors identified in review:
#   1. The sample-size formula for log-OR omitted a factor of 2.
#   2. The observed (diluted) log-OR in the enrolled mixed cohort must be
#      derived from the pooled observed event rates, not as a weighted average
#      of subgroup log-ORs.
# Minor defensive/clarity edits are also included (see "FIX:" comments below).
# -----------------------------------------------------------------------------

# --- 1. SETUP ----------------------------------------------------------------

pacman::p_load(tidyverse, glue, readr, scales)

theme_set(theme_minimal())

# --- 2. SCENARIOS & INPUT DATA ----------------------------------------------

# Odds-ratio scenarios (same as in 01_run_simulations.R)
or_arrest        <- c(A = 1.00, B = 18.8, C = 0.79, D = 1.40, E = 0.30)
or_conservative  <- c(A = 1.00, B =  2.0, C = 0.70, D = 1.20, E = 0.80)

scenario_definitions <- tribble(
  ~scenario_name, ~or_vector,
  "ARREST",        or_arrest,
  "Conservative",  or_conservative
)

# Prevalence and baseline risks (from 01_run_simulations.R)
freq_arrest <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)
p0_arrest_raw      <- c(A = 7/60, B = 0/52, C = 11/138, D = 11/69, E = 1/69)

# --- Continuity correction helper -------------------------------------------
# FIX: Generalize the previously hard-coded continuity correction that was only
# applied to Group B. We replicate the original correction form (0.5/(n+0.5))
# whenever a 0 event count would produce p == 0.
apply_zero_cc <- function(p_vec, freq_vec, add = 0.5) {
  # p_vec: vector of raw risks (events/n)
  # freq_vec: denominators (n) on the same scale/order as p_vec
  # For elements where p_vec == 0, replace with add / (n + add)
  out <- p_vec
  zero_idx <- which(p_vec == 0)
  if (length(zero_idx)) {
    out[zero_idx] <- add / (freq_vec[zero_idx] + add)
  }
  out
}

# Apply the continuity correction
p0_arrest_adjusted <- apply_zero_cc(p0_arrest_raw, freq_arrest)

# Sensitivity/specificity test grid (including the “Perfect” test)
sens_spec_scenarios <- tribble(
  ~test_type,            ~sensitivity, ~specificity,
  "Perfect (100%)",            1.00,         1.00,
  "Near-Perfect",       0.99,         0.99,
  "High Sens/High Spec",0.95,         0.95,
  "High Sens/Low Spec", 0.95,         0.70,
  "Low Sens/High Spec", 0.70,         0.95,
  "Balanced/Moderate",  0.80,         0.80
)

target_groups <- c("B", "C", "D", "E")

# Constants for power calculation
alpha          <- 0.05               # two-sided
z_alpha_half   <- qnorm(1 - alpha/2) # 1.96
z_power_target <- qnorm(0.80)        # 0.84

# --- 3. HELPER FUNCTIONS -----------------------------------------------------

get_or_p1 <- function(or, p0) { (or * p0) / (1 - p0 + (or * p0)) }

# FIX: calculate_mixed_cohort_properties now returns beta_obs computed from the
# pooled observed risks (p0_obs, p1_obs), *not* a weighted average of log-ORs.
# This yields an internally consistent effect size with the variance calculation.
calculate_mixed_cohort_properties <- function(target_group, sensitivity, specificity,
                                              or_vector, p0_vector, freq_vector) {

  p_vec <- freq_vector / sum(freq_vector)
  p_target <- p_vec[target_group]
  p_non_target <- 1 - p_target
  non_target_groups <- names(p_vec)[names(p_vec) != target_group]

  # --- 1. Calculate properties of the two components: True Target vs. Contaminants
  # True Target
  beta_target <- log(or_vector[target_group])
  p0_target <- p0_vector[target_group]
  p1_target <- get_or_p1(or_vector[target_group], p0_target)

  # Contaminant Mix -----------------------------------------------------------
  if (p_non_target > 0) {
    # Weight each non-target group by its prevalence within the non-target pool
    weights_mix <- p_vec[non_target_groups] / p_non_target

    # Weighted average properties of the contaminant mix
    p0_mix   <- sum(p0_vector[non_target_groups] * weights_mix)
    beta_mix <- sum(log(or_vector[non_target_groups]) * weights_mix)
    p1_mix   <- sum(get_or_p1(or_vector[non_target_groups], p0_vector[non_target_groups]) * weights_mix)
  } else {
    # Degenerate case: no contaminants
    p0_mix   <- NA_real_
    beta_mix <- NA_real_
    p1_mix   <- NA_real_
  }

  # --- 2. Calculate properties of the final enrolled (mixed) cohort ----------
  # Proportion of enrolled cohort that is true target vs. contaminant
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

  # Observed (diluted) baseline & treated risks in the final enrolled cohort
  if (p_non_target > 0) {
    p0_obs   <- p_true_in_enrolled * p0_target + p_mix_in_enrolled * p0_mix
    p1_obs   <- p_true_in_enrolled * p1_target + p_mix_in_enrolled * p1_mix
  } else {
    # All enrolled are true targets
    p0_obs <- p0_target
    p1_obs <- p1_target
  }

  # FIX: derive observed log-OR from pooled risks so that
  # beta_obs == log( (p1_obs/(1-p1_obs)) / (p0_obs/(1-p0_obs)) ).
  if (p0_obs %in% c(0,1) || p1_obs %in% c(0,1)) {
    beta_obs <- 0  # will be trapped upstream in sample-size fn
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

  # Guard against degenerate inputs ------------------------------------------------
  if (abs(beta_obs) < 1e-6 || is.na(beta_obs) || enrol_rate == 0 ||
      p0_obs == 0 || p0_obs == 1 || p1_obs == 0 || p1_obs == 1) {
    return(list(n_total = Inf, nns = Inf))
  }

  # Variance for log-OR based on the OBSERVED event rates in the mixed cohort
  # var_term is the sum of inverse cell probabilities; true Wald variance adds
  # a factor 2/n_total under balanced allocation. (FIX: multiply by 2 below.)
  var_term <- (1 / (p0_obs * (1 - p0_obs))) + (1 / (p1_obs * (1 - p1_obs)))
  if (!is.finite(var_term)) return(list(n_total = Inf, nns = Inf))

  # Standard sample size formula, assuming 1:1 allocation.
  # n_total = 2 * (z_{alpha/2} + z_{power})^2 * var_term / beta^2
  n_total <- 2 * (z_alpha_half + z_power_target)^2 * var_term / (beta_obs^2)  # FIX

  nns <- ceiling(n_total / enrol_rate)

  list(n_total = ceiling(n_total), nns = nns)
}

# --- 4. MAIN CALCULATION -----------------------------------------------------

# Define multipliers
event_rate_multipliers <- seq(1, 2, by = 0.1)
proportion_multipliers <- seq(1, 2, by = 0.1)

# Create all combinations for the new analysis
all_cases_new <- expand_grid(
  scenario_name = "ARREST", # Focus on ARREST scenario
  sens_spec_scenarios,
  target_group = "E", # Focus on subgroup E
  event_rate_multiplier = event_rate_multipliers,
  proportion_multiplier = proportion_multipliers
)


results_new <- pmap_dfr(all_cases_new, function(scenario_name, test_type, sensitivity, specificity,
                                                target_group, event_rate_multiplier, proportion_multiplier) {

  # Get the base or_vector for the scenario
  or_vector <- scenario_definitions$or_vector[scenario_definitions$scenario_name == scenario_name][[1]]

  # Apply multipliers to a fresh copy of the base p0 vector each time
  p0_modified <- p0_arrest_adjusted
  p0_modified <- p0_modified * event_rate_multiplier
  p0_modified[target_group] <- p0_modified[target_group] * proportion_multiplier
  p0_modified <- pmin(p0_modified, 0.99) # Cap probabilities at 0.99


  # Calculate diluted effect size AND diluted event rates for the enrolled cohort
  cohort_props <- calculate_mixed_cohort_properties(target_group, sensitivity, specificity,
                                                    or_vector, p0_modified, freq_arrest)

  # Use these observed properties to calculate required NNS/NNR
  reqs <- calc_required_nns(cohort_props$beta_obs,
                            cohort_props$p0_obs,
                            cohort_props$p1_obs,
                            cohort_props$enrol_rate,
                            z_alpha_half, z_power_target)

  # Also calculate the expected bias relative to the *true* target-group log-OR
  beta_target <- cohort_props$beta_target
  beta_obs <- cohort_props$beta_obs
  bias <- beta_obs - beta_target

  tibble(test_type, sensitivity, specificity, target_group,
         event_rate_multiplier, proportion_multiplier,
         nns_closed_form = reqs$nns, nnr_closed_form = reqs$n_total,
         true_beta = beta_target, predicted_beta = beta_obs, bias_closed_form = bias)
})


# --- 5. OUTPUT ---------------------------------------------------------------

# Ensure output directory exists -------------------------------------------------
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

write_tsv(results_new,
          "results/tables/aim4_closed_form_subgroupE_summary.tsv")

cat("New Aim 4 (Subgroup E) closed-form results saved to results/tables/aim4_closed_form_subgroupE_summary.tsv\n")

results_new %>% gt::gt()
