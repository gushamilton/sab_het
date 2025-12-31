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
source("code/common_parameters.R")
scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}
output_root <- file.path("results", scenario_set)
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)

# --- 2. SCENARIOS & INPUT DATA ----------------------------------------------

params <- get_common_params()

# Odds-ratio scenarios (same as in 01_run_simulations.R)
or_arrest_raw    <- params$or_arrest
or_arrest_shrunk <- params$or_arrest_shrunk_k05
or_conservative  <- params$or_conservative

scenario_definitions <- if (scenario_set == "main") {
  tribble(
    ~scenario_name,   ~or_vector,
    "ARREST_raw",     or_arrest_raw
  )
} else {
  tribble(
    ~scenario_name,   ~or_vector,
    "ARREST_shrunk",  or_arrest_shrunk
  )
}

# Prevalence and baseline risks (overall mortality across both arms)
freq_arrest        <- params$freq_arrest
p0_arrest_adjusted <- params$p0_overall

# Sensitivity/specificity test grid (including the “Perfect” test)
sens_spec_scenarios <- params$sens_spec_scenarios
target_groups <- params$target_groups_aim3

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

  # --- PREP: Calculate population-level parameters ---
  p_vec <- freq_vector / sum(freq_vector) # Normalize prevalences
  p_target <- p_vec[target_group]         # Prevalence of the target subgroup
  p_non_target <- 1 - p_target            # Prevalence of all other (contaminant) subgroups

  # Define which groups are contaminants
  non_target_groups <- names(p_vec)[names(p_vec) != target_group]

  # --- 1. Calculate properties of the two conceptual components: TRUE TARGET vs. CONTAMINANTS ---
  
  # A. Properties of the TRUE TARGET subgroup
  beta_target <- log(or_vector[target_group]) # True log(OR) in the target group
  p0_target <- p0_vector[target_group]        # True baseline risk in the target group
  p1_target <- get_or_p1(or_vector[target_group], p0_target) # True risk under treatment in the target group

  # B. Properties of the CONTAMINANT MIX (a weighted average of all non-target groups)
  if (p_non_target > 0) {
    # Weight each non-target group by its prevalence *within the contaminant pool*
    weights_mix <- p_vec[non_target_groups] / p_non_target

    # Calculate weighted average baseline risk (p0) and treated risk (p1) for the contaminant mix
    p0_mix   <- sum(p0_vector[non_target_groups] * weights_mix)
    p1_mix   <- sum(get_or_p1(or_vector[non_target_groups], p0_vector[non_target_groups]) * weights_mix)
  } else {
    # Handle the edge case where the target group is 100% prevalent (no contaminants)
    p0_mix   <- NA_real_
    p1_mix   <- NA_real_
  }

  # --- 2. Calculate properties of the final ENROLLED (mixed) cohort ---
  # The enrolled cohort consists of true positives and false positives.

  # Calculate the proportion of the *entire screened population* that are true positives (TP) and false positives (FP)
  tp_rate <- sensitivity * p_target     # Correctly identified targets
  fp_rate <- (1 - specificity) * p_non_target # Incorrectly identified non-targets
  
  # The overall enrollment rate is the sum of TP and FP rates
  enrol_rate <- tp_rate + fp_rate

  # Handle the case of a perfect test that perfectly excludes everyone (enroll rate is 0)
  if (enrol_rate == 0) {
    return(list(beta_obs = 0, p0_obs = 0, p1_obs = 0, enrol_rate = 0,
                p_true_in_enrolled = NA_real_,
                beta_target = beta_target))
  }

  # Calculate the composition of the *enrolled cohort*
  p_true_in_enrolled <- tp_rate / enrol_rate # Proportion of enrolled patients who are true targets
  p_mix_in_enrolled  <- fp_rate / enrol_rate # Proportion of enrolled patients who are contaminants

  # Calculate the OBSERVED (diluted) event rates in the final enrolled cohort
  # This is a weighted average of the rates from the true targets and the contaminant mix
  if (p_non_target > 0) {
    p0_obs   <- p_true_in_enrolled * p0_target + p_mix_in_enrolled * p0_mix # Observed baseline risk
    p1_obs   <- p_true_in_enrolled * p1_target + p_mix_in_enrolled * p1_mix # Observed treated risk
  } else {
    # If no contaminants, the observed rates are just the true target rates
    p0_obs <- p0_target
    p1_obs <- p1_target
  }

  # FIX: derive observed log-OR from pooled risks so that
  # beta_obs == log( (p1_obs/(1-p1_obs)) / (p0_obs/(1-p0_obs)) ).
  # This is the diluted effect size we expect to observe in the trial.
  if (p0_obs %in% c(0,1) || p1_obs %in% c(0,1)) {
    beta_obs <- 0  # Handle edge cases where OR is undefined; will be trapped by next function
  } else {
    beta_obs <- log( (p1_obs / (1 - p1_obs)) / (p0_obs / (1 - p0_obs)) )
  }

  # Return all calculated properties of the enrolled cohort
  list(beta_obs = beta_obs,           # The OBSERVED log(OR) in the enrolled cohort
       p0_obs = p0_obs,               # The OBSERVED baseline risk
       p1_obs = p1_obs,               # The OBSERVED risk under treatment
       enrol_rate = enrol_rate,       # The proportion of the screened population that is enrolled
       p_true_in_enrolled = p_true_in_enrolled, # The purity of the enrolled cohort
       beta_target = beta_target)     # The TRUE log(OR) for the target group (for bias calculation)
}


calc_required_nns <- function(beta_obs, p0_obs, p1_obs, enrol_rate,
                              z_alpha_half, z_power_target) {

  # --- INPUT GUARD: Check for degenerate cases ---
  # If the observed effect size is zero, or rates are 0/1, power is undefined, so N is infinite.
  if (abs(beta_obs) < 1e-6 || is.na(beta_obs) || enrol_rate == 0 ||
      p0_obs == 0 || p0_obs == 1 || p1_obs == 0 || p1_obs == 1) {
    return(list(n_total = Inf, nns = Inf))
  }

  # --- SAMPLE SIZE CALCULATION ---
  # This uses the standard formula for comparing two proportions (via log-OR).
  # The variance of the log(OR) is approximated by the sum of the inverse variances of the two binomials.
  
  # 1. Calculate the variance term for the log(OR) based on the OBSERVED event rates in the mixed cohort.
  # This term is 1/p0(1-p0) + 1/p1(1-p1).
  var_term <- (1 / (p0_obs * (1 - p0_obs))) + (1 / (p1_obs * (1 - p1_obs)))
  if (!is.finite(var_term)) return(list(n_total = Inf, nns = Inf))

  # 2. Apply the standard sample size formula for a two-group comparison with 1:1 allocation.
  # N_total = 2 * (Z_alpha/2 + Z_power)^2 * (variance_term) / (log(OR))^2
  # The factor of 2 at the start accounts for the two groups (treatment and control).
  n_total <- 2 * (z_alpha_half + z_power_target)^2 * var_term / (beta_obs^2)  # FIX

  # 3. Calculate Number Needed to Screen (NNS).
  # This is the total required sample size (n_total, or NNR) divided by the proportion of
  # screened patients who are actually enrolled.
  nns <- ceiling(n_total / enrol_rate)

  list(n_total = ceiling(n_total), nns = nns)
}

# --- 4. MAIN CALCULATION -----------------------------------------------------

all_cases <- expand_grid(scenario_definitions,
                         sens_spec_scenarios,
                         target_group = target_groups)

results_closed <- pmap_dfr(all_cases, function(scenario_name, or_vector,
                                               test_type, sensitivity, specificity,
                                               target_group) {

  # Calculate diluted effect size AND diluted event rates for the enrolled cohort
  cohort_props <- calculate_mixed_cohort_properties(target_group, sensitivity, specificity,
                                                    or_vector, p0_arrest_adjusted, freq_arrest)

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

  tibble(scenario_name, test_type, sensitivity, specificity, target_group,
         nns_closed_form = reqs$nns, nnr_closed_form = reqs$n_total,
         true_beta = beta_target, predicted_beta = beta_obs, bias_closed_form = bias)
})

# --- 5. OUTPUT ---------------------------------------------------------------

results_closed <- results_closed %>% mutate(scenario_set = scenario_set)
write_tsv(
  results_closed,
  file.path(output_root, "tables/aim3_closed_form_summary.tsv")
)

cat("Closed-form Aim 3 results saved to", file.path(output_root, "tables/aim3_closed_form_summary.tsv"), "\n")

if (Sys.getenv("PRINT_GT", unset = "1") == "1") {
  results_closed %>% gt::gt()
}
