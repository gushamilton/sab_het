# -----------------------------------------------------------------------------
# 02_closed_form_aim3.R
#
# Analytic (closed-form) approximation for Aim 3: compute the minimum number of
# patients to screen (NNS) required to achieve 80 % power for an enrichment
# trial. This version uses a more accurate variance calculation based on the
# expected event rates within the *enrolled (mixed)* cohort.
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
# Correct for zero events in Group B
p0_B_corrected     <- 0.5 / (52 + 0.5)
p0_arrest_adjusted <- p0_arrest_raw
p0_arrest_adjusted["B"] <- p0_B_corrected


# Sensitivity/specificity test grid (including the “Perfect” test)
sens_spec_scenarios <- tribble(
  ~test_type,            ~sensitivity, ~specificity,
  "Perfect",            1.00,         1.00,
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

  # Contaminant Mix
  # Weight each non-target group by its prevalence within the non-target pool
  weights_mix <- p_vec[non_target_groups] / p_non_target
  
  # Weighted average properties of the contaminant mix
  p0_mix   <- sum(p0_vector[non_target_groups] * weights_mix)
  beta_mix <- sum(log(or_vector[non_target_groups]) * weights_mix)
  p1_mix   <- sum(get_or_p1(or_vector[non_target_groups], p0_vector[non_target_groups]) * weights_mix)

  # --- 2. Calculate properties of the final enrolled (mixed) cohort
  # Proportion of enrolled cohort that is true target vs. contaminant
  tp_rate <- sensitivity * p_target
  fp_rate <- (1 - specificity) * p_non_target
  enrol_rate <- tp_rate + fp_rate
  
  if (enrol_rate == 0) return(list(beta_obs = 0, p0_obs = 0, p1_obs = 0, enrol_rate = 0))

  p_true_in_enrolled <- tp_rate / enrol_rate
  p_mix_in_enrolled  <- fp_rate / enrol_rate

  # Observed (diluted) properties in the final enrolled cohort
  beta_obs <- p_true_in_enrolled * beta_target + p_mix_in_enrolled * beta_mix
  p0_obs   <- p_true_in_enrolled * p0_target   + p_mix_in_enrolled * p0_mix
  p1_obs   <- p_true_in_enrolled * p1_target   + p_mix_in_enrolled * p1_mix

  list(beta_obs = beta_obs, p0_obs = p0_obs, p1_obs = p1_obs, enrol_rate = enrol_rate)
}


calc_required_nns <- function(beta_obs, p0_obs, p1_obs, enrol_rate,
                              z_alpha_half, z_power_target) {

  if (abs(beta_obs) < 1e-6 || is.na(beta_obs) || enrol_rate == 0 ||
      p0_obs %in% c(0, 1) || p1_obs %in% c(0, 1)) {
    return(list(n_total = Inf, nns = Inf))
  }

  # Variance for log-OR based on the OBSERVED event rates in the mixed, enrolled cohort
  var_term <- (1 / (p0_obs * (1 - p0_obs))) + (1 / (p1_obs * (1 - p1_obs)))
  if (is.infinite(var_term)) return(list(n_total = Inf, nns = Inf))

  # Standard sample size formula, assuming balanced arms (n_total / 2 per arm)
  n_total <- (z_alpha_half + z_power_target)^2 * var_term / (beta_obs^2)
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

  tibble(scenario_name, test_type, sensitivity, specificity, target_group,
         nns_closed_form = reqs$nns, nnr_closed_form = reqs$n_total)
})

# --- 5. OUTPUT ---------------------------------------------------------------

write_tsv(results_closed,
          "results/tables/aim3_closed_form_summary.tsv")

cat("Closed-form Aim 3 results saved to results/tables/aim3_closed_form_summary.tsv\n")

results_closed
