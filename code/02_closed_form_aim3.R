# -----------------------------------------------------------------------------
# 02_closed_form_aim3.R
#
# Analytic (closed-form) approximation for Aim 3: compute the minimum number of
# patients to screen (NNS) required to achieve 80 % power for an enrichment
# trial, without Monte-Carlo simulation. This version uses a more accurate
# variance calculation based on event rates, not the simplified "4/N" formula.
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

calc_observed_beta <- function(target_group, sensitivity, specificity,
                               or_vector, freq_vector) {
  # Prevalence of each group
  p_vec <- freq_vector / sum(freq_vector)
  p_target <- p_vec[target_group]
  p_non_target <- 1 - p_target

  # True-/false-positive rates per screened patient
  tp_rate <- sensitivity * p_target
  fp_rate <- (1 - specificity) * p_non_target
  enrol_rate <- tp_rate + fp_rate      # Pr(test +) per screen

  # Proportion of enrolled who are true target
  p_true <- tp_rate / enrol_rate

  # Average log-OR among contaminating patients (weighted by prevalence)
  log_or_vec <- log(or_vector)
  avg_beta_mix <- sum(log_or_vec[names(p_vec) != target_group] *
                      (p_vec[names(p_vec) != target_group] / p_non_target))

  beta_target <- log_or_vec[target_group]
  beta_obs    <- p_true * beta_target + (1 - p_true) * avg_beta_mix

  list(beta_obs = beta_obs,
       enrol_rate = enrol_rate)
}

calc_required_nns <- function(beta_obs, enrol_rate,
                              p0_target, or_target,
                              z_alpha_half, z_power_target) {
  
  if (abs(beta_obs) < 1e-6 || is.na(beta_obs) || enrol_rate == 0) {
    return(list(n_total = Inf, nns = Inf))
  }

  # Calculate event rate in the treatment arm for the *true* target subgroup
  p1_target <- (or_target * p0_target) / (1 - p0_target + (or_target * p0_target))
  
  # Check for p=0 or p=1 which would make variance infinite
  if (p0_target %in% c(0, 1) || p1_target %in% c(0, 1)) {
     return(list(n_total = Inf, nns = Inf))
  }

  # More accurate variance formula for log-OR from a cohort study/RCT.
  # Var(logOR) = 1/(N_t*p_t) + 1/(N_t*(1-p_t)) + 1/(N_c*p_c) + 1/(N_c*(1-p_c))
  # This uses event rates from the true target group as a proxy for the enrolled group's variance.
  var_term <- (1 / (p0_target * (1 - p0_target))) + (1 / (p1_target * (1 - p1_target)))
  
  if (is.infinite(var_term)) return(list(n_total = Inf, nns = Inf))
  
  # Assuming balanced arms (n_total / 2 per arm)
  # Rearranging the standard sample size formula for log-OR to solve for n_total:
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
  # Calculate diluted effect size and enrolment rate
  vals <- calc_observed_beta(target_group, sensitivity, specificity,
                             or_vector, freq_arrest)
  beta_obs   <- vals$beta_obs
  enrol_rate <- vals$enrol_rate

  # Calculate required NNS/NNR using the more accurate variance formula
  reqs <- calc_required_nns(beta_obs, enrol_rate,
                            p0_target = p0_arrest_adjusted[target_group],
                            or_target = or_vector[target_group],
                            z_alpha_half, z_power_target)

  tibble(scenario_name, test_type, sensitivity, specificity, target_group,
         nns_closed_form = reqs$nns, nnr_closed_form = reqs$n_total)
})

# --- 5. OUTPUT ---------------------------------------------------------------

write_tsv(results_closed,
          "results/tables/aim3_closed_form_summary.tsv")

cat("Closed-form Aim 3 results saved to results/tables/aim3_closed_form_summary.tsv\n") 

results_closed
