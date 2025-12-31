## =====================================================================
##  common_parameters.R  -- Single source of shared trial parameters  
##  (Swets et al. CID data + helper grids)                            
## =====================================================================

# This file is loaded (via `source()`) by every analysis / simulation script.
# It defines:
#   * Control-arm baseline risks (with continuity correction where needed)
#   * Subgroup prevalences in the source population
#   * Two OR scenarios (ARREST & Conservative)
#   * Sensitivity / specificity design grid
#   * Convenience helper `get_common_params()` that returns a list
#     so callers can opt for `params <- get_common_params()` and avoid
#     polluting the global namespace.

library(tibble)

## ---------------------------------------------------------------------
##  Raw counts from Swets et al. supplementary table 2
## ---------------------------------------------------------------------
# Control-arm (standard care) cohort
n_ctrl          <- c(A = 60,  B = 52,  C = 138, D = 69,  E = 69)
deaths_ctrl     <- c(A = 13, B =  0, C = 29,  D = 11, E =  3)

# Treatment-arm cohort (for completeness – not used in baseline risks)
n_treat         <- c(A = 55,  B = 55,  C = 138, D = 52,  E = 70)
deaths_treat    <- c(A = 12, B =  8, C = 24,  D = 11, E =  1)

## ---------------------------------------------------------------------
##  Baseline risks (probability of death if UNTREATED)
##  – continuity-corrected for the single zero cell in subgroup B
## ---------------------------------------------------------------------
raw_p0 <- deaths_ctrl / n_ctrl              # 0 for B
apply_zero_cc <- function(p_vec, n_vec, add = 0.5) {
  zero_idx <- which(p_vec == 0)
  if (length(zero_idx)) p_vec[zero_idx] <- add / (n_vec[zero_idx] + add)
  p_vec
}

p0_ctrl <- apply_zero_cc(raw_p0, n_ctrl)    # final baseline-risk vector

## ---------------------------------------------------------------------
##  Overall risks (probability of death across both arms)
##  – used as the default baseline for simulations/closed-form per spec
## ---------------------------------------------------------------------
n_total     <- n_ctrl + n_treat
deaths_total <- deaths_ctrl + deaths_treat
raw_p0_overall <- deaths_total / n_total
p0_overall <- apply_zero_cc(raw_p0_overall, n_total)

## ---------------------------------------------------------------------
##  Subgroup prevalences in the screened population
## ---------------------------------------------------------------------
freq_raw    <- n_ctrl + n_treat
freq_arrest <- freq_raw / sum(freq_raw)     # normalised to 1
freq_ctrl   <- n_ctrl / sum(n_ctrl)

## ---------------------------------------------------------------------
##  Treatment effect scenarios (log-ORs come from Swets paper)
## ---------------------------------------------------------------------
or_arrest       <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
or_conservative <- c(A = 1.0, B =  2.0, C = 0.70, D = 1.2, E = 0.8)
or_arrest_shrunk_k05 <- exp(0.5 * log(or_arrest))

## ---------------------------------------------------------------------
##  Sensitivity / specificity grid used across Aim 3
## ---------------------------------------------------------------------
sens_spec_scenarios <- tribble(
  ~test_type,              ~sensitivity, ~specificity,
  "Perfect (100%)",              1.00,        1.00,
  "Near-Perfect",                0.99,        0.99,
  "High Sens/High Spec",         0.95,        0.95,
  "High Sens/Low Spec",          0.95,        0.70,
  "Low Sens/High Spec",          0.70,        0.95,
  "Balanced/Moderate",           0.80,        0.80
)

target_groups_aim3 <- c("B", "C", "D", "E")

## ---------------------------------------------------------------------
##  Helper accessor (so callers can attach or pull named list)         
## ---------------------------------------------------------------------
get_common_params <- function() {
  list(
    p0_ctrl             = p0_ctrl,
    p0_overall          = p0_overall,
    freq_arrest         = freq_arrest,
    freq_ctrl           = freq_ctrl,
    or_arrest           = or_arrest,
    or_arrest_shrunk_k05 = or_arrest_shrunk_k05,
    or_conservative     = or_conservative,
    sens_spec_scenarios = sens_spec_scenarios,
    target_groups_aim3  = target_groups_aim3
  )
}
