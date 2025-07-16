# -----------------------------------------------------------------------------
# 16_final_clipped_ipw_simulation.R
# -----------------------------------------------------------------------------
# This script simulates and compares four clinical trial designs with a
# fixed final sample size. The fourth design uses a "clipped" probabilistic
# enrollment to improve the stability of the IPW analysis.
# -----------------------------------------------------------------------------

# --- 1. SET-UP -------------------------------------------------------------
pacman::p_load(tidyverse, broom, furrr)
# Use a parallel-safe random number generator for reproducibility
set.seed(2025, kind = "L'Ecuyer-CMRG")

# Detect number of cores allocated by SLURM, with a fallback for local use
slurm_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
n_cores <- if (nzchar(slurm_cores)) as.integer(slurm_cores) else max(1, parallel::detectCores() - 1)
cat("Setting up", n_cores, "parallel workers.\n")
plan(multisession, workers = n_cores)


# --- 2. SCENARIO & USER PARAMETERS ----------------------------------------

# Fixed biological parameters
or_vec   <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
freq_vec <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)
p0_raw <- c(A = 7/60,  B = 0/52,  C = 11/138, D = 11/69, E = 1/69)
p0_raw["B"] <- 0.5 / (52 + 0.5)

# Different discretisations of the score
score_versions <- c("binary", "quartile", "continuous")

# Define the 4 test configurations to run for EACH target group
base_scenarios <- tidyr::crossing(
  threshold_type = c("sens", "spec"),
  snr = c(1, 3)
) %>%
  mutate(
    sens_target = if_else(threshold_type == "sens", 0.95, NA_real_),
    spec_target = if_else(threshold_type == "spec", 0.95, NA_real_)
  )

# Define target groups, N to randomize, and cross with all other parameters
scenario_definitions <- tidyr::crossing(
  group_info = tribble(
    ~target_group, ~n_to_randomize,
    "B",           300,
    "C",           6900,
    "D",           1400,
    "E",           750
  ),
  base_scenarios,
  score_version = score_versions
) %>%
  unnest(group_info) %>%
  select(target_group, n_to_randomize, score_version, everything())
# Simulation settings
n_reps <- 300
alpha  <- 0.05
PILOT_POOL_SIZE <- 1E5
BATCH_SIZE <- 200


# --- 3. HELPER & ANALYSIS FUNCTIONS --------------------------------------

convert_or <- function(or, p0) (or * p0) / (1 - p0 + or * p0)

draw_patient_batch <- function(n, target_group, snr) {
  mu1 <- snr / 2; mu0 <- -snr / 2
  grp <- sample(names(freq_vec), n, TRUE, freq_vec)
  is_target <- grp == target_group
  z <- rnorm(n, mean = ifelse(is_target, mu1, mu0), sd = 1)
  p <- plogis(z)
  tibble(group = grp, is_target = is_target, p = p)
}

compute_tau <- function(pool, threshold_type, sens_target, spec_target) {
  if (threshold_type == "youden") {
    cuts <- sort(unique(pool$p))
    youden_values <- sapply(cuts, function(t) {
      sens <- mean(pool$p[pool$is_target] > t)
      spec <- mean(pool$p[!pool$is_target] <= t)
      sens + spec - 1
    })
    return(cuts[which.max(youden_values)])
  } else if (threshold_type == "sens") {
    return(quantile(pool$p[pool$is_target], probs = 1 - sens_target, type = 1))
  } else if (threshold_type == "spec") {
    return(quantile(pool$p[!pool$is_target], probs = spec_target, type = 1))
  }
}

analyze_trial <- function(cohort, weights_vec, true_beta) {
  n_enrolled <- nrow(cohort)
  if (n_enrolled < 10) return(list(power = NA, bias = NA, mse = NA))

  treat <- rbinom(n_enrolled, 1, 0.5)
  cohort$treat <- treat

  p0 <- p0_raw[cohort$group]
  p1 <- convert_or(or_vec[cohort$group], p0)
  y  <- rbinom(n_enrolled, 1, ifelse(cohort$treat == 1, p1, p0))

  fit <- tryCatch(
    glm(y ~ treat, data = cohort, family = binomial(), weights = weights_vec),
    error = function(e) NULL
  )

  if (is.null(fit) || nrow(coef(summary(fit))) < 2 || anyNA(coef(summary(fit))[2,])) {
    return(list(power = 0, bias = NA, mse = NA))
  }

  coef_row <- coef(summary(fit))[2, ]
  beta_hat <- as.numeric(coef_row["Estimate"])
  pval <- as.numeric(coef_row["Pr(>|z|)"])

  list(power = as.integer(pval < alpha), bias = beta_hat - true_beta, mse = (beta_hat - true_beta)^2)
}


# --- 4. MAIN SIMULATION FUNCTION ------------------------------------------

run_full_scenario <- function(target_group, n_to_randomize, score_version, threshold_type, snr, sens_target, spec_target) {
  
  cat("Running Scenario: Group", target_group, "| N_to_randomize", n_to_randomize,
      "| Score", score_version, "| SNR", snr, "| Type", threshold_type, "\n")

  true_beta <- log(or_vec[target_group])
  
  run_one_rep <- function(rep_id) {
    pilot_pool <- draw_patient_batch(PILOT_POOL_SIZE, target_group, snr)
    tau <- compute_tau(pilot_pool, threshold_type, sens_target, spec_target)
    
    sens_emp <- mean(pilot_pool$p[pilot_pool$is_target] > tau)
    spec_emp <- mean(pilot_pool$p[!pilot_pool$is_target] <= tau)

    # Recruitment for Hard/Weighted designs
    n_screened_hard <- 0
    enrolled_hard <- list()
    while(nrow(bind_rows(enrolled_hard)) < n_to_randomize) {
      new_batch <- draw_patient_batch(BATCH_SIZE, target_group, snr)
      n_screened_hard <- n_screened_hard + BATCH_SIZE
      eligible <- new_batch %>% filter(p > tau)
      if(nrow(eligible) > 0) enrolled_hard[[length(enrolled_hard) + 1]] <- eligible
      if(n_screened_hard > 10 * PILOT_POOL_SIZE) break
    }
    cohort_hard <- bind_rows(enrolled_hard) %>% head(n_to_randomize)

    # Recruitment for Probabilistic design
    n_screened_prob <- 0
    enrolled_prob <- list()
    while(nrow(bind_rows(enrolled_prob)) < n_to_randomize) {
      new_batch <- draw_patient_batch(BATCH_SIZE, target_group, snr)
      n_screened_prob <- n_screened_prob + BATCH_SIZE
      
      p_enroll <- switch(score_version,
          continuous = new_batch$p,
          quartile   = ntile(new_batch$p, 4)/4,
          binary     = as.numeric(new_batch$p > tau))
      
      # *** CLIPPING LOGIC ***
      # Exclude anyone with less than a 10% chance of enrollment
      p_enroll[p_enroll < 0.1] <- 0
      
      new_batch$p_enroll <- p_enroll
      eligible_indices <- rbinom(nrow(new_batch), 1, new_batch$p_enroll) == 1
      if(sum(eligible_indices) > 0) enrolled_prob[[length(enrolled_prob) + 1]] <- new_batch[eligible_indices, ]
      if(n_screened_prob > 10 * PILOT_POOL_SIZE) break
    }
    cohort_prob <- bind_rows(enrolled_prob) %>% head(n_to_randomize)

    # --- Analysis Phase ---
    weights_for_B <- switch(score_version,
      continuous = cohort_hard$p,
      quartile   = ntile(cohort_hard$p, 4)/4,
      binary     = NULL) # Corrected to NULL for unweighted analysis
      
    weights_for_D <- 1 / ifelse(cohort_prob$p_enroll > 0, cohort_prob$p_enroll, 1e-9)
    
    res_A_hard_unweighted  <- analyze_trial(cohort_hard, weights_vec = NULL, true_beta)
    res_B_soft_weighted    <- analyze_trial(cohort_hard, weights_vec = weights_for_B, true_beta)
    res_C_prob_unweighted  <- analyze_trial(cohort_prob, weights_vec = NULL, true_beta)
    res_D_prob_ipw         <- analyze_trial(cohort_prob, weights_vec = weights_for_D, true_beta)

    # --- Collate Results ---
    tibble(
      sens_emp = sens_emp, spec_emp = spec_emp, nns_hard = n_screened_hard, nns_prob = n_screened_prob,
      power_hard = res_A_hard_unweighted$power, power_soft_weighted = res_B_soft_weighted$power,
      power_soft_prob_unw = res_C_prob_unweighted$power, power_soft_prob_ipw = res_D_prob_ipw$power,
      bias_hard = res_A_hard_unweighted$bias, bias_soft_weighted = res_B_soft_weighted$bias,
      bias_soft_prob_unw = res_C_prob_unweighted$bias, bias_soft_prob_ipw = res_D_prob_ipw$bias,
      mse_hard = res_A_hard_unweighted$mse, mse_soft_weighted = res_B_soft_weighted$mse,
      mse_soft_prob_unw = res_C_prob_unweighted$mse, mse_soft_prob_ipw = res_D_prob_ipw$mse
    )
  }

  reps_results <- future_map_dfr(1:n_reps, run_one_rep, .options = furrr_options(seed = TRUE))
  reps_results %>% summarise(across(everything(), ~mean(.x, na.rm = TRUE)))
}

# --- 5. MAIN EXECUTION LOOP -----------------------------------------------
full_results <- pmap_dfr(scenario_definitions, run_full_scenario)

# Bind scenario definitions to the results for a full table
final_output <- bind_cols(scenario_definitions, full_results)


# --- 6. SAVE AND PRINT ----------------------------------------------------
print(final_output, n = nrow(final_output), width = Inf)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
file_out <- "results/tables/aim4_final_clipped_ipw.tsv"
write_tsv(final_output, file_out)
cat("\nSaved to", file_out, "\n")