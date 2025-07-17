# -----------------------------------------------------------------------------
# 09_compare_trial_designs.R
# -----------------------------------------------------------------------------
# This script simulates and compares three different clinical trial designs
# that use a biomarker score for patient selection and analysis.
#
# The three designs are:
#   1) Hard Enrichment: Enroll only patients with a score above a sharp
#      threshold. The analysis is unweighted.
#   2) Soft Weighted Analysis: Enroll the same group as the Hard design, but
#      use the continuous score as a weight in the final statistical analysis.
#   3) Soft Probabilistic Enrollment: Probabilistically enroll patients based
#      on their score, using an adaptive scaling method to match the sample
#      size of the Hard design. The final analysis is unweighted.
# -----------------------------------------------------------------------------

# --- 1. SET-UP -------------------------------------------------------------
pacman::p_load(tidyverse, broom, furrr)
# Use a parallel-safe random number generator for reproducibility
set.seed(2025, kind = "L'Ecuyer-CMRG")

# Detect number of cores allocated by SLURM, with a fallback for local use
slurm_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
n_cores <- if (nzchar(slurm_cores)) {
  as.integer(slurm_cores)
} else {
  max(1, parallel::detectCores() - 1)
}
cat("Setting up", n_cores, "parallel workers.\n")
plan(multisession, workers = n_cores)


# --- 2. SCENARIO & USER PARAMETERS ----------------------------------------

# Fixed biological parameters
or_vec   <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
freq_vec <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)
p0_raw <- c(A = 7/60,  B = 0/52,  C = 11/138, D = 11/69, E = 1/69)
p0_raw["B"] <- 0.5 / (52 + 0.5)  # Continuity correction

# Define the 4 test configurations to run for EACH target group
base_scenarios <- tidyr::crossing(
  threshold_type = c("sens", "spec"),
  snr = c(1, 3)
) %>%
  mutate(
    sens_target = if_else(threshold_type == "sens", 0.95, NA_real_),
    spec_target = if_else(threshold_type == "spec", 0.95, NA_real_)
  )

# Define target groups and their sample sizes
scenario_definitions <- tidyr::crossing(
  group_info = tribble(
    ~target_group, ~n_screened,
    "B",           1500,
    "C",           15000,
    "D",           6000,
    "E",           9000
  ),
  base_scenarios
) %>%
  unnest(group_info) %>%
  select(target_group, n_screened, everything())


# Different discretisations of the score available for analysis
score_versions <- c("binary", "tertile", "quartile", "decile", "continuous")

# Simulation settings
n_reps <- 300
alpha  <- 0.05


# --- 3. Helper Functions ---------------------------------------------------

convert_or <- function(or, p0) {
  (or * p0) / (1 - p0 + or * p0)
}

# Determine threshold tau for a given score vector
compute_tau <- function(p, is_target, threshold_type, sens_target, spec_target) {
  if (threshold_type == "sens") {
    if (is.na(sens_target)) stop("sens_target must be provided for threshold_type 'sens'")
    return(quantile(p[is_target], probs = 1 - sens_target, type = 1))
  } else if (threshold_type == "spec") {
    if (is.na(spec_target)) stop("spec_target must be provided for threshold_type 'spec'")
    return(quantile(p[!is_target], probs = 1 - spec_target, type = 1))
  } else {
    stop("Unknown threshold_type")
  }
}

# Generate one pool of simulated patients
make_pool <- function(snr, score_version, target_group, n_screened,
                      threshold_type, sens_target, spec_target) {
  mu1 <-  snr/2
  mu0 <- -snr/2
  grp <- sample(names(freq_vec), n_screened, TRUE, freq_vec)
  is_target <- grp == target_group
  z <- rnorm(n_screened, mean = ifelse(is_target, mu1, mu0), sd = 1)
  p <- plogis(z)

  tau <- compute_tau(p, is_target, threshold_type, sens_target, spec_target)
  test_positive <- as.integer(p > tau)

  p_exposed <- switch(score_version,
    continuous = p,
    decile     = ntile(p, 10)/10,
    quartile   = ntile(p, 4)/4,
    tertile    = ntile(p, 3)/3,
    binary     = test_positive,
    stop("bad score_version"))

  tibble(group = grp, is_target, test_positive, p_target = p_exposed, p, tau)
}


# --- 4. Trial Design Functions ----------------------------------------------

## Design 1: Hard Enrichment
hard_enrichment <- function(pool) {
  enrol <- pool %>% filter(test_positive == 1)
  n_enrolled <- nrow(enrol)
  if (n_enrolled < 10) {
    return(list(pval = 1, n_enrolled = n_enrolled, ess = n_enrolled, beta_hat = NA))
  }
  treat <- rbinom(n_enrolled, 1, 0.5)
  p0 <- p0_raw[enrol$group]
  p1 <- convert_or(or_vec[enrol$group], p0)
  y  <- rbinom(n_enrolled, 1, ifelse(treat == 1, p1, p0))

  fit <- tryCatch(glm(y ~ treat, family = binomial()), error = function(e) NULL)

  if (is.null(fit) || nrow(coef(summary(fit))) < 2 || anyNA(coef(summary(fit))[2,])) {
    return(list(pval = 1, n_enrolled = n_enrolled, ess = n_enrolled, beta_hat = NA))
  }

  coef_row <- coef(summary(fit))[2, ]
  list(
    pval = as.numeric(coef_row["Pr(>|z|)"]),
    n_enrolled = n_enrolled,
    ess = n_enrolled,
    beta_hat = as.numeric(coef_row["Estimate"])
  )
}


## Design 2: Soft Weighted Analysis
soft_weighted_analysis <- function(pool) {
  enrol <- pool %>% filter(test_positive == 1)
  n_enrolled <- nrow(enrol)
  if (n_enrolled < 10) {
    return(list(pval = 1, n_enrolled = n_enrolled, ess = NA, beta_hat = NA))
  }
  treat <- rbinom(n_enrolled, 1, 0.5)
  enrol_data <- enrol %>% mutate(treat = treat)

  p0 <- p0_raw[enrol_data$group]
  p1 <- convert_or(or_vec[enrol_data$group], p0)
  y  <- rbinom(n_enrolled, 1, ifelse(enrol_data$treat == 1, p1, p0))

  w <- enrol_data$p_target
  ess <- sum(w)

  fit <- tryCatch(
    glm(y ~ treat, data = enrol_data, family = binomial(), weights = w),
    error = function(e) NULL
  )

  if (is.null(fit) || nrow(coef(summary(fit))) < 2 || anyNA(coef(summary(fit))[2,])) {
    return(list(pval = 1, n_enrolled = n_enrolled, ess = ess, beta_hat = NA))
  }

  coef_row <- coef(summary(fit))[2, ]
  list(
    pval = as.numeric(coef_row["Pr(>|z|)"]),
    n_enrolled = n_enrolled,
    ess = ess,
    beta_hat = as.numeric(coef_row["Estimate"])
   )
}


## Design 3: Soft Probabilistic Enrollment (with Adaptive Scaling)
soft_probabilistic_enrollment <- function(pool) {
  n_target <- sum(pool$test_positive)

  if (n_target < 10) {
      return(list(pval = 1, n_enrolled = 0, ess = 0, beta_hat = NA))
  }
  
  sum_of_scores <- sum(pool$p)
  k <- if (sum_of_scores > 0) n_target / sum_of_scores else 0
  enroll_prob <- pmin(1, k * pool$p)

  enroll_indices <- rbinom(nrow(pool), 1, enroll_prob) == 1
  enrol <- pool[enroll_indices, ]
  n_enrolled <- nrow(enrol)

  if (n_enrolled < 10) {
    return(list(pval = 1, n_enrolled = n_enrolled, ess = n_enrolled, beta_hat = NA))
  }

  treat <- rbinom(n_enrolled, 1, 0.5)
  enrol_data <- enrol %>% mutate(treat = treat)

  p0 <- p0_raw[enrol_data$group]
  p1 <- convert_or(or_vec[enrol_data$group], p0)
  y  <- rbinom(n_enrolled, 1, ifelse(enrol_data$treat == 1, p1, p0))

  fit <- tryCatch(
    glm(y ~ treat, data = enrol_data, family = binomial()),
    error = function(e) NULL
  )

  if (is.null(fit) || nrow(coef(summary(fit))) < 2 || anyNA(coef(summary(fit))[2,])) {
    return(list(pval = 1, n_enrolled = n_enrolled, ess = n_enrolled, beta_hat = NA))
  }

  coef_row <- coef(summary(fit))[2, ]
  list(
    pval = as.numeric(coef_row["Pr(>|z|)"]),
    n_enrolled = n_enrolled,
    ess = n_enrolled,
    beta_hat = as.numeric(coef_row["Estimate"])
  )
}


# --- 5. Simulation Wrappers ------------------------------------------------

# This function runs n_reps for ONE specific scenario and ONE score_version
run_reps_for_combo <- function(score_ver, target_group, n_screened,
                               threshold_type, sens_target, spec_target, snr) {

  true_beta <- log(or_vec[target_group])

  reps <- future_map_dfr(1:n_reps, function(i) {
    pool <- make_pool(snr, score_ver, target_group, n_screened,
                      threshold_type, sens_target, spec_target)

    res_hard <- hard_enrichment(pool)
    res_weighted <- soft_weighted_analysis(pool)
    res_prob <- soft_probabilistic_enrollment(pool)

    tibble(
      power_hard = res_hard$pval < alpha,
      power_weighted = res_weighted$pval < alpha,
      power_prob = res_prob$pval < alpha,
      ess_hard = res_hard$ess,
      ess_weighted = res_weighted$ess,
      ess_prob = res_prob$ess,
      beta_hard = res_hard$beta_hat,
      beta_weighted = res_weighted$beta_hat,
      beta_prob = res_prob$beta_hat,
      sens_emp = mean(pool$test_positive[pool$is_target]),
      spec_emp = mean(1 - pool$test_positive[!pool$is_target])
    )
  }, .options = furrr_options(seed = TRUE))

  # Summarize results for all three designs
  summarise(reps,
            power_hard = mean(power_hard),
            power_soft_weighted = mean(power_weighted),
            power_soft_prob_enroll = mean(power_prob),
            mean_ess_hard = mean(ess_hard, na.rm = TRUE),
            mean_ess_soft_weighted = mean(ess_weighted, na.rm = TRUE),
            mean_ess_soft_prob_enroll = mean(ess_prob, na.rm = TRUE),
            bias_hard = mean(beta_hard - true_beta, na.rm = TRUE),
            bias_soft_weighted = mean(beta_weighted - true_beta, na.rm = TRUE),
            bias_soft_prob_enroll = mean(beta_prob - true_beta, na.rm = TRUE),
            mse_hard  = mean((beta_hard - true_beta)^2, na.rm = TRUE),
            mse_soft_weighted  = mean((beta_weighted - true_beta)^2, na.rm = TRUE),
            mse_soft_prob_enroll  = mean((beta_prob - true_beta)^2, na.rm = TRUE),
            sens_emp = mean(sens_emp, na.rm = TRUE),
            spec_emp = mean(spec_emp, na.rm = TRUE)
  )
}

# This function takes one row from the scenario table and loops over score_versions
run_one_scenario <- function(target_group, n_screened, threshold_type,
                             sens_target, spec_target, snr) {
  cat("Running scenario: Group", target_group, "| SNR", snr, "| Type", threshold_type, "\n")
  purrr::map_dfr(score_versions, function(ver) {
    results <- run_reps_for_combo(
      ver, target_group, n_screened,
      threshold_type, sens_target, spec_target, snr
    )
    # Add scenario parameters to the results row
    mutate(results,
           score_version = ver,
           .before = 1
    )
  })
}


# --- 6. MAIN LOOP ----------------------------------------------------------
# Use pmap to iterate over the rows of the scenario_definitions tibble
results <- pmap_dfr(scenario_definitions, .f = run_one_scenario)


# --- 7. SAVE and PRINT -----------------------------------------------------
print(results, n = nrow(results))
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
file_out <- "results/tables/aim4_three_designs_comparison.tsv"
write_tsv(results, file_out)
cat("\nSaved to", file_out, "\n")