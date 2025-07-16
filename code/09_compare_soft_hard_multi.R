# -----------------------------------------------------------------------------
# 08_compare_soft_hard.R (updated – multi-scenario version)
# -----------------------------------------------------------------------------
# This script runs a series of pre-defined simulation scenarios to compare
# the power and bias of HARD vs. SOFT enrichment designs.
#
# Each scenario defines a target subgroup, a sample size, a latent predictor
# strength (SNR), and a thresholding policy for the diagnostic test.
#
# The main loop iterates through these scenarios, and for each one, it
# evaluates the performance across multiple discretisations of the SOFT score
# (from binary to continuous).
# -----------------------------------------------------------------------------

# --- 1. SET-UP -------------------------------------------------------------
pacman::p_load(tidyverse, broom, furrr)
# Use a parallel-safe random number generator for reproducibility
set.seed(2025, kind = "L'Ecuyer-CMRG")
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

# Define the 4 scenarios to run for EACH target group
base_scenarios <- tidyr::crossing(
  threshold_type = c("sens", "spec"),
  snr = c(1, 3)
) %>%
  mutate(
    sens_target = if_else(threshold_type == "sens", 0.9, NA_real_),
    spec_target = if_else(threshold_type == "spec", 0.9, NA_real_)
  )

# Define target groups and their estimated sample size for ~60% power
scenario_definitions <- tidyr::crossing(
  group_info = tribble(
    ~target_group, ~n_screened,
    "B",           400,
    "C",           11500,
    "D",           4000,
    "E",           2200
  ),
  base_scenarios
) %>%
  unnest(group_info) %>%
  select(target_group, n_screened, everything())


# Different discretisations of the score available to the SOFT design
score_versions <- c("binary", "tertile", "quartile", "decile", "continuous")

# Simulation settings
n_reps <- 50
alpha  <- 0.05


# --- 3. Helper Functions ---------------------------------------------------

convert_or <- function(or, p0) {
  (or * p0) / (1 - p0 + or * p0)
}

# Determine threshold tau for a given score vector (now takes arguments)
compute_tau <- function(p, is_target, threshold_type, sens_target, spec_target) {
  if (threshold_type == "youden") {
    cuts <- sort(unique(p))
    youden <- sapply(cuts, function(t) {
      se <- mean(p[is_target] > t)
      sp <- mean(p[!is_target] <= t)
      se + sp - 1
    })
    return(cuts[which.max(youden)])
  } else if (threshold_type == "sens") {
    if (is.na(sens_target)) stop("sens_target must be provided for threshold_type 'sens'")
    # Use quantile for robustness with step functions
    return(quantile(p[is_target], probs = 1 - sens_target, type = 1))
  } else if (threshold_type == "spec") {
    if (is.na(spec_target)) stop("spec_target must be provided for threshold_type 'spec'")
    # Use quantile for robustness with step functions
    return(quantile(p[!is_target], probs = 1 - spec_target, type = 1))
  } else {
    stop("Unknown threshold_type")
  }
}

# Generate one pool (now takes arguments)
make_pool <- function(snr, score_version, target_group, n_screened,
                      threshold_type, sens_target, spec_target) {
  # Means chosen so that mu1 - mu0 = snr, centred at 0
  mu1 <-  snr/2
  mu0 <- -snr/2

  grp <- sample(names(freq_vec), n_screened, TRUE, freq_vec)
  is_target <- grp == target_group
  z <- rnorm(n_screened, mean = ifelse(is_target, mu1, mu0), sd = 1)
  p <- plogis(z)

  tau <- compute_tau(p, is_target, threshold_type, sens_target, spec_target)
  test_positive <- as.integer(p > tau)

  # Expose chosen discretisation
  p_exposed <- switch(score_version,
    continuous = p,
    decile     = ntile(p, 10)/10,
    quartile   = ntile(p, 4)/4,
    tertile    = ntile(p, 3)/3,
    binary     = test_positive,
    stop("bad score_version"))

  # Return pool
  tibble(group = grp, is_target, test_positive, p_target = p_exposed, p, tau)
}

# HARD analysis (unchanged)
hard_design <- function(pool) {
  enrol <- pool %>% filter(test_positive == 1)
  n_enrolled <- nrow(enrol)
  if (n_enrolled < 10) { # Not enough patients to run a trial
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

# SOFT analysis (unchanged)
soft_design <- function(pool) {
  treat_prob <- 0.5 * pool$p_target
  treat <- rbinom(nrow(pool), 1, treat_prob)
  w <- pool$p_target
  ess <- sum(w)

  if (ess < 10) { # Not enough effective sample size
    return(list(pval = 1, n_enrolled = sum(pool$p_target > 0), ess = ess, beta_hat = NA))
  }

  p0 <- p0_raw[pool$group]
  p1 <- convert_or(or_vec[pool$group], p0)
  y  <- rbinom(nrow(pool), 1, ifelse(treat == 1, p1, p0))

  fit <- tryCatch(glm(y ~ treat, family = binomial(), weights = w), error = function(e) NULL)

  if (is.null(fit) || nrow(coef(summary(fit))) < 2 || anyNA(coef(summary(fit))[2,])) {
    return(list(pval = 1, n_enrolled = sum(pool$p_target > 0), ess = ess, beta_hat = NA))
  }

  coef_row <- coef(summary(fit))[2, ]
  list(
    pval = as.numeric(coef_row["Pr(>|z|)"]),
    n_enrolled = sum(pool$p_target > 0),
    ess = ess,
    beta_hat = as.numeric(coef_row["Estimate"])
   )
}


# This function runs n_reps for ONE specific scenario and ONE score_version
run_reps_for_combo <- function(score_ver, target_group, n_screened,
                               threshold_type, sens_target, spec_target, snr) {

  true_beta <- log(or_vec[target_group])

  reps <- future_map_dfr(1:n_reps, function(i) {
    pool <- make_pool(snr, score_ver, target_group, n_screened,
                      threshold_type, sens_target, spec_target)
    hard_res <- hard_design(pool)
    soft_res <- soft_design(pool)
    tibble(
      hard = hard_res$pval < alpha,
      soft = soft_res$pval < alpha,
      ess_hard = hard_res$ess,
      ess_soft = soft_res$ess,
      beta_hard = hard_res$beta_hat,
      beta_soft = soft_res$beta_hat,
      sens_emp = mean(pool$test_positive[pool$is_target]),
      spec_emp = mean(1 - pool$test_positive[!pool$is_target])
    )
  }, .options = furrr_options(seed = TRUE))

  summarise(reps,
            power_hard = mean(hard),
            power_soft = mean(soft),
            mean_ess_hard = mean(ess_hard, na.rm = TRUE),
            mean_ess_soft = mean(ess_soft, na.rm = TRUE),
            bias_hard = mean(beta_hard - true_beta, na.rm = TRUE),
            bias_soft = mean(beta_soft - true_beta, na.rm = TRUE),
            mse_hard  = mean((beta_hard - true_beta)^2, na.rm = TRUE),
            mse_soft  = mean((beta_soft - true_beta)^2, na.rm = TRUE),
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


# --- 4. MAIN LOOP ----------------------------------------------------------
# Use pmap to iterate over the rows of the scenario_definitions tibble
results <- pmap_dfr(scenario_definitions, .f = run_one_scenario)


# --- 5. SAVE and PRINT -----------------------------------------------------
print(results, n = 80)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
file_out <- "results/tables/aim4_logistic_multi_group.tsv"
write_tsv(results, file_out)
cat("Saved to", file_out, "\n")