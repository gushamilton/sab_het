# -----------------------------------------------------------------------------
# 17_final_sandwich_latent_simulation.R
# -----------------------------------------------------------------------------
# Upgrades:
#   • robust (sandwich) SE in analyse_trial()
#   • latent-class mixture estimator (Design E)
# -----------------------------------------------------------------------------

# --- 1. SET-UP -------------------------------------------------------------
pacman::p_load(
  tidyverse, broom, furrr,
  sandwich, lmtest,       # NEW: robust variance
  flexmix                 # NEW: latent-class mixture
)

set.seed(2025, kind = "L'Ecuyer-CMRG")

slurm_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
n_cores <- if (nzchar(slurm_cores)) as.integer(slurm_cores) else max(1, parallel::detectCores() - 1)
cat("Setting up", n_cores, "parallel workers.\n")
plan(multisession, workers = n_cores)

# --- 2. SCENARIO & USER PARAMETERS  (unchanged) ----------------------------
or_vec   <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
freq_vec <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)
p0_raw   <- c(A = 7/60,  B = 0/52,  C = 11/138, D = 11/69, E = 1/69)
p0_raw["B"] <- 0.5 / (52 + 0.5)

score_versions  <- c("binary", "quartile", "continuous")

base_scenarios <- tidyr::crossing(
  threshold_type = c("sens", "spec"),
  snr = c(1, 3)
) %>%
  mutate(
    sens_target = if_else(threshold_type == "sens", 0.95, NA_real_),
    spec_target = if_else(threshold_type == "spec", 0.95, NA_real_)
  )

scenario_definitions <- tidyr::crossing(
  group_info = tribble(
    ~target_group, ~n_to_randomize,
    "B",           500,
    "C",           6900
  ),
  base_scenarios,
  score_version = score_versions
) %>%
  unnest(group_info) %>%
  select(target_group, n_to_randomize, score_version, everything())

n_reps <- 300
alpha  <- 0.05
PILOT_POOL_SIZE <- 1e5
BATCH_SIZE      <- 200

# --- 3. HELPER & ANALYSIS FUNCTIONS ---------------------------------------
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
    quantile(pool$p[pool$is_target], probs = 1 - sens_target, type = 1)
  } else if (threshold_type == "spec") {
    quantile(pool$p[!pool$is_target], probs = spec_target, type = 1)
  }
}

# ---------- 3a.  Core analysis with robust SE  -----------------------------
analyse_trial <- function(cohort, weights_vec, true_beta) {
  n <- nrow(cohort)
  if (n < 10) return(list(power = NA, bias = NA, mse = NA))

  cohort$treat <- rbinom(n, 1, 0.5)

  p0 <- p0_raw[cohort$group]
  p1 <- convert_or(or_vec[cohort$group], p0)
  cohort$y <- rbinom(n, 1, ifelse(cohort$treat == 1, p1, p0))

  fit <- tryCatch(
    glm(y ~ treat, data = cohort, family = binomial(),
        weights = weights_vec),
    error = function(e) NULL
  )
  if (is.null(fit) || length(coef(fit)) < 2) {
    return(list(power = 0, bias = NA, mse = NA))
  }

  # --- robust z, p-value ---------------------------------------------------
  vcov_hc  <- sandwich::vcovHC(fit, type = "HC0")
  se_rob   <- sqrt(vcov_hc["treat","treat"])
  beta_hat <- coef(fit)["treat"]
  z_rob    <- beta_hat / se_rob
  pval_rob <- 2 * (1 - pnorm(abs(z_rob)))

  list(
    power = as.integer(pval_rob < alpha),
    bias  = beta_hat - true_beta,
    mse   = (beta_hat - true_beta)^2
  )
}

# ---------- 3b.  Latent-class mixture estimator (Design E) -----------------
analyse_trial_mixture <- function(cohort, true_beta) {
  n <- nrow(cohort)
  if (n < 10) return(list(power = NA, bias = NA, mse = NA))

  # outcome + treatment must already exist
  dat <- cohort %>% select(y, treat, p)

  res <- tryCatch({
    mix <- flexmix(y ~ treat | p, data = dat, k = 2,
                   model = FLXMRglm(family = "binomial"))
    # identify component with higher mean biomarker score
    comp_means <- tapply(dat$p, clusters(mix), mean)
    target_comp <- names(which.max(comp_means))
    beta_hat <- parameters(mix)["treat", target_comp]

    vc <- vcov(mix)[[as.integer(target_comp)]]
    se   <- sqrt(diag(vc))["treat"]
    z    <- beta_hat / se
    pval <- 2 * (1 - pnorm(abs(z)))

    list(
      power = as.integer(pval < alpha),
      bias  = beta_hat - true_beta,
      mse   = (beta_hat - true_beta)^2
    )
  }, error = function(e) list(power = NA, bias = NA, mse = NA))

  res
}

# --- 4. MAIN SIMULATION FUNCTION ------------------------------------------
run_full_scenario <- function(target_group, n_to_randomize, score_version,
                              threshold_type, snr, sens_target, spec_target) {

  cat("Running Scenario: Group", target_group, "| N_to_randomize", n_to_randomize,
      "| Score", score_version, "| SNR", snr, "| Type", threshold_type, "\n")

  true_beta <- log(or_vec[target_group])

  run_one_rep <- function(rep_id) {

    pilot_pool <- draw_patient_batch(PILOT_POOL_SIZE, target_group, snr)
    tau <- compute_tau(pilot_pool, threshold_type, sens_target, spec_target)

    # empirical sens/spec for reporting only
    sens_emp <- mean(pilot_pool$p[pilot_pool$is_target] > tau)
    spec_emp <- mean(pilot_pool$p[!pilot_pool$is_target] <= tau)

    # -------- Recruitment: Hard ------------------------------------------------
    n_screened_hard <- 0; enrolled_hard <- list()
    while (nrow(bind_rows(enrolled_hard)) < n_to_randomize) {
      new_batch <- draw_patient_batch(BATCH_SIZE, target_group, snr)
      n_screened_hard <- n_screened_hard + BATCH_SIZE
      eligible <- new_batch %>% filter(p > tau)
      if (nrow(eligible) > 0) enrolled_hard[[length(enrolled_hard)+1]] <- eligible
      if (n_screened_hard > 10 * PILOT_POOL_SIZE &&
          n_screened_hard > 2 * n_to_randomize / 0.01) break
    }
    cohort_hard <- bind_rows(enrolled_hard) %>% head(n_to_randomize)

    # -------- Recruitment: Probabilistic --------------------------------------
    n_screened_prob <- 0; enrolled_prob <- list()
    while (nrow(bind_rows(enrolled_prob)) < n_to_randomize) {
      new_batch <- draw_patient_batch(BATCH_SIZE, target_group, snr)
      n_screened_prob <- n_screened_prob + BATCH_SIZE

      p_enroll <- switch(score_version,
        continuous = new_batch$p,
        quartile   = ntile(new_batch$p, 4) / 4,
        binary     = as.numeric(new_batch$p > tau)
      )
      p_enroll[p_enroll < 0.1] <- 0  # keeps positivity but trims tails
      new_batch$p_enroll <- p_enroll

      picked <- rbinom(nrow(new_batch), 1, p_enroll) == 1
      if (sum(picked) > 0) enrolled_prob[[length(enrolled_prob)+1]] <- new_batch[picked,]
      if (n_screened_prob > 10 * PILOT_POOL_SIZE &&
          n_screened_prob > 2 * n_to_randomize / 0.01) break
    }
    cohort_prob <- bind_rows(enrolled_prob) %>% head(n_to_randomize)

    # -------- Analysis ---------------------------------------------------------
    weights_B <- switch(score_version,
      continuous = cohort_hard$p,
      quartile   = ntile(cohort_hard$p, 4) / 4,
      binary     = NULL)

    weights_D <- 1 / pmax(cohort_prob$p_enroll, 1e-8)  # stabilise

    res_A <- analyse_trial(cohort_hard, weights_vec = NULL,       true_beta)
    res_B <- analyse_trial(cohort_hard, weights_vec = weights_B,  true_beta)
    res_C <- analyse_trial(cohort_prob, weights_vec = NULL,       true_beta)
    res_D <- analyse_trial(cohort_prob, weights_vec = weights_D,  true_beta)

    # Mixture needs outcome + treat first; simulate once and reuse
    cohort_prob$treat <- rbinom(nrow(cohort_prob), 1, 0.5)
    p0_prob <- p0_raw[cohort_prob$group]
    p1_prob <- convert_or(or_vec[cohort_prob$group], p0_prob)
    cohort_prob$y <- rbinom(nrow(cohort_prob), 1,
                            ifelse(cohort_prob$treat == 1, p1_prob, p0_prob))

    res_E <- analyse_trial_mixture(cohort_prob, true_beta)

    tibble(
      sens_emp, spec_emp, nns_hard = n_screened_hard, nns_prob = n_screened_prob,
      power_hard = res_A$power,     power_soft_weighted = res_B$power,
      power_soft_prob_unw = res_C$power, power_soft_prob_ipw = res_D$power,
      power_latent = res_E$power,
      bias_hard = res_A$bias,       bias_soft_weighted = res_B$bias,
      bias_soft_prob_unw = res_C$bias, bias_soft_prob_ipw = res_D$bias,
      bias_latent = res_E$bias,
      mse_hard = res_A$mse,         mse_soft_weighted = res_B$mse,
      mse_soft_prob_unw = res_C$mse, mse_soft_prob_ipw = res_D$mse,
      mse_latent = res_E$mse
    )
  }

  future_map_dfr(1:n_reps, run_one_rep, .options = furrr_options(seed = TRUE)) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
}

# --- 5. MAIN EXECUTION LOOP -----------------------------------------------
full_results <- pmap_dfr(scenario_definitions, run_full_scenario)

final_output <- bind_cols(scenario_definitions, full_results)

# --- 6. SAVE AND PRINT -----------------------------------------------------
print(final_output, n = nrow(final_output), width = Inf)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
file_out <- "results/tables/aim5_sandwich_latent.tsv"
write_tsv(final_output, file_out)
cat("\nSaved to", file_out, "\n")
