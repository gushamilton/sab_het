# -----------------------------------------------------------------------------
# 08_compare_soft_hard.R  (updated – logistic-score version)
# -----------------------------------------------------------------------------
# This version replaces the Beta/κ mechanism with a *logistic latent predictor*:
#   z | target     ~ Normal(mu1, 1)
#   z | non-target ~ Normal(mu0, 1)             (signal-to-noise ratio SNR = mu1−mu0)
#   p = plogis(z)                               (0–1 probability score)
# A single threshold τ on p creates the diagnostic test’s binary flag.
#   • threshold_type = "sens" fixes sensitivity to sens_target.
#   • threshold_type = "spec" fixes specificity to spec_target.
#   • threshold_type = "youden" maximises Youden index (sens+spec-1).
# The HARD design sees only the binary flag; the SOFT design sees a coarsened
# or continuous version of p.
# The script outputs empirical power_hard / power_soft *and* the realised mean
# sensitivity & specificity for each SNR setting.
# -----------------------------------------------------------------------------

# --- 1. SET-UP -------------------------------------------------------------
pacman::p_load(tidyverse, broom, furrr)
set.seed(20295)
plan(multisession, workers = max(1, parallel::detectCores() - 1))

# --- 2. SCENARIO & USER PARAMETERS ----------------------------------------

or_vec   <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
freq_vec <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)

p0_raw <- c(A = 7/60,  B = 0/52,  C = 11/138, D = 11/69, E = 1/69)
p0_raw["B"] <- 0.5 / (52 + 0.5)  # continuity correction

target_group  <- "C"
n_screened    <- 3284

# Latent predictor settings --------------------------------------------------
# Define multiple predictor strengths (signal-to-noise ratios)
snr_levels <- c(1, 3)  # e.g. 1 = modest, 3 = strong predictor

# Different discretisations of the score available to the SOFT design
score_versions <- c("binary", "tertile", "quartile", "decile", "continuous")

# Threshold policy -----------------------------------------------------------
threshold_type <- "sens"   # "sens", "spec", "youden"
sens_target    <- 0.90
spec_target    <- 0.95

# Simulation ---------------------------------------------------------------
n_reps <- 500
alpha  <- 0.05

# --- 3. Helpers ------------------------------------------------------------
convert_or <- function(or, p0) (or * p0) / (1 - p0 + or * p0)

# Determine threshold tau for a given score vector --------------------------
compute_tau <- function(p, is_target) {
  if (threshold_type == "youden") {
    # grid search over unique p values
    cuts <- sort(unique(p))
    youden <- sapply(cuts, function(t) {
      se <- mean(p[is_target] > t)
      sp <- mean(p[!is_target] <= t)
      se + sp - 1
    })
    return(cuts[which.max(youden)])
  } else if (threshold_type == "sens") {
    rootfun <- function(t) mean(p[is_target] > t) - sens_target
    return(uniroot(rootfun, c(0,1))$root)
  } else if (threshold_type == "spec") {
    rootfun <- function(t) mean(p[!is_target] <= t) - spec_target
    return(uniroot(rootfun, c(0,1))$root)
  } else stop("Unknown threshold_type")
}

# Generate one pool ---------------------------------------------------------
make_pool <- function(snr, score_version) {
  # Means chosen so that mu1 - mu0 = snr, centred at 0
  mu1 <-  snr/2; mu0 <- -snr/2
  grp <- sample(names(freq_vec), n_screened, TRUE, freq_vec)
  is_target <- grp == target_group
  z <- rnorm(n_screened, mean = ifelse(is_target, mu1, mu0), sd = 1)
  p <- plogis(z)

  tau <- compute_tau(p, is_target)
  test_positive <- as.integer(p > tau)

# Expose chosen discretisation ----------------------------------------------
  p_exposed <- switch(score_version,
    continuous = p,
    decile     = ntile(p, 10)/10,
    quartile   = ntile(p, 4)/4,
    tertile    = ntile(p, 3)/3,
    binary     = test_positive,
    stop("bad score_version"))

# Return pool ---------------------------------------------------------------
  tibble(group = grp, is_target, test_positive, p_target = p_exposed, p, tau)
}

# HARD analysis -------------------------------------------------------------
hard_design <- function(pool) {
  enrol <- pool %>% filter(test_positive == 1)
  n_enrolled <- nrow(enrol)
  treat <- rbinom(n_enrolled, 1, 0.5)
  p0 <- p0_raw[enrol$group]
  p1 <- convert_or(or_vec[enrol$group], p0)
  y  <- rbinom(n_enrolled, 1, ifelse(treat == 1, p1, p0))

  # Fit GLM with error handling
  fit <- tryCatch(glm(y ~ treat, family = binomial()), error = function(e) NULL)

  if (is.null(fit) || nrow(coef(summary(fit))) < 2 || anyNA(coef(summary(fit))[2,])) {
    return(list(pval = 1, n_enrolled = n_enrolled, n_treated = sum(treat==1),
                n_control = sum(treat==0), ess = n_enrolled, beta_hat = NA))
  }

  coef_row <- coef(summary(fit))[2, ]
  beta_hat <- as.numeric(coef_row["Estimate"])
  pval <- as.numeric(coef_row["Pr(>|z|)"])
  ess <- n_enrolled  # weight =1 each
  list(pval = pval, n_enrolled = n_enrolled, n_treated = sum(treat==1),
       n_control = sum(treat==0), ess = ess, beta_hat = beta_hat)
}

# SOFT analysis -------------------------------------------------------------
soft_design <- function(pool) {
  # --------------------------------------------------------------
  # Soft-enrichment weighting
  #   Use *normalised importance weights* so that, in expectation,
  #   non-targets contribute zero information and the estimator is
  #   unbiased.  See discussion in comments.
  #   w_i = p_i / mean(p)
  # --------------------------------------------------------------
  treat_prob <- 0.5 * pool$p_target
  treat <- rbinom(nrow(pool), 1, treat_prob)
  n_enrolled <- sum(pool$p_target > 0)
  n_treated  <- sum(treat == 1)
  n_control  <- sum(treat == 0 & pool$p_target > 0)
  # Use natural 0-1 weights. This is the fair, non-cheating approach.
  # The binary-score soft design will now be identical to the hard design.
  w <- pool$p_target
  ess <- sum(w)
  p0 <- p0_raw[pool$group]
  p1 <- convert_or(or_vec[pool$group], p0)
  y  <- rbinom(nrow(pool), 1, ifelse(treat == 1, p1, p0))

  # Fit GLM with error handling
  fit <- tryCatch(glm(y ~ treat, family = binomial(), weights = w), error = function(e) NULL)

  if (is.null(fit) || nrow(coef(summary(fit))) < 2 || anyNA(coef(summary(fit))[2,])) {
    return(list(pval = 1, n_enrolled = n_enrolled, n_treated = n_treated,
                n_control = n_control, ess = ess, beta_hat = NA))
  }

  coef_row <- coef(summary(fit))[2, ]
  beta_hat <- as.numeric(coef_row["Estimate"])
  pval <- as.numeric(coef_row["Pr(>|z|)"])
  list(pval = pval, n_enrolled = n_enrolled, n_treated = n_treated,
       n_control = n_control, ess = ess, beta_hat = beta_hat)
}

# simulate for one combination of snr and discretisation
sim_one_combo <- function(snr, ver) {
  reps <- future_map_dfr(1:n_reps, function(i) {
    pool <- make_pool(snr, ver)
    hard_res <- hard_design(pool)
    soft_res <- soft_design(pool)
    tibble(hard = hard_res$pval < alpha,
           soft = soft_res$pval < alpha,
           n_enrol_hard = hard_res$n_enrolled,
           n_enrol_soft = soft_res$n_enrolled,
           n_treat_hard = hard_res$n_treated,
           n_treat_soft = soft_res$n_treated,
           ess_hard = hard_res$ess,
           ess_soft = soft_res$ess,
           beta_hard = hard_res$beta_hat,
           beta_soft = soft_res$beta_hat,
           sens_emp = mean(pool$test_positive[pool$is_target]),
           spec_emp = mean(1 - pool$test_positive[!pool$is_target]))
  }, .options = furrr_options(seed = TRUE))
  summarise(reps,
            power_hard = mean(hard),
            power_soft = mean(soft),
            mean_enrol_hard = mean(n_enrol_hard),
            se_enrol_hard   = sd(n_enrol_hard)/sqrt(n_reps),
            mean_enrol_soft = mean(n_enrol_soft),
            se_enrol_soft   = sd(n_enrol_soft)/sqrt(n_reps),
            mean_treat_hard = mean(n_treat_hard, na.rm=TRUE),
            se_treat_hard   = sd(n_treat_hard, na.rm=TRUE)/sqrt(n_reps),
            mean_treat_soft = mean(n_treat_soft),
            se_treat_soft   = sd(n_treat_soft)/sqrt(n_reps),
            mean_ess_hard  = mean(ess_hard),
            se_ess_hard    = sd(ess_hard) / sqrt(n_reps),
            mean_ess_soft  = mean(ess_soft),
            se_ess_soft    = sd(ess_soft) / sqrt(n_reps),
            mean_beta_hard = mean(beta_hard, na.rm = TRUE),
            se_beta_hard   = sd(beta_hard, na.rm = TRUE) / sqrt(n_reps),
            mean_beta_soft = mean(beta_soft, na.rm = TRUE),
            se_beta_soft   = sd(beta_soft, na.rm = TRUE) / sqrt(n_reps),
            bias_hard = mean(beta_hard - log(or_vec[target_group]), na.rm = TRUE),
            mse_hard  = mean((beta_hard - log(or_vec[target_group]))^2, na.rm = TRUE),
            bias_soft = mean(beta_soft - log(or_vec[target_group]), na.rm = TRUE),
            mse_soft  = mean((beta_soft - log(or_vec[target_group]))^2, na.rm = TRUE),
            true_beta = log(or_vec[target_group]),
            sens_emp   = mean(sens_emp),
            spec_emp   = mean(spec_emp),
            score_version = ver,
            snr = snr)
}

# --- 4. MAIN LOOP ----------------------------------------------------------
results <- purrr::map_dfr(snr_levels, function(snr) {
  purrr::map_dfr(score_versions, function(ver) {
    sim_one_combo(snr, ver)
  })
})
print(results)
view(results)
# --- 5. SAVE ---------------------------------------------------------------
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
file_out <- sprintf("results/tables/aim4_logistic_%s_multi.tsv", target_group)
write_tsv(results, file_out)
cat("Saved to", file_out, "\n") 
