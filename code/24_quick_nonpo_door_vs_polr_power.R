## =====================================================================
##  24_quick_nonpo_door_vs_polr_power.R
##  Quick power check under the non-PO DGM:
##   - treatment affects death only
##   - compare Ordinal (polr) vs DOOR
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")
source("code/functions/01_simulation_helpers.R")
source("code/functions/02_metrics.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table

or_vector <- setNames(subphenotypes$or_arrest_shrunk, subphenotypes$subphenotype)
prevalence <- setNames(subphenotypes$prevalence, subphenotypes$subphenotype)
baseline <- setNames(subphenotypes$baseline_mortality, subphenotypes$subphenotype)

n <- as.integer(Sys.getenv("SAMPLE_SIZE", unset = "2500"))
n_reps <- as.integer(Sys.getenv("N_REPS", unset = "100"))
accuracy_grid <- as.numeric(strsplit(Sys.getenv("ACCURACY_GRID", unset = "0.7,1.0"), ",")[[1]])
seed_base <- as.integer(Sys.getenv("SEED_BASE", unset = "24001"))
alpha <- params$alpha_primary
or_upper_levels <- as.numeric(Sys.getenv("OR_UPPER_LEVELS", unset = "1.0"))
run_tag <- Sys.getenv("RUN_TAG", unset = "quick")
tag_suffix <- if (nzchar(run_tag)) paste0("_", run_tag) else ""

survivor_split <- c(0.08, 0.12, 0.20, 0.25, 0.35)

output_root <- "results/supp/door_comparison"
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)

build_cum_control <- function(p_dead, survivor_split) {
  nondeath <- cumsum((1 - p_dead) * survivor_split)
  c(p_dead, p_dead + nondeath[1:4])
}

simulate_nonPO_trial <- function(or_vector, prevalence, baseline,
                                 n, seed,
                                 or_upper_levels = 1.0) {
  set.seed(seed)
  groups <- names(or_vector)

  dat <- tibble(
    id = seq_len(n),
    group = sample(groups, size = n, replace = TRUE, prob = prevalence),
    treatment = rbinom(n, 1, 0.5)
  )

  out <- map2_int(dat$group, dat$treatment, function(g, t) {
    p0 <- baseline[[g]]
    cum_ctrl <- build_cum_control(p0, survivor_split)
    theta <- qlogis(cum_ctrl)
    betas <- c(log(or_vector[[g]]), rep(log(or_upper_levels), 4))
    cum_t <- plogis(theta + betas * t)
    probs <- diff(c(0, cum_t, 1))
    probs <- pmax(probs, 0)
    probs <- probs / sum(probs)
    sample(1:6, size = 1, prob = probs)
  })

  dat %>% mutate(outcome6 = out)
}

fit_door_win_odds <- function(data) {
  trt <- data %>% filter(treatment == 1) %>% pull(outcome)
  ctrl <- data %>% filter(treatment == 0) %>% pull(outcome)
  n1 <- length(trt)
  n0 <- length(ctrl)
  if (n1 < 5 || n0 < 5) {
    return(tibble(log_or_hat = NA_real_, se = NA_real_))
  }

  wins <- sum(outer(trt, ctrl, `>`))
  losses <- sum(outer(trt, ctrl, `<`))
  if (wins == 0 || losses == 0) {
    return(tibble(log_or_hat = NA_real_, se = NA_real_))
  }

  # Same direction convention as log-OR: positive = harmful.
  log_hat <- log(losses / wins)
  se <- sqrt(1 / wins + 1 / losses)
  tibble(log_or_hat = log_hat, se = se)
}

compute_nonPO_polr_truth <- function(p_dead, survivor_split,
                                     or_death, or_upper = 1.0) {
  cum_ctrl <- build_cum_control(p_dead, survivor_split)
  betas <- c(log(or_death), rep(log(or_upper), 4))
  cum_c <- c(cum_ctrl, 1)
  wts <- cum_c[1:5] * (1 - cum_c[1:5])
  wts <- wts / sum(wts)
  sum(wts * betas)
}

compute_true_log_door_effect_nonpo <- function(p_dead, survivor_split,
                                               or_death, or_upper = 1.0) {
  cum_ctrl <- build_cum_control(p_dead, survivor_split)
  theta <- qlogis(cum_ctrl)
  betas <- c(log(or_death), rep(log(or_upper), 4))
  cum_trt <- plogis(theta + betas)
  p_ctrl <- diff(c(0, cum_ctrl, 1))
  p_trt <- diff(c(0, cum_trt, 1))
  idx <- seq_along(p_ctrl)
  n_win <- sum(outer(idx, idx, function(t, c) (t > c) * p_trt[t] * p_ctrl[c]))
  n_loss <- sum(outer(idx, idx, function(t, c) (t < c) * p_trt[t] * p_ctrl[c]))
  log(n_loss / n_win)
}

true_log_polr <- setNames(
  map_dbl(names(or_vector), function(g) {
    compute_nonPO_polr_truth(baseline[[g]], survivor_split, or_vector[[g]], or_upper_levels)
  }),
  names(or_vector)
)

true_log_door <- setNames(
  map_dbl(names(or_vector), function(g) {
    compute_true_log_door_effect_nonpo(baseline[[g]], survivor_split, or_vector[[g]], or_upper_levels)
  }),
  names(or_vector)
)

run_accuracy <- function(accuracy) {
  map_dfr(seq_len(n_reps), function(k) {
    sim <- simulate_nonPO_trial(
      or_vector = or_vector,
      prevalence = prevalence,
      baseline = baseline,
      n = n,
      seed = seed_base + k,
      or_upper_levels = or_upper_levels
    )

    sim$observed_group <- misclassify_groups(
      true_group = sim$group,
      prevalence = prevalence,
      accuracy = accuracy,
      seed = seed_base + 10000 + k
    )

    groups <- sort(unique(sim$observed_group))
    map_dfr(groups, function(g) {
      d <- sim %>% filter(observed_group == g)
      bind_rows(
        fit_polr_or(d %>% transmute(treatment, outcome = outcome6)) %>%
          mutate(model_type = "Ordinal (polr)", log_or_true = true_log_polr[[g]]),
        fit_door_win_odds(d %>% transmute(treatment, outcome = outcome6)) %>%
          mutate(model_type = "DOOR", log_or_true = true_log_door[[g]])
      ) %>%
        mutate(group = g, rep_id = k, accuracy = accuracy, sample_size = n)
    })
  })
}

raw_results <- map_dfr(accuracy_grid, run_accuracy)

metric_results <- raw_results %>%
  group_by(sample_size, accuracy, group, model_type) %>%
  group_modify(~ {
    truth <- unique(.x$log_or_true)
    is_null <- abs(truth) < 1e-6
    summarise_metrics(.x, log_or_true = truth, alpha = alpha, is_null = is_null)
  }) %>%
  ungroup()

out_path <- file.path(output_root, "tables", paste0("quick_nonpo_door_vs_polr_power", tag_suffix, ".tsv"))
write_tsv(metric_results, out_path)
message("Wrote ", out_path)
