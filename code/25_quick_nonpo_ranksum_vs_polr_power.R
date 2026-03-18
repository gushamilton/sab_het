## =====================================================================
##  25_quick_nonpo_ranksum_vs_polr_power.R
##  Quick power-only comparison under non-PO:
##   - Ordinal (polr)
##   - Rank-sum (Wilcoxon / Mann-Whitney)
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
seed_base <- as.integer(Sys.getenv("SEED_BASE", unset = "25001"))
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

compute_nonPO_polr_truth <- function(p_dead, survivor_split,
                                     or_death, or_upper = 1.0) {
  cum_ctrl <- build_cum_control(p_dead, survivor_split)
  betas <- c(log(or_death), rep(log(or_upper), 4))
  cum_c <- c(cum_ctrl, 1)
  wts <- cum_c[1:5] * (1 - cum_c[1:5])
  wts <- wts / sum(wts)
  sum(wts * betas)
}

true_log_polr <- setNames(
  map_dbl(names(or_vector), function(g) {
    compute_nonPO_polr_truth(baseline[[g]], survivor_split, or_vector[[g]], or_upper_levels)
  }),
  names(or_vector)
)

fit_ranksum_test <- function(data) {
  if (nrow(data) < 10 || length(unique(data$treatment)) < 2 || length(unique(data$outcome)) < 2) {
    return(tibble(p_value = NA_real_, sign_hat = NA_real_))
  }

  test <- tryCatch(
    suppressWarnings(wilcox.test(outcome ~ treatment, data = data, exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(test)) {
    return(tibble(p_value = NA_real_, sign_hat = NA_real_))
  }

  trt_mean <- mean(data$outcome[data$treatment == 1], na.rm = TRUE)
  ctrl_mean <- mean(data$outcome[data$treatment == 0], na.rm = TRUE)

  tibble(
    p_value = test$p.value,
    sign_hat = sign(ctrl_mean - trt_mean) # positive = harmful, aligned with log-OR sign
  )
}

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

      polr_fit <- fit_polr_or(d %>% transmute(treatment, outcome = outcome6)) %>%
        mutate(model_type = "Ordinal (polr)", log_or_true = true_log_polr[[g]])

      rank_fit <- fit_ranksum_test(d %>% transmute(treatment, outcome = outcome6)) %>%
        mutate(model_type = "Rank-sum", log_or_true = log(or_vector[[g]]))

      bind_rows(
        polr_fit %>% transmute(
          model_type,
          p_value = {
            z <- qnorm(1 - alpha / 2)
            ci_lower <- log_or_hat - z * se
            ci_upper <- log_or_hat + z * se
            significant <- !is.na(log_or_hat) & !is.na(se) & (ci_lower > 0 | ci_upper < 0)
            if_else(significant, alpha / 2, 1)
          },
          sign_hat = sign(log_or_hat),
          log_or_true
        ),
        rank_fit
      ) %>%
        mutate(group = g, rep_id = k, accuracy = accuracy, sample_size = n)
    })
  })
}

raw_results <- map_dfr(accuracy_grid, run_accuracy)

metric_results <- raw_results %>%
  group_by(sample_size, accuracy, group, model_type) %>%
  summarise(
    type1 = if (first(group) == "A") mean(p_value < alpha, na.rm = TRUE) else NA_real_,
    power = if (first(group) != "A") mean((p_value < alpha) & (sign_hat == sign(first(log_or_true))), na.rm = TRUE) else NA_real_,
    type_s = if (first(group) != "A") mean((p_value < alpha) & (sign_hat != sign(first(log_or_true))), na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

out_path <- file.path(output_root, "tables", paste0("quick_nonpo_ranksum_vs_polr_power", tag_suffix, ".tsv"))
write_tsv(metric_results, out_path)
message("Wrote ", out_path)
