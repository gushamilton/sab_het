## =====================================================================
##  14_run_ordinal6_same_effect_comparison.R
##  Simulate 6-point ordinal outcomes, collapse to binary from same data,
##  and compare power/bias for ordinal vs binary models.
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

sample_sizes <- params$sample_sizes
sample_override <- Sys.getenv("SAMPLE_SIZES", unset = "")
if (nzchar(sample_override)) {
  sample_sizes <- as.integer(strsplit(sample_override, ",")[[1]])
}

accuracy_grid <- c(0.70, 1.00)
n_reps <- as.integer(Sys.getenv("N_REPS_ORD6", unset = "1000"))
seed_base <- as.integer(Sys.getenv("SEED_BASE", unset = "14001"))
rep_offset <- as.integer(Sys.getenv("REP_OFFSET", unset = "0"))
alpha <- params$alpha_primary
calibration_n <- as.integer(Sys.getenv("CALIBRATION_N", unset = "200000"))
run_tag <- Sys.getenv("RUN_TAG", unset = "")
skip_calibration <- tolower(Sys.getenv("SKIP_CALIBRATION", unset = "0")) %in% c("1", "true", "yes")
tag_suffix <- if (nzchar(run_tag)) paste0("_", run_tag) else ""

# Five non-death levels for the 6-point total outcome.
# Dead is level 1; levels 2-6 split the survivor mass.
survivor_split <- c(0.08, 0.12, 0.20, 0.25, 0.35)
if (abs(sum(survivor_split) - 1) > 1e-8) {
  stop("survivor_split must sum to 1.")
}

output_root <- "results/supp/ordinal6_sameeffect"
dir.create(file.path(output_root, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)

build_cum_control <- function(p_dead, survivor_split) {
  # Cumulative probabilities for levels <=1..<=5 under control.
  nondeath <- cumsum((1 - p_dead) * survivor_split)
  c(p_dead, p_dead + nondeath[1:4])
}

simulate_ordinal6_trial <- function(or_vector, prevalence, baseline, n, seed) {
  set.seed(seed)
  groups <- names(or_vector)

  dat <- tibble(
    id = seq_len(n),
    group = sample(groups, size = n, replace = TRUE, prob = prevalence),
    treatment = rbinom(n, 1, 0.5)
  )

  # Simulate from proportional-odds latent cumulative logits.
  # logit(P(Y<=k|t)) = theta_k + beta*t
  # theta_1 is calibrated so the death OR matches or_vector exactly.
  out <- map2_int(dat$group, dat$treatment, function(g, t) {
    p0 <- baseline[[g]]
    beta <- log(or_vector[[g]])
    cum_ctrl <- build_cum_control(p0, survivor_split)
    theta <- qlogis(cum_ctrl)
    cum_t <- plogis(theta + beta * t)
    probs <- diff(c(0, cum_t, 1))
    probs <- pmax(probs, 0)
    probs <- probs / sum(probs)
    sample(1:6, size = 1, prob = probs)
  })

  dat %>% mutate(outcome6 = out)
}

fit_both_models_by_group <- function(sim_data) {
  groups <- sort(unique(sim_data$observed_group))
  map_dfr(groups, function(g) {
    d <- sim_data %>% filter(observed_group == g)

    ord_fit <- fit_polr_or(d %>% transmute(treatment, outcome = outcome6))
    bin_fit <- fit_logistic_or(d %>% transmute(treatment, outcome = as.integer(outcome6 == 1)))

    bind_rows(
      ord_fit %>% mutate(model_type = "Ordinal", group = g),
      bin_fit %>% mutate(model_type = "Binary", group = g)
    )
  })
}

run_scenario <- function(n, accuracy, n_reps, seed_base, rep_offset = 0L) {
  map_dfr(seq_len(n_reps), function(k) {
    rep_id <- rep_offset + k
    sim <- simulate_ordinal6_trial(
      or_vector = or_vector,
      prevalence = prevalence,
      baseline = baseline,
      n = n,
      seed = seed_base + rep_id
    )
    sim$observed_group <- misclassify_groups(
      true_group = sim$group,
      prevalence = prevalence,
      accuracy = accuracy,
      seed = seed_base + 10000 + rep_id
    )

    fit_both_models_by_group(sim) %>%
      mutate(rep_id = rep_id, sample_size = n, accuracy = accuracy)
  })
}

calibrate_once <- function(n, seed = 999) {
  sim <- simulate_ordinal6_trial(or_vector, prevalence, baseline, n = n, seed = seed)
  sim$observed_group <- sim$group
  fits <- fit_both_models_by_group(sim)
  fits %>%
    mutate(
      log_or_true = log(or_vector[group]),
      or_hat = exp(log_or_hat),
      or_true = exp(log_or_true),
      abs_err_log_or = abs(log_or_hat - log_or_true)
    ) %>%
    select(group, model_type, or_true, or_hat, log_or_true, log_or_hat, abs_err_log_or)
}

if (!skip_calibration) {
  message("Running calibration check...")
  calibration <- calibrate_once(calibration_n)
  write_tsv(calibration, file.path(output_root, "tables", paste0("ordinal6_sameeffect_calibration", tag_suffix, ".tsv")))
}

scenario_grid <- expand_grid(sample_size = sample_sizes, accuracy = accuracy_grid)
message("Running simulation grid for 6-point ordinal vs collapsed binary...")

raw_results <- pmap_dfr(scenario_grid, function(sample_size, accuracy) {
  message(sprintf("  N=%d, accuracy=%.2f", sample_size, accuracy))
  run_scenario(sample_size, accuracy, n_reps, seed_base, rep_offset = rep_offset)
}) %>%
  mutate(log_or_true = log(or_vector[group]))

write_tsv(raw_results, file.path(output_root, "data", paste0("ordinal6_sameeffect_raw", tag_suffix, ".tsv.gz")))

metric_results <- raw_results %>%
  group_by(sample_size, accuracy, group, model_type) %>%
  group_modify(~ {
    truth <- log(or_vector[[.y$group]])
    is_null <- abs(truth) < 1e-6
    summarise_metrics(.x, log_or_true = truth, alpha = alpha, is_null = is_null)
  }) %>%
  ungroup()

write_tsv(metric_results, file.path(output_root, "tables", paste0("ordinal6_sameeffect_metrics", tag_suffix, ".tsv")))

bias_results <- raw_results %>%
  filter(group != "A") %>%
  group_by(sample_size, accuracy, group, model_type) %>%
  summarise(
    mean_bias_log_or = mean(log_or_hat - log_or_true, na.rm = TRUE),
    mean_abs_bias_log_or = mean(abs(log_or_hat - log_or_true), na.rm = TRUE),
    rmse_log_or = sqrt(mean((log_or_hat - log_or_true)^2, na.rm = TRUE)),
    .groups = "drop"
  )

write_tsv(bias_results, file.path(output_root, "tables", paste0("ordinal6_sameeffect_bias", tag_suffix, ".tsv")))

power_plot <- metric_results %>%
  filter(group != "A") %>%
  ggplot(aes(x = sample_size, y = power, color = model_type)) +
  geom_line() +
  geom_point(size = 1.4) +
  facet_grid(group ~ accuracy, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "6-point ordinal vs collapsed binary from same simulated data",
    subtitle = "Power (correct-direction significant effects only)",
    x = "Sample size",
    y = "Power",
    color = "Model"
  ) +
  theme_minimal()

ggsave(
  file.path(output_root, "plots", paste0("ordinal6_sameeffect_power", tag_suffix, ".pdf")),
  power_plot,
  width = 12,
  height = 9
)

bias_plot <- bias_results %>%
  ggplot(aes(x = sample_size, y = mean_abs_bias_log_or, color = model_type)) +
  geom_line() +
  geom_point(size = 1.4) +
  facet_grid(group ~ accuracy, scales = "free_y") +
  labs(
    title = "6-point ordinal vs collapsed binary from same simulated data",
    subtitle = "Mean absolute bias in log(OR)",
    x = "Sample size",
    y = "Mean absolute bias",
    color = "Model"
  ) +
  theme_minimal()

ggsave(
  file.path(output_root, "plots", paste0("ordinal6_sameeffect_bias", tag_suffix, ".pdf")),
  bias_plot,
  width = 12,
  height = 9
)

message("Done.")
