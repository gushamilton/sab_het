## =====================================================================
##  19_run_door_comparison.R
##  Compare ordinal (polr), binary (death cut-point), and DOOR (win odds)
##  from the same 6-point ordinal DGM. Each method is assessed against
##  its own true estimand; proportional bias enables cross-method
##  comparison on a common scale.
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

# Runtime options
sample_sizes <- params$sample_sizes
sample_override <- Sys.getenv("SAMPLE_SIZES", unset = "")
if (nzchar(sample_override)) {
  sample_sizes <- as.integer(strsplit(sample_override, ",")[[1]])
}

accuracy_grid <- c(0.70, 1.00)
n_reps <- as.integer(Sys.getenv("N_REPS_DOOR", unset = "1000"))
seed_base <- as.integer(Sys.getenv("SEED_BASE", unset = "15001"))
rep_offset <- as.integer(Sys.getenv("REP_OFFSET", unset = "0"))
alpha <- params$alpha_primary
calibration_n <- as.integer(Sys.getenv("CALIBRATION_N", unset = "200000"))
skip_cal <- tolower(Sys.getenv("SKIP_CALIBRATION", unset = "0")) %in% c("1", "true", "yes")
run_tag <- Sys.getenv("RUN_TAG", unset = "")
tag_suffix <- if (nzchar(run_tag)) paste0("_", run_tag) else ""

# Survivor mass split across levels 2-6 (level 1 = dead)
survivor_split <- c(0.08, 0.12, 0.20, 0.25, 0.35)
stopifnot(abs(sum(survivor_split) - 1) < 1e-8)

# Output directories
output_root <- "results/supp/door_comparison"
for (d in c("data", "tables", "plots")) {
  dir.create(file.path(output_root, d), showWarnings = FALSE, recursive = TRUE)
}

# =====================================================================
# 1. DGM helpers
# =====================================================================

build_cum_control <- function(p_dead, survivor_split) {
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

# =====================================================================
# 2. True estimands
# =====================================================================

compute_true_log_door_effect <- function(p_dead, survivor_split, or_value) {
  cum_ctrl <- build_cum_control(p_dead, survivor_split)
  theta <- qlogis(cum_ctrl)
  beta <- log(or_value)
  cum_trt <- plogis(theta + beta)

  p_ctrl <- diff(c(0, cum_ctrl, 1))
  p_trt <- diff(c(0, cum_trt, 1))

  k <- length(p_ctrl)
  idx <- seq_len(k)
  n_win <- sum(outer(idx, idx, function(t, c) (t > c) * p_trt[t] * p_ctrl[c]))
  n_loss <- sum(outer(idx, idx, function(t, c) (t < c) * p_trt[t] * p_ctrl[c]))

  # Direction aligned to log-OR convention used elsewhere:
  # positive = harmful, negative = protective.
  # Since higher ordinal levels are better, harmful treatment yields more losses.
  log(n_loss / n_win)
}

true_log_win_odds <- setNames(
  map_dbl(names(or_vector), function(g) {
    compute_true_log_door_effect(baseline[[g]], survivor_split, or_vector[[g]])
  }),
  names(or_vector)
)

message("True estimands:")
walk(names(or_vector), function(g) {
  message(sprintf(
    "  %s  log-OR (polr/binary) = %.3f   log-win-odds (DOOR) = %.3f",
    g, log(or_vector[[g]]), true_log_win_odds[[g]]
  ))
})

# =====================================================================
# 3. Model fitters
# =====================================================================

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

  # Align DOOR estimate direction to log-OR convention:
  # positive = harmful, negative = protective.
  log_wo <- log(losses / wins)
  se <- sqrt(1 / wins + 1 / losses)

  tibble(log_or_hat = log_wo, se = se)
}

# =====================================================================
# 4. Fit all three models from the same data
# =====================================================================

fit_all_models_by_group <- function(sim_data) {
  groups <- sort(unique(sim_data$observed_group))

  map_dfr(groups, function(g) {
    d <- sim_data %>% filter(observed_group == g)

    polr_fit <- fit_polr_or(
      d %>% transmute(treatment, outcome = outcome6)
    ) %>% mutate(model_type = "Ordinal (polr)")

    bin_fit <- fit_logistic_or(
      d %>% transmute(treatment, outcome = as.integer(outcome6 == 1))
    ) %>% mutate(model_type = "Binary (death)")

    door_fit <- fit_door_win_odds(
      d %>% transmute(treatment, outcome = outcome6)
    ) %>% mutate(model_type = "DOOR")

    bind_rows(polr_fit, bin_fit, door_fit) %>% mutate(group = g)
  })
}

# =====================================================================
# 5. Simulation runner
# =====================================================================

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
      seed = seed_base + 10000L + rep_id
    )

    fit_all_models_by_group(sim) %>%
      mutate(rep_id = rep_id, sample_size = n, accuracy = accuracy)
  })
}

# =====================================================================
# 6. Calibration check
# =====================================================================

calibrate_once <- function(n, seed = 999) {
  sim <- simulate_ordinal6_trial(or_vector, prevalence, baseline, n = n, seed = seed)
  sim$observed_group <- sim$group

  fits <- fit_all_models_by_group(sim)

  fits %>%
    mutate(
      log_or_true = case_when(
        model_type == "DOOR" ~ true_log_win_odds[group],
        TRUE ~ log(or_vector[group])
      ),
      bias = log_or_hat - log_or_true,
      prop_bias = bias / abs(log_or_true)
    ) %>%
    select(group, model_type, log_or_true, log_or_hat, bias, prop_bias)
}

if (!skip_cal) {
  message("Running calibration check (n = ", calibration_n, ")...")
  calibration <- calibrate_once(calibration_n)
  write_tsv(
    calibration,
    file.path(output_root, "tables", paste0("door_comparison_calibration", tag_suffix, ".tsv"))
  )
}

# =====================================================================
# 7. Main simulation grid
# =====================================================================

scenario_grid <- expand_grid(sample_size = sample_sizes, accuracy = accuracy_grid)
message("Running simulation grid...")

raw_results <- pmap_dfr(scenario_grid, function(sample_size, accuracy) {
  message(sprintf("  N = %d, accuracy = %.2f", sample_size, accuracy))
  run_scenario(sample_size, accuracy, n_reps, seed_base, rep_offset = rep_offset)
}) %>%
  mutate(
    log_or_true = case_when(
      model_type == "DOOR" ~ true_log_win_odds[group],
      TRUE ~ log(or_vector[group])
    )
  )

write_tsv(
  raw_results,
  file.path(output_root, "data", paste0("door_comparison_raw", tag_suffix, ".tsv.gz"))
)

# =====================================================================
# 8. Metrics
# =====================================================================

metric_results <- raw_results %>%
  group_by(sample_size, accuracy, group, model_type) %>%
  group_modify(~ {
    truth <- unique(.x$log_or_true)
    is_null <- abs(truth) < 1e-6
    summarise_metrics(.x, log_or_true = truth, alpha = alpha, is_null = is_null)
  }) %>%
  ungroup()

write_tsv(
  metric_results,
  file.path(output_root, "tables", paste0("door_comparison_metrics", tag_suffix, ".tsv"))
)

# =====================================================================
# 9. Bias summaries
# =====================================================================

bias_results <- raw_results %>%
  filter(group != "A") %>%
  mutate(
    raw_bias = log_or_hat - log_or_true,
    prop_bias = raw_bias / abs(log_or_true)
  ) %>%
  group_by(sample_size, accuracy, group, model_type) %>%
  summarise(
    median_prop_bias = median(prop_bias, na.rm = TRUE),
    mean_prop_bias = mean(prop_bias, na.rm = TRUE),
    median_raw_bias = median(raw_bias, na.rm = TRUE),
    rmse_log = sqrt(mean(raw_bias^2, na.rm = TRUE)),
    .groups = "drop"
  )

write_tsv(
  bias_results,
  file.path(output_root, "tables", paste0("door_comparison_bias", tag_suffix, ".tsv"))
)

# =====================================================================
# 10. Plots
# =====================================================================

model_colours <- c(
  "Ordinal (polr)" = "#2166ac",
  "Binary (death)" = "#d6604d",
  "DOOR" = "#4dac26"
)

accuracy_labels <- c(`0.7` = "Accuracy = 70%", `1` = "Accuracy = 100%")

p_power <- metric_results %>%
  filter(group != "A") %>%
  ggplot(aes(x = sample_size, y = power, colour = model_type)) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(
    group ~ accuracy,
    scales = "free_y",
    labeller = labeller(accuracy = accuracy_labels)
  ) +
  scale_colour_manual(values = model_colours) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title = "Power: ordinal (polr) vs binary (death) vs DOOR",
    subtitle = "Same 6-point ordinal DGM; each method fitted to same simulated data",
    x = "Total trial size",
    y = "Power",
    colour = "Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_root, "plots", paste0("door_comparison_power", tag_suffix, ".pdf")),
  p_power, width = 12, height = 9
)

p_bias <- bias_results %>%
  ggplot(aes(x = sample_size, y = median_prop_bias, colour = model_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(
    group ~ accuracy,
    scales = "free_y",
    labeller = labeller(accuracy = accuracy_labels)
  ) +
  scale_colour_manual(values = model_colours) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Proportional bias: (estimate - truth) / |truth|",
    subtitle = "Each method assessed against its own estimand; scale comparable across methods",
    x = "Total trial size",
    y = "Median proportional bias",
    colour = "Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_root, "plots", paste0("door_comparison_bias", tag_suffix, ".pdf")),
  p_bias, width = 12, height = 9
)

p_dist <- raw_results %>%
  filter(sample_size == max(sample_sizes), group != "A") %>%
  mutate(
    prop_bias = (log_or_hat - log_or_true) / abs(log_or_true),
    accuracy_label = accuracy_labels[as.character(accuracy)]
  ) %>%
  ggplot(aes(x = model_type, y = prop_bias, fill = model_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_grid(group ~ accuracy_label, scales = "free_y") +
  scale_fill_manual(values = model_colours) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = paste0("Replicate distribution of proportional bias (N = ", max(sample_sizes), ")"),
    subtitle = "Whiskers = 1.5 x IQR",
    x = NULL,
    y = "Proportional bias",
    fill = "Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave(
  file.path(output_root, "plots", paste0("door_comparison_dist", tag_suffix, ".pdf")),
  p_dist, width = 12, height = 9
)

estimand_tbl <- tibble(
  group = names(or_vector),
  log_or_polr = log(or_vector),
  log_wo_door = true_log_win_odds
) %>%
  pivot_longer(
    cols = c(log_or_polr, log_wo_door),
    names_to = "estimand",
    values_to = "value"
  ) %>%
  mutate(
    estimand = recode(estimand,
      log_or_polr = "log-OR (polr / binary)",
      log_wo_door = "log-win-odds (DOOR)"
    )
  )

p_estimand <- estimand_tbl %>%
  ggplot(aes(x = group, y = value, fill = estimand)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5) +
  geom_hline(yintercept = 0, colour = "grey30") +
  scale_fill_manual(values = c(
    "log-OR (polr / binary)" = "#2166ac",
    "log-win-odds (DOOR)" = "#4dac26"
  )) +
  labs(
    title = "True estimand by subphenotype and method",
    subtitle = "DOOR win-odds shrunk toward null for rare-event subphenotypes (B, E)",
    x = "Subphenotype",
    y = "True log-effect",
    fill = "Estimand"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_root, "plots", paste0("door_comparison_estimands", tag_suffix, ".pdf")),
  p_estimand, width = 8, height = 5
)

message("Done. Results in: ", output_root)
