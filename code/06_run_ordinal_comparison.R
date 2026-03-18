## =====================================================================
##  06_run_ordinal_comparison.R -- Supplement 3 ordinal comparison
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

ordinal_baseline <- params$ordinal_baseline$prevalence
ordinal_labels <- params$ordinal_baseline$label

accuracy <- as.numeric(Sys.getenv("ACCURACY", unset = "0.80"))
sample_sizes <- params$ordinal_sample_sizes
sample_override <- Sys.getenv("ORDINAL_SAMPLE_SIZES", unset = "")
if (nzchar(sample_override)) {
  sample_sizes <- as.integer(strsplit(sample_override, ",")[[1]])
} else {
  sample_override <- Sys.getenv("SAMPLE_SIZES", unset = "")
  if (nzchar(sample_override)) {
    sample_sizes <- as.integer(strsplit(sample_override, ",")[[1]])
  }
}

n_reps <- as.integer(Sys.getenv("N_REPS_ORDINAL", unset = "1000"))
seed_base <- as.integer(Sys.getenv("SEED_BASE", unset = "321"))
run_order <- Sys.getenv("RUN_ORDER", unset = "binary,ordinal")
run_parts <- strsplit(run_order, ",")[[1]]
run_tag <- Sys.getenv("RUN_TAG", unset = "")

output_root <- file.path("results", "supp")
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "data"), showWarnings = FALSE, recursive = TRUE)
acc_label <- sprintf("acc%03d", as.integer(round(accuracy * 100)))
tag_suffix <- if (nzchar(run_tag)) paste0("_", run_tag) else ""

get_expected_ordinal_probs <- function(or_value, baseline_probs) {
  cumulative_probs <- cumsum(baseline_probs)[-length(baseline_probs)]
  intercepts <- qlogis(cumulative_probs)
  eta_control <- intercepts
  eta_treated <- intercepts + log(or_value)

  to_probs <- function(eta) {
    cumulative <- plogis(eta)
    probs <- diff(c(0, cumulative, 1))
    probs <- pmax(probs, 0)
    probs / sum(probs)
  }

  tibble(
    level = seq_along(baseline_probs),
    prob_control = to_probs(eta_control),
    prob_treated = to_probs(eta_treated)
  )
}

run_binary_scenario <- function(n, n_reps, seed_base) {
  groups <- names(or_vector)
  map_dfr(seq_len(n_reps), function(rep_id) {
    sim_data <- simulate_binary_trial(or_vector, prevalence, baseline, n, seed = seed_base + rep_id)
    sim_data$observed_group <- misclassify_groups(sim_data$group, prevalence, accuracy, seed = seed_base + rep_id + 10000)

    map_dfr(groups, function(group_label) {
      fit <- fit_logistic_or(sim_data %>% filter(observed_group == group_label))
      fit %>% mutate(group = group_label, rep_id = rep_id)
    })
  })
}

run_ordinal_scenario <- function(n, n_reps, seed_base) {
  groups <- names(or_vector)
  map_dfr(seq_len(n_reps), function(rep_id) {
    sim_data <- simulate_ordinal_trial(or_vector, prevalence, ordinal_baseline, n, seed = seed_base + rep_id)
    sim_data$observed_group <- misclassify_groups(sim_data$group, prevalence, accuracy, seed = seed_base + rep_id + 20000)

    map_dfr(groups, function(group_label) {
      fit <- fit_polr_or(sim_data %>% filter(observed_group == group_label))
      fit %>% mutate(group = group_label, rep_id = rep_id)
    })
  })
}

binary_results <- tibble()
ordinal_results <- tibble()

if ("binary" %in% run_parts) {
  binary_results <- map_dfr(sample_sizes, function(sample_size) {
    message(sprintf("Binary comparison N=%d", sample_size))
    run_binary_scenario(sample_size, n_reps, seed_base) %>%
      mutate(sample_size = sample_size, outcome_type = "Binary")
  })
}

if ("ordinal" %in% run_parts) {
  ordinal_results <- map_dfr(sample_sizes, function(sample_size) {
    message(sprintf("Ordinal comparison N=%d", sample_size))
    run_ordinal_scenario(sample_size, n_reps, seed_base) %>%
      mutate(sample_size = sample_size, outcome_type = "Ordinal")
  })
}

combined_results <- bind_rows(binary_results, ordinal_results)

raw_results <- combined_results %>%
  mutate(
    log_or_true = log(or_vector[group]),
    bias_log_or = log_or_hat - log_or_true,
    abs_bias_log_or = abs(bias_log_or),
    accuracy = accuracy
  )

raw_path <- file.path(output_root, "data", paste0("supp3_binary_vs_ordinal_raw_", acc_label, tag_suffix, ".tsv.gz"))
write_tsv(raw_results, raw_path)

metric_results <- raw_results %>%
  group_by(sample_size, outcome_type, group) %>%
  group_modify(~ {
    log_or_true <- log(or_vector[[.y$group]])
    is_null <- abs(log_or_true) < 1e-6
    summarise_metrics(.x, log_or_true, params$alpha_primary, is_null = is_null)
  }) %>%
  ungroup() %>%
  mutate(accuracy = accuracy)

metrics_path <- file.path(output_root, "tables", paste0("supp3_binary_vs_ordinal_metrics_", acc_label, tag_suffix, ".tsv"))
write_tsv(metric_results, metrics_path)

table_s4 <- metric_results %>%
  filter(sample_size == 5000, group != "A") %>%
  select(sample_size, outcome_type, group, power, type_s, type_m)

table_path <- file.path(output_root, "tables", paste0("supp3_power_binary_vs_ordinal_N5000_", acc_label, tag_suffix, ".tsv"))
write_tsv(table_s4, table_path)

power_gain <- metric_results %>%
  filter(group != "A") %>%
  select(sample_size, group, outcome_type, power) %>%
  tidyr::pivot_wider(names_from = outcome_type, values_from = power) %>%
  mutate(
    power_gain_abs = Ordinal - Binary,
    power_gain_rel = if_else(Binary > 0, (Ordinal / Binary) - 1, NA_real_)
  )

gain_path <- file.path(output_root, "tables", paste0("supp3_power_gain_", acc_label, tag_suffix, ".tsv"))
write_tsv(power_gain, gain_path)

type1_results <- metric_results %>%
  filter(group == "A") %>%
  select(sample_size, outcome_type, type1, accuracy)

type1_path <- file.path(output_root, "tables", paste0("supp3_type1_binary_vs_ordinal_", acc_label, tag_suffix, ".tsv"))
write_tsv(type1_results, type1_path)

type1_plot <- type1_results %>%
  ggplot(aes(x = sample_size, y = type1, color = outcome_type)) +
  geom_line() +
  geom_point(size = 1.5) +
  geom_hline(yintercept = params$alpha_primary, linetype = "dashed", color = "grey50") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Supplement 3: Type I error (null subgroup A)",
    x = "Sample size",
    y = "Type I error",
    color = "Outcome"
  ) +
  theme_minimal()

type1_plot_path <- file.path(output_root, "plots", paste0("figureS5_type1_binary_vs_ordinal_", acc_label, tag_suffix, ".pdf"))
ggsave(
  type1_plot_path,
  type1_plot,
  width = 9,
  height = 6
)

bias_summary <- raw_results %>%
  filter(group != "A") %>%
  group_by(sample_size, outcome_type, group) %>%
  summarise(
    mean_bias_log_or = mean(bias_log_or, na.rm = TRUE),
    mean_abs_bias_log_or = mean(abs_bias_log_or, na.rm = TRUE),
    rmse_log_or = sqrt(mean(bias_log_or^2, na.rm = TRUE)),
    .groups = "drop"
  )

bias_summary_path <- file.path(output_root, "tables", paste0("supp3_bias_binary_vs_ordinal_", acc_label, tag_suffix, ".tsv"))
write_tsv(bias_summary, bias_summary_path)

bias_plot <- bias_summary %>%
  ggplot(aes(x = sample_size, y = mean_abs_bias_log_or, color = outcome_type)) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_wrap(~ group, scales = "free_y") +
  labs(
    title = "Supplement 3: Mean absolute bias (binary vs ordinal)",
    x = "Sample size",
    y = "Mean absolute bias in log-OR",
    color = "Outcome"
  ) +
  theme_minimal()

bias_plot_path <- file.path(output_root, "plots", paste0("figureS5_bias_binary_vs_ordinal_", acc_label, tag_suffix, ".pdf"))
ggsave(
  bias_plot_path,
  bias_plot,
  width = 11,
  height = 7
)

bias_wide <- bias_summary %>%
  select(sample_size, group, outcome_type, mean_abs_bias_log_or, mean_bias_log_or) %>%
  tidyr::pivot_wider(
    names_from = outcome_type,
    values_from = c(mean_abs_bias_log_or, mean_bias_log_or)
  ) %>%
  mutate(
    delta_abs_bias = mean_abs_bias_log_or_Ordinal - mean_abs_bias_log_or_Binary,
    delta_signed_bias = mean_bias_log_or_Ordinal - mean_bias_log_or_Binary
  )

power_bias_tradeoff <- power_gain %>%
  left_join(bias_wide, by = c("sample_size", "group"))

tradeoff_path <- file.path(output_root, "tables", paste0("supp3_power_vs_bias_tradeoff_", acc_label, tag_suffix, ".tsv"))
write_tsv(power_bias_tradeoff, tradeoff_path)

tradeoff_plot <- power_bias_tradeoff %>%
  ggplot(aes(x = power_gain_abs, y = delta_abs_bias, color = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(alpha = 0.8, size = 2) +
  labs(
    title = "Supplement 3: Power gain vs change in absolute bias (ordinal - binary)",
    x = "Absolute power gain (Ordinal - Binary)",
    y = "Change in absolute bias (Ordinal - Binary)",
    color = "Subphenotype"
  ) +
  theme_minimal()

tradeoff_plot_path <- file.path(output_root, "plots", paste0("figureS5_power_vs_bias_tradeoff_", acc_label, tag_suffix, ".pdf"))
ggsave(
  tradeoff_plot_path,
  tradeoff_plot,
  width = 9,
  height = 6
)

power_plot <- metric_results %>%
  filter(group != "A") %>%
  ggplot(aes(x = sample_size, y = power, color = outcome_type)) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_wrap(~ group, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = sprintf("Supplement 3: Power curves (binary vs ordinal, %.0f%% accuracy)", accuracy * 100),
    x = "Sample size",
    y = "Power",
    color = "Outcome"
  ) +
  theme_minimal()

plot_path <- file.path(output_root, "plots", paste0("figureS3_power_binary_vs_ordinal_", acc_label, tag_suffix, ".pdf"))
ggsave(
  plot_path,
  power_plot,
  width = 11,
  height = 7
)

ordinal_shift <- map_dfr(names(or_vector), function(group_label) {
  probs <- get_expected_ordinal_probs(or_vector[[group_label]], ordinal_baseline)
  bind_rows(
    probs %>%
      transmute(
        group = group_label,
        level,
        label = ordinal_labels[level],
        treatment = "Control",
        probability = prob_control
      ),
    probs %>%
      transmute(
        group = group_label,
        level,
        label = ordinal_labels[level],
        treatment = "Treatment",
        probability = prob_treated
      )
  )
}) %>%
  mutate(
    group_label = subphenotypes$label[match(group, subphenotypes$subphenotype)],
    label = factor(label, levels = ordinal_labels)
  )

shift_path <- file.path(output_root, "tables", paste0("supp3_ordinal_shift_expected_", acc_label, tag_suffix, ".tsv"))
write_tsv(ordinal_shift, shift_path)

shift_plot <- ordinal_shift %>%
  ggplot(aes(x = treatment, y = probability, fill = label)) +
  geom_col(position = "fill", width = 0.65) +
  facet_wrap(~ group_label, ncol = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = sprintf("Supplement 3: Expected ordinal outcome distribution shift (%.0f%% accuracy context)", accuracy * 100),
    x = "",
    y = "Category proportion",
    fill = "Ordinal level"
  ) +
  theme_minimal()

shift_plot_path <- file.path(output_root, "plots", paste0("figureS4_ordinal_shift_stacked_", acc_label, tag_suffix, ".pdf"))
ggsave(
  shift_plot_path,
  shift_plot,
  width = 12,
  height = 8
)

shift_delta <- ordinal_shift %>%
  select(group, group_label, level, label, treatment, probability) %>%
  tidyr::pivot_wider(names_from = treatment, values_from = probability) %>%
  mutate(delta_treatment_minus_control = Treatment - Control)

shift_delta_path <- file.path(output_root, "tables", paste0("supp3_ordinal_shift_delta_", acc_label, tag_suffix, ".tsv"))
write_tsv(shift_delta, shift_delta_path)

delta_plot <- shift_delta %>%
  ggplot(aes(x = label, y = delta_treatment_minus_control, color = group_label, group = group_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line() +
  geom_point(size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Supplement 3: Treatment-driven shift by ordinal level",
    x = "Ordinal level",
    y = "Absolute change in probability",
    color = "Subphenotype"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

delta_plot_path <- file.path(output_root, "plots", paste0("figureS4_ordinal_shift_delta_", acc_label, tag_suffix, ".pdf"))
ggsave(
  delta_plot_path,
  delta_plot,
  width = 11,
  height = 7
)

# Keep legacy file names for 80% accuracy run to preserve downstream links.
if (abs(accuracy - 0.80) < 1e-8 && !nzchar(run_tag)) {
  write_tsv(raw_results, file.path(output_root, "data", "supp3_binary_vs_ordinal_raw.tsv.gz"))
  write_tsv(metric_results, file.path(output_root, "tables", "supp3_binary_vs_ordinal_metrics.tsv"))
  write_tsv(table_s4, file.path(output_root, "tables", "supp3_power_binary_vs_ordinal_N5000.tsv"))
  write_tsv(type1_results, file.path(output_root, "tables", "supp3_type1_binary_vs_ordinal.tsv"))
  write_tsv(bias_summary, file.path(output_root, "tables", "supp3_bias_binary_vs_ordinal.tsv"))
  write_tsv(power_bias_tradeoff, file.path(output_root, "tables", "supp3_power_vs_bias_tradeoff.tsv"))
  write_tsv(ordinal_shift, file.path(output_root, "tables", "supp3_ordinal_shift_expected.tsv"))
  write_tsv(shift_delta, file.path(output_root, "tables", "supp3_ordinal_shift_delta.tsv"))
  ggsave(
    file.path(output_root, "plots", "figureS3_power_binary_vs_ordinal.pdf"),
    power_plot,
    width = 11,
    height = 7
  )
  ggsave(
    file.path(output_root, "plots", "figureS5_type1_binary_vs_ordinal.pdf"),
    type1_plot,
    width = 9,
    height = 6
  )
  ggsave(
    file.path(output_root, "plots", "figureS5_bias_binary_vs_ordinal.pdf"),
    bias_plot,
    width = 11,
    height = 7
  )
  ggsave(
    file.path(output_root, "plots", "figureS5_power_vs_bias_tradeoff.pdf"),
    tradeoff_plot,
    width = 9,
    height = 6
  )
  ggsave(
    file.path(output_root, "plots", "figureS4_ordinal_shift_stacked.pdf"),
    shift_plot,
    width = 12,
    height = 8
  )
  ggsave(
    file.path(output_root, "plots", "figureS4_ordinal_shift_delta.pdf"),
    delta_plot,
    width = 11,
    height = 7
  )
}
