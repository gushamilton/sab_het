## =====================================================================
##  11_run_figure1_accuracy_scan.R -- Core Figure 1 (accuracy scan)
## =====================================================================

pacman::p_load(tidyverse, patchwork, ggdist)

source("code/01_parameters.R")
source("code/functions/01_simulation_helpers.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table

or_vector <- setNames(subphenotypes$or_arrest_shrunk, subphenotypes$subphenotype)
prevalence <- setNames(subphenotypes$prevalence, subphenotypes$subphenotype)
baseline <- setNames(subphenotypes$baseline_mortality, subphenotypes$subphenotype)

accuracy_grid <- params$accuracy_grid
n_total <- as.integer(Sys.getenv("FIG1_N", unset = "20000"))
n_reps <- as.integer(Sys.getenv("FIG1_N_REPS", unset = "1000"))
seed_base <- as.integer(Sys.getenv("SEED_BASE", unset = "9001"))
alpha <- params$alpha_primary
z_alpha <- qnorm(1 - alpha / 2)
trim_outliers <- tolower(Sys.getenv("FIG1_TRIM_OUTLIERS", unset = "0")) %in% c("1", "true", "yes")
trim_prob <- as.numeric(Sys.getenv("FIG1_TRIM_PROB", unset = "0.01"))

output_root <- "results/core_figures"
dir.create(file.path(output_root, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)

group_palette <- c(
  A = "#F8766D",
  B = "#A3A500",
  C = "#00BF7D",
  D = "#00B0F6",
  E = "#E76BF3"
)

simulate_accuracy_slice <- function(acc) {
  groups <- names(or_vector)
  map_dfr(seq_len(n_reps), function(rep_id) {
    sim_data <- simulate_binary_trial(
      or_vector = or_vector,
      prevalence = prevalence,
      baseline_mortality = baseline,
      n = n_total,
      seed = seed_base + rep_id
    )
    sim_data$observed_group <- misclassify_groups(
      true_group = sim_data$group,
      prevalence = prevalence,
      accuracy = acc,
      seed = seed_base + rep_id + 10000
    )

    map_dfr(groups, function(group_label) {
      fit <- fit_logistic_or(sim_data %>% filter(observed_group == group_label))
      fit %>%
        mutate(
          group = group_label,
          rep_id = rep_id,
          accuracy = acc
        )
    })
  })
}

raw_results <- map_dfr(accuracy_grid, function(acc) {
  message(sprintf("Figure1 simulation: N=%d, accuracy=%.2f", n_total, acc))
  simulate_accuracy_slice(acc)
}) %>%
  mutate(
    log_or_true = log(or_vector[group]),
    ci_lower = log_or_hat - z_alpha * se,
    ci_upper = log_or_hat + z_alpha * se,
    significant = !is.na(log_or_hat) & !is.na(se) & (ci_lower > 0 | ci_upper < 0),
    power_flag = significant & (sign(log_or_hat) == sign(log_or_true)),
    bias_toward_null = log_or_true - log_or_hat
  )

write_tsv(
  raw_results,
  file.path(output_root, "data", "figure1_accuracy_scan_raw.tsv.gz")
)

summary_curve <- raw_results %>%
  group_by(accuracy, group) %>%
  summarise(
    power = mean(power_flag, na.rm = TRUE),
    bias_toward_null = mean(bias_toward_null, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  summary_curve,
  file.path(output_root, "data", "figure1_accuracy_scan_summary.tsv")
)

panel_a <- summary_curve %>%
  ggplot(aes(x = accuracy, y = power, color = group)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = group_palette) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = accuracy_grid) +
  labs(
    x = "Accuracy",
    y = "Power",
    color = "Group"
  ) +
  ggdist::theme_ggdist(base_size = 11)

panel_b <- summary_curve %>%
  ggplot(aes(x = accuracy, y = bias_toward_null, color = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#F8766D") +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = group_palette) +
  scale_x_continuous(breaks = accuracy_grid) +
  labs(
    x = "Accuracy",
    y = "Bias (log OR)",
    color = "Group"
  ) +
  ggdist::theme_ggdist(base_size = 11)

panel_c_data <- raw_results %>%
  mutate(accuracy_round = round(accuracy, 2)) %>%
  filter(accuracy_round %in% c(1.00, 0.90, 0.80, 0.70)) %>%
  mutate(
    accuracy_label = factor(
      paste0("Accuracy = ", sprintf("%.1f", accuracy_round)),
      levels = c("Accuracy = 1.0", "Accuracy = 0.9", "Accuracy = 0.8", "Accuracy = 0.7")
    ),
    group = factor(group, levels = c("A", "B", "C", "D", "E"))
  ) %>%
  select(-accuracy_round)

if (trim_outliers) {
  panel_c_before <- nrow(panel_c_data)
  panel_c_data <- panel_c_data %>%
    group_by(accuracy_label, group) %>%
    mutate(
      q_low = quantile(log_or_hat, probs = trim_prob, na.rm = TRUE),
      q_high = quantile(log_or_hat, probs = 1 - trim_prob, na.rm = TRUE)
    ) %>%
    filter(log_or_hat >= q_low, log_or_hat <= q_high) %>%
    ungroup() %>%
    select(-q_low, -q_high)
  message(sprintf(
    "Panel C outlier trim enabled: removed %d points (%.2f%%)",
    panel_c_before - nrow(panel_c_data),
    100 * (panel_c_before - nrow(panel_c_data)) / panel_c_before
  ))
}

true_effect_bars <- tidyr::expand_grid(
  accuracy_label = levels(panel_c_data$accuracy_label),
  group = factor(c("A", "B", "C", "D", "E"), levels = c("A", "B", "C", "D", "E"))
) %>%
  mutate(log_or_true = log(or_vector[as.character(group)]))

panel_c <- panel_c_data %>%
  ggplot(aes(x = group, y = log_or_hat, color = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#F8766D") +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.25) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 0.7) +
  geom_crossbar(
    data = true_effect_bars,
    aes(x = group, y = log_or_true, ymin = log_or_true, ymax = log_or_true),
    inherit.aes = FALSE,
    color = "black",
    width = 0.55,
    linewidth = 0.45
  ) +
  scale_color_manual(values = group_palette) +
  facet_wrap(~ accuracy_label, ncol = 2) +
  labs(
    x = "Subphenotype",
    y = "Beta (log OR)",
    color = "Group"
  ) +
  ggdist::theme_ggdist(base_size = 11)

figure_1 <- ((panel_a + panel_b) / panel_c) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", heights = c(1, 1.6)) &
  theme(legend.position = "bottom")

ggsave(
  file.path(output_root, "plots", "figure1_accuracy_dilution_patchwork.pdf"),
  figure_1,
  width = 11,
  height = 12
)
