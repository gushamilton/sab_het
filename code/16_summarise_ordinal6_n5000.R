## =====================================================================
##  16_summarise_ordinal6_n5000.R
##  Focused summaries/plots for ordinal6 sensitivity at N=5000.
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")
source("code/functions/02_metrics.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table
or_vector <- setNames(subphenotypes$or_arrest_shrunk, subphenotypes$subphenotype)
alpha <- params$alpha_primary

output_root <- "results/supp/ordinal6_sameeffect"
raw_path <- file.path(output_root, "data", "ordinal6_sameeffect_raw.tsv.gz")
if (!file.exists(raw_path)) {
  stop("Missing raw file: ", raw_path)
}

raw_results <- read_tsv(raw_path, show_col_types = FALSE) %>%
  filter(sample_size == 5000) %>%
  mutate(
    group = factor(group, levels = subphenotypes$subphenotype),
    model_type = factor(model_type, levels = c("Binary", "Ordinal"))
  )

if (nrow(raw_results) == 0) {
  stop("No rows for sample_size == 5000 in raw results.")
}

metrics <- raw_results %>%
  group_by(sample_size, accuracy, group, model_type) %>%
  group_modify(~ {
    g <- as.character(.y$group)
    truth <- log(or_vector[[g]])
    is_null <- abs(truth) < 1e-6
    metric_row <- summarise_metrics(.x, log_or_true = truth, alpha = alpha, is_null = is_null)

    tibble(
      type1 = metric_row$type1,
      power = metric_row$power,
      type_s = metric_row$type_s,
      type_m = metric_row$type_m,
      median_signed_bias_log_or = median(.x$log_or_hat - .x$log_or_true, na.rm = TRUE),
      median_abs_error_log_or = median(abs(.x$log_or_hat - .x$log_or_true), na.rm = TRUE),
      rmse_log_or = sqrt(mean((.x$log_or_hat - .x$log_or_true)^2, na.rm = TRUE))
    )
  }) %>%
  ungroup() %>%
  mutate(group = as.character(group))

write_tsv(metrics, file.path(output_root, "tables", "ordinal6_sameeffect_metrics_n5000_focus.tsv"))

truth_df <- subphenotypes %>%
  transmute(
    group = factor(subphenotype, levels = subphenotypes$subphenotype),
    log_or_true = log(or_arrest_shrunk)
  )

replicate_plot <- raw_results %>%
  ggplot(aes(x = group, y = log_or_hat, color = model_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_jitter(width = 0.16, height = 0, alpha = 0.18, size = 0.7, stroke = 0) +
  geom_errorbar(
    data = truth_df,
    aes(x = group, ymin = log_or_true, ymax = log_or_true),
    inherit.aes = FALSE,
    width = 0.55,
    linewidth = 0.9,
    color = "black"
  ) +
  facet_grid(model_type ~ accuracy) +
  labs(
    title = "N=5,000: replicate log(OR) estimates vs true subgroup effect",
    subtitle = "Black bars show true log(OR); points are replicate estimates",
    x = "Subphenotype",
    y = "Estimated log(OR)",
    color = "Model"
  ) +
  theme_minimal()

ggsave(
  file.path(output_root, "plots", "ordinal6_sameeffect_replicates_n5000.pdf"),
  replicate_plot,
  width = 12,
  height = 8
)

metric_long <- metrics %>%
  mutate(group = factor(group, levels = subphenotypes$subphenotype)) %>%
  select(sample_size, accuracy, group, model_type, power, type_s, type1, median_signed_bias_log_or) %>%
  pivot_longer(
    cols = c(power, type_s, type1, median_signed_bias_log_or),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      power = "Power",
      type_s = "Type S",
      type1 = "Type I (A only)",
      median_signed_bias_log_or = "Median Signed Bias"
    )
  )

metric_plot <- metric_long %>%
  ggplot(aes(x = group, y = value, color = model_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65") +
  geom_point(position = position_dodge(width = 0.35), size = 2) +
  facet_grid(metric ~ accuracy, scales = "free_y") +
  labs(
    title = "N=5,000: key operating characteristics",
    x = "Subphenotype",
    y = "Value",
    color = "Model"
  ) +
  theme_minimal()

ggsave(
  file.path(output_root, "plots", "ordinal6_sameeffect_key_metrics_n5000.pdf"),
  metric_plot,
  width = 12,
  height = 10
)

power_scatter <- metrics %>%
  filter(group != "A") %>%
  select(accuracy, group, model_type, power) %>%
  pivot_wider(names_from = model_type, values_from = power) %>%
  ggplot(aes(x = Binary, y = Ordinal, color = group, shape = factor(accuracy))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.8, alpha = 0.9) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "N=5,000: Power gain from ordinal vs binary model",
    subtitle = "Points above diagonal indicate higher power for ordinal",
    x = "Binary power",
    y = "Ordinal power",
    color = "Subphenotype",
    shape = "Accuracy"
  ) +
  theme_minimal()

ggsave(
  file.path(output_root, "plots", "ordinal6_sameeffect_power_scatter_n5000.pdf"),
  power_scatter,
  width = 8,
  height = 7
)

bias_style_plot <- raw_results %>%
  filter(group != "A") %>%
  ggplot(aes(x = group, y = log_or_hat, color = model_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_jitter(width = 0.16, height = 0, alpha = 0.20, size = 0.75, stroke = 0) +
  geom_errorbar(
    data = truth_df %>% filter(group != "A"),
    aes(x = group, ymin = log_or_true, ymax = log_or_true),
    inherit.aes = FALSE,
    width = 0.58,
    linewidth = 1.0,
    color = "black"
  ) +
  facet_wrap(~ accuracy, ncol = 2) +
  labs(
    title = "N=5,000: replicate estimates vs true subgroup log(OR)",
    subtitle = "Black horizontal bars are true subgroup log(OR)",
    x = "Subphenotype",
    y = "Estimated log(OR)",
    color = "Model"
  ) +
  theme_minimal()

ggsave(
  file.path(output_root, "plots", "ordinal6_sameeffect_bias_replicates_n5000.pdf"),
  bias_style_plot,
  width = 11,
  height = 7
)

message("Wrote focused N=5000 summaries and plots.")
