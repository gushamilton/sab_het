## =====================================================================
##  26_build_core_figure4_ordinal.R -- Core Figure 4 (ordinal analysis)
## =====================================================================

library(tidyverse)
library(patchwork)

accuracy <- as.numeric(Sys.getenv("ACCURACY", unset = "1.00"))
acc_label <- sprintf("acc%03d", as.integer(round(accuracy * 100)))
run_tag <- Sys.getenv("RUN_TAG", unset = "")
tag_suffix <- if (nzchar(run_tag)) paste0("_", run_tag) else ""

tables_dir <- file.path("results", "supp", "tables")
plots_dir <- file.path("results", "core_figures", "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

metrics_path <- file.path(tables_dir, paste0("supp3_binary_vs_ordinal_metrics_", acc_label, tag_suffix, ".tsv"))
type1_path <- file.path(tables_dir, paste0("supp3_type1_binary_vs_ordinal_", acc_label, tag_suffix, ".tsv"))
bias_path <- file.path(tables_dir, paste0("supp3_bias_binary_vs_ordinal_", acc_label, tag_suffix, ".tsv"))
shift_path <- file.path(tables_dir, paste0("supp3_ordinal_shift_expected_", acc_label, tag_suffix, ".tsv"))

required <- c(metrics_path, type1_path, bias_path, shift_path)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop(
    "Cannot build Figure 4. Missing inputs:\n",
    paste0(" - ", missing, collapse = "\n"),
    call. = FALSE
  )
}

metrics <- readr::read_tsv(metrics_path, show_col_types = FALSE)
type1 <- readr::read_tsv(type1_path, show_col_types = FALSE)
bias <- readr::read_tsv(bias_path, show_col_types = FALSE)
shift <- readr::read_tsv(shift_path, show_col_types = FALSE)

panel_a <- metrics %>%
  filter(group != "A") %>%
  ggplot(aes(x = sample_size, y = power, color = outcome_type)) +
  geom_line() +
  geom_point(size = 1.4) +
  facet_wrap(~ group, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Power by sample size",
    x = "Sample size",
    y = "Power",
    color = "Outcome"
  ) +
  theme_minimal(base_size = 11)

panel_b <- type1 %>%
  ggplot(aes(x = sample_size, y = type1, color = outcome_type)) +
  geom_line() +
  geom_point(size = 1.4) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Type I error (null subgroup A)",
    x = "Sample size",
    y = "Type I error",
    color = "Outcome"
  ) +
  theme_minimal(base_size = 11)

panel_c <- bias %>%
  filter(group != "A") %>%
  ggplot(aes(x = sample_size, y = mean_abs_bias_log_or, color = outcome_type)) +
  geom_line() +
  geom_point(size = 1.4) +
  facet_wrap(~ group, scales = "free_y") +
  labs(
    title = "Mean absolute bias in log-OR",
    x = "Sample size",
    y = "Mean absolute bias",
    color = "Outcome"
  ) +
  theme_minimal(base_size = 11)

panel_d <- shift %>%
  mutate(
    treatment = factor(treatment, levels = c("Control", "Treatment"))
  ) %>%
  ggplot(aes(x = treatment, y = probability, fill = label)) +
  geom_col(position = "fill", width = 0.65) +
  facet_wrap(~ group_label, ncol = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Expected ordinal distribution shift",
    x = NULL,
    y = "Category proportion",
    fill = "Ordinal level"
  ) +
  theme_minimal(base_size = 11)

figure_4 <- (panel_a / panel_b / panel_c / panel_d) +
  plot_annotation(
    title = sprintf("Figure 4: Ordinal vs binary outcomes (%.0f%% accuracy)", accuracy * 100),
    tag_levels = "A"
  ) +
  plot_layout(guides = "collect", heights = c(1, 0.8, 1, 1.2)) &
  theme(legend.position = "bottom")

ggsave(
  file.path(plots_dir, "figure4.pdf"),
  figure_4,
  width = 12,
  height = 16
)
