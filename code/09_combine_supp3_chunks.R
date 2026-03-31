## =====================================================================
##  09_combine_supp3_chunks.R -- Combine chunked Supplement 3 outputs
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table
ordinal_points <- as.integer(Sys.getenv("ORDINAL_POINTS", unset = "6"))
if (ordinal_points == 6) {
  ordinal_labels <- c(
    "Dead",
    "ICU/ventilated",
    "Still hospitalised",
    "Discharged to rehab",
    "Discharged with complications",
    "Discharged well"
  )
} else if (ordinal_points == 5) {
  ordinal_labels <- params$ordinal_baseline$label
} else {
  stop("ORDINAL_POINTS must be 5 or 6.", call. = FALSE)
}

accuracy <- as.numeric(Sys.getenv("ACCURACY", unset = "1.00"))
acc_label <- sprintf("acc%03d", as.integer(round(accuracy * 100)))
chunk_tags <- Sys.getenv("CHUNK_TAGS", unset = "")
chunk_count <- as.integer(Sys.getenv("CHUNK_COUNT", unset = "0"))

if (nzchar(chunk_tags)) {
  tags <- strsplit(chunk_tags, ",")[[1]]
} else if (!is.na(chunk_count) && chunk_count > 0) {
  tags <- paste0("chunk", seq_len(chunk_count))
} else {
  stop("Set CHUNK_TAGS or CHUNK_COUNT for chunked Supplement 3 combine.")
}

output_root <- file.path("results", "supp")
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "data"), showWarnings = FALSE, recursive = TRUE)

read_chunk <- function(prefix, tag, ext = "tsv") {
  path <- file.path(output_root, if (prefix == "supp3_binary_vs_ordinal_raw") "data" else "tables",
                    paste0(prefix, "_", acc_label, "_", tag, ".", ext))
  if (!file.exists(path)) stop("Missing chunk file: ", path)
  readr::read_tsv(path, show_col_types = FALSE)
}

raw_all <- map_dfr(tags, ~ read_chunk("supp3_binary_vs_ordinal_raw", .x, "tsv.gz")) %>%
  distinct(sample_size, outcome_type, group, rep_id, .keep_all = TRUE)

metric_all <- map_dfr(tags, ~ read_chunk("supp3_binary_vs_ordinal_metrics", .x, "tsv")) %>%
  distinct(sample_size, outcome_type, group, .keep_all = TRUE) %>%
  arrange(sample_size, outcome_type, group)

write_tsv(raw_all, file.path(output_root, "data", paste0("supp3_binary_vs_ordinal_raw_", acc_label, ".tsv.gz")))
write_tsv(metric_all, file.path(output_root, "tables", paste0("supp3_binary_vs_ordinal_metrics_", acc_label, ".tsv")))

table_s4 <- metric_all %>%
  filter(sample_size == 5000, group != "A") %>%
  select(sample_size, outcome_type, group, power, type_s, type_m)
write_tsv(table_s4, file.path(output_root, "tables", paste0("supp3_power_binary_vs_ordinal_N5000_", acc_label, ".tsv")))

power_gain <- metric_all %>%
  filter(group != "A") %>%
  select(sample_size, group, outcome_type, power) %>%
  tidyr::pivot_wider(names_from = outcome_type, values_from = power) %>%
  mutate(
    power_gain_abs = Ordinal - Binary,
    power_gain_rel = if_else(Binary > 0, (Ordinal / Binary) - 1, NA_real_)
  )
write_tsv(power_gain, file.path(output_root, "tables", paste0("supp3_power_gain_", acc_label, ".tsv")))

type1_results <- metric_all %>%
  filter(group == "A") %>%
  select(sample_size, outcome_type, type1, accuracy)
write_tsv(type1_results, file.path(output_root, "tables", paste0("supp3_type1_binary_vs_ordinal_", acc_label, ".tsv")))

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
ggsave(file.path(output_root, "plots", paste0("figureS5_type1_binary_vs_ordinal_", acc_label, ".pdf")), type1_plot, width = 9, height = 6)

bias_summary <- raw_all %>%
  filter(group != "A") %>%
  group_by(sample_size, outcome_type, group) %>%
  summarise(
    mean_bias_log_or = mean(bias_log_or, na.rm = TRUE),
    mean_abs_bias_log_or = mean(abs_bias_log_or, na.rm = TRUE),
    rmse_log_or = sqrt(mean(bias_log_or^2, na.rm = TRUE)),
    .groups = "drop"
  )
write_tsv(bias_summary, file.path(output_root, "tables", paste0("supp3_bias_binary_vs_ordinal_", acc_label, ".tsv")))

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
ggsave(file.path(output_root, "plots", paste0("figureS5_bias_binary_vs_ordinal_", acc_label, ".pdf")), bias_plot, width = 11, height = 7)

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
write_tsv(power_bias_tradeoff, file.path(output_root, "tables", paste0("supp3_power_vs_bias_tradeoff_", acc_label, ".tsv")))

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
ggsave(file.path(output_root, "plots", paste0("figureS5_power_vs_bias_tradeoff_", acc_label, ".pdf")), tradeoff_plot, width = 9, height = 6)

power_plot <- metric_all %>%
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
ggsave(file.path(output_root, "plots", paste0("figureS3_power_binary_vs_ordinal_", acc_label, ".pdf")), power_plot, width = 11, height = 7)

ordinal_shift <- read_chunk("supp3_ordinal_shift_expected", tags[[1]], "tsv") %>%
  mutate(label = factor(label, levels = ordinal_labels)) %>%
  mutate(
    group_label = subphenotypes$label[match(group, subphenotypes$subphenotype)],
    label = factor(label, levels = ordinal_labels)
  )
write_tsv(ordinal_shift, file.path(output_root, "tables", paste0("supp3_ordinal_shift_expected_", acc_label, ".tsv")))

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
ggsave(file.path(output_root, "plots", paste0("figureS4_ordinal_shift_stacked_", acc_label, ".pdf")), shift_plot, width = 12, height = 8)

shift_delta <- ordinal_shift %>%
  select(group, group_label, level, label, treatment, probability) %>%
  tidyr::pivot_wider(names_from = treatment, values_from = probability) %>%
  mutate(delta_treatment_minus_control = Treatment - Control)
write_tsv(shift_delta, file.path(output_root, "tables", paste0("supp3_ordinal_shift_delta_", acc_label, ".tsv")))

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
ggsave(file.path(output_root, "plots", paste0("figureS4_ordinal_shift_delta_", acc_label, ".pdf")), delta_plot, width = 11, height = 7)
