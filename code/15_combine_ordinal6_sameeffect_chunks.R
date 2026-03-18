## =====================================================================
##  15_combine_ordinal6_sameeffect_chunks.R
##  Combine chunked outputs from 14_run_ordinal6_same_effect_comparison.R
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")
source("code/functions/02_metrics.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table
or_vector <- setNames(subphenotypes$or_arrest_shrunk, subphenotypes$subphenotype)
alpha <- params$alpha_primary

output_root <- "results/supp/ordinal6_sameeffect"
dir.create(file.path(output_root, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)

sample_override <- Sys.getenv("CHUNK_SAMPLE_SIZES", unset = "")
if (nzchar(sample_override)) {
  sample_sizes <- as.integer(strsplit(sample_override, ",")[[1]])
} else {
  sample_sizes <- params$sample_sizes
}

n_reps_total <- as.integer(Sys.getenv("N_REPS_TOTAL", unset = "1000"))
reps_per_chunk <- as.integer(Sys.getenv("REPS_PER_CHUNK", unset = "200"))

if (n_reps_total %% reps_per_chunk != 0) {
  stop("N_REPS_TOTAL must be divisible by REPS_PER_CHUNK.")
}
n_chunks <- n_reps_total / reps_per_chunk

chunk_tags <- expand_grid(sample_size = sample_sizes, chunk = seq_len(n_chunks)) %>%
  mutate(tag = paste0("n", sample_size, "_r", chunk)) %>%
  pull(tag)

read_chunk_raw <- function(tag) {
  path <- file.path(output_root, "data", paste0("ordinal6_sameeffect_raw_", tag, ".tsv.gz"))
  if (!file.exists(path)) stop("Missing chunk file: ", path)
  read_tsv(path, show_col_types = FALSE)
}

raw_all <- map_dfr(chunk_tags, read_chunk_raw) %>%
  distinct(sample_size, accuracy, group, model_type, rep_id, .keep_all = TRUE) %>%
  arrange(sample_size, accuracy, group, model_type, rep_id)

write_tsv(raw_all, file.path(output_root, "data", "ordinal6_sameeffect_raw.tsv.gz"))

metric_results <- raw_all %>%
  group_by(sample_size, accuracy, group, model_type) %>%
  group_modify(~ {
    truth <- log(or_vector[[.y$group]])
    is_null <- abs(truth) < 1e-6
    summarise_metrics(.x, log_or_true = truth, alpha = alpha, is_null = is_null)
  }) %>%
  ungroup()

write_tsv(metric_results, file.path(output_root, "tables", "ordinal6_sameeffect_metrics.tsv"))

bias_results <- raw_all %>%
  filter(group != "A") %>%
  group_by(sample_size, accuracy, group, model_type) %>%
  summarise(
    mean_bias_log_or = mean(log_or_hat - log_or_true, na.rm = TRUE),
    mean_abs_bias_log_or = mean(abs(log_or_hat - log_or_true), na.rm = TRUE),
    rmse_log_or = sqrt(mean((log_or_hat - log_or_true)^2, na.rm = TRUE)),
    .groups = "drop"
  )

write_tsv(bias_results, file.path(output_root, "tables", "ordinal6_sameeffect_bias.tsv"))

calibration_paths <- list.files(
  file.path(output_root, "tables"),
  pattern = "^ordinal6_sameeffect_calibration_.*\\.tsv$",
  full.names = TRUE
)
if (length(calibration_paths) > 0) {
  calibration <- map_dfr(calibration_paths, ~ read_tsv(.x, show_col_types = FALSE)) %>%
    distinct(group, model_type, .keep_all = TRUE) %>%
    arrange(group, model_type)
  write_tsv(calibration, file.path(output_root, "tables", "ordinal6_sameeffect_calibration.tsv"))
}

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
  file.path(output_root, "plots", "ordinal6_sameeffect_power.pdf"),
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
  file.path(output_root, "plots", "ordinal6_sameeffect_bias.pdf"),
  bias_plot,
  width = 12,
  height = 9
)

message("Combined chunk outputs written to ", output_root)
