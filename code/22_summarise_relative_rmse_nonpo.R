## =====================================================================
##  22_summarise_relative_rmse_nonpo.R
##  Compute relative RMSE from an existing ordinal_nonPO comparison run.
## =====================================================================

library(tidyverse)

run_tag <- Sys.getenv("RUN_TAG", unset = "")
if (!nzchar(run_tag)) {
  stop("Set RUN_TAG to an existing ordinal_nonPO comparison run.")
}

input_path <- file.path(
  "results", "supp", "ordinal_nonPO_comparison", "data",
  paste0("ordinal_nonPO_comparison_raw_", run_tag, ".tsv.gz")
)
output_root <- file.path("results", "supp", "ordinal_nonPO_comparison")
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)

raw_results <- read_tsv(input_path, show_col_types = FALSE) %>%
  mutate(
    raw_error = log_or_hat - log_or_true,
    rel_error = raw_error / abs(log_or_true)
  )

summary_tbl <- raw_results %>%
  filter(group != "A") %>%
  group_by(sample_size, accuracy, dgm, group, model_type) %>%
  summarise(
    rmse_log = sqrt(mean(raw_error^2, na.rm = TRUE)),
    relative_rmse = sqrt(mean(rel_error^2, na.rm = TRUE)),
    median_abs_rel_error = median(abs(rel_error), na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  summary_tbl,
  file.path(output_root, "tables", paste0("ordinal_nonPO_comparison_relative_rmse_", run_tag, ".tsv"))
)

plot_data <- summary_tbl %>%
  mutate(
    panel = factor(
      paste0(if_else(accuracy == 0.7, "Accuracy = 70%", "Accuracy = 100%"), "\nDGM = ", dgm),
      levels = c(
        "Accuracy = 70%\nDGM = PO",
        "Accuracy = 70%\nDGM = nonPO",
        "Accuracy = 100%\nDGM = PO",
        "Accuracy = 100%\nDGM = nonPO"
      )
    )
  )

p <- ggplot(plot_data, aes(x = sample_size, y = relative_rmse, colour = model_type)) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(group ~ panel, scales = "free_y") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_colour_manual(values = c(
    "Ordinal (polr)" = "#2166ac",
    "Binary (death)" = "#d6604d"
  )) +
  labs(
    title = "Relative RMSE by subgroup, DGM, and accuracy",
    subtitle = "RMSE divided by absolute true effect",
    x = "Total trial size",
    y = "Relative RMSE",
    colour = "Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_root, "plots", paste0("ordinal_nonPO_comparison_relative_rmse_", run_tag, ".pdf")),
  p,
  width = 14,
  height = 9
)

message("Wrote relative RMSE outputs for run tag: ", run_tag)
