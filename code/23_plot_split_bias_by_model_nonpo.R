## =====================================================================
##  23_plot_split_bias_by_model_nonpo.R
##  View-only diagnostic plots:
##   - separate PO and nonPO replicate plots
##   - facets by model type and accuracy
##   - free y-axis
##   - black target bars for method-specific truth
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table

run_tag <- Sys.getenv("RUN_TAG", unset = "")
if (!nzchar(run_tag)) {
  stop("Set RUN_TAG to an existing ordinal_nonPO comparison run.")
}

input_path <- file.path(
  "results", "supp", "ordinal_nonPO_comparison", "data",
  paste0("ordinal_nonPO_comparison_raw_", run_tag, ".tsv.gz")
)
output_root <- file.path("results", "supp", "ordinal_nonPO_comparison", "plots")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)

raw_results <- read_tsv(input_path, show_col_types = FALSE) %>%
  mutate(
    group = factor(group, levels = subphenotypes$subphenotype),
    model_type = factor(model_type, levels = c("Binary (death)", "Ordinal (polr)"))
  )

truth_df <- raw_results %>%
  distinct(sample_size, dgm, group, model_type, accuracy, log_or_true)

model_colours <- c(
  "Ordinal (polr)" = "#2166ac",
  "Binary (death)" = "#d6604d"
)

make_plot <- function(dgm_name, file_stub, title_text) {
  plot_data <- raw_results %>%
    filter(dgm == dgm_name, sample_size == 20000, group != "A") %>%
    group_by(model_type, accuracy, group) %>%
    mutate(
      lo = quantile(log_or_hat, 0.01, na.rm = TRUE),
      hi = quantile(log_or_hat, 0.99, na.rm = TRUE),
      keep = log_or_hat >= lo & log_or_hat <= hi
    ) %>%
    ungroup() %>%
    filter(keep)

  p <- ggplot(plot_data, aes(x = group, y = log_or_hat, color = model_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_jitter(width = 0.14, height = 0, alpha = 0.18, size = 0.7, stroke = 0) +
    geom_errorbar(
      data = truth_df %>% filter(dgm == dgm_name, sample_size == 20000, group != "A"),
      aes(x = group, ymin = log_or_true, ymax = log_or_true),
      inherit.aes = FALSE,
      width = 0.55,
      linewidth = 1.0,
      color = "black"
    ) +
    facet_grid(
      model_type ~ accuracy,
      scales = "free_y",
      labeller = labeller(accuracy = c(`0.7` = "Accuracy = 70%", `1` = "Accuracy = 100%"))
    ) +
    scale_color_manual(values = model_colours) +
    labs(
      title = paste0(title_text, " at N = 20,000"),
      subtitle = "Black bars show method-specific targets; trimmed to 1st-99th percentile",
      x = "Subphenotype",
      y = "Estimated log-effect",
      color = "Model"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(
    file.path(output_root, paste0(file_stub, "_", run_tag, ".pdf")),
    p,
    width = 12,
    height = 8
  )
}

make_plot("PO", "ordinal_nonPO_split_bias_PO", "PO: replicate estimates by model")
make_plot("nonPO", "ordinal_nonPO_split_bias_nonPO", "nonPO: replicate estimates by model")

message("Wrote split PO/nonPO bias plots to ", output_root)
