## =====================================================================
##  21_plot_ordinal_po_nonpo_centered_bias_patchwork.R
##  Trial figure: centered bias panels instead of raw-effect panels.
##  Uses existing ordinal_nonPO comparison output.
## =====================================================================

library(tidyverse)
library(patchwork)

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
    model_type = factor(model_type, levels = c("Binary (death)", "Ordinal (polr)")),
    centered_bias = log_or_hat - log_or_true
  )

available_n <- sort(unique(raw_results$sample_size))
power_n <- if (2500 %in% available_n) 2500 else min(available_n)

model_colours <- c(
  "Ordinal (polr)" = "#2166ac",
  "Binary (death)" = "#d6604d"
)

make_centered_bias_plot <- function(dgm_name, panel_title) {
  plot_data <- raw_results %>%
    filter(sample_size == 20000, dgm == dgm_name, group != "A") %>%
    group_by(accuracy, group, model_type) %>%
    mutate(
      lo = quantile(centered_bias, 0.01, na.rm = TRUE),
      hi = quantile(centered_bias, 0.99, na.rm = TRUE),
      keep = centered_bias >= lo & centered_bias <= hi
    ) %>%
    ungroup() %>%
    filter(keep)

  ggplot(plot_data, aes(x = group, y = centered_bias, color = model_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.6),
      alpha = 0.2,
      size = 0.75,
      stroke = 0
    ) +
    facet_wrap(
      ~ accuracy,
      ncol = 2,
      labeller = labeller(accuracy = c(`0.7` = "Accuracy = 70%", `1` = "Accuracy = 100%"))
    ) +
    scale_color_manual(values = model_colours) +
    labs(
      title = panel_title,
      subtitle = "Centered at method-specific truth; trimmed to 1st-99th percentile",
      x = "Subphenotype",
      y = "Estimate - true effect",
      color = "Model"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

make_power_scatter <- function(dgm_name, panel_title) {
  plot_data <- raw_results %>%
    filter(sample_size == power_n, dgm == dgm_name, group != "A") %>%
    group_by(accuracy, group, model_type) %>%
    summarise(
      power = mean(
        {
          z <- qnorm(1 - params$alpha_primary / 2)
          ci_lower <- log_or_hat - z * se
          ci_upper <- log_or_hat + z * se
          significant <- !is.na(log_or_hat) & !is.na(se) & (ci_lower > 0 | ci_upper < 0)
          significant & (sign(log_or_hat) == sign(log_or_true))
        },
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = model_type, values_from = power)

  if (!("Binary (death)" %in% names(plot_data))) {
    plot_data[["Binary (death)"]] <- NA_real_
  }
  if (!("Ordinal (polr)" %in% names(plot_data))) {
    plot_data[["Ordinal (polr)"]] <- NA_real_
  }

  ggplot(plot_data, aes(x = `Binary (death)`, y = `Ordinal (polr)`, color = group, shape = factor(accuracy))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 2.8, alpha = 0.9) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = panel_title,
      subtitle = "Points above diagonal indicate ordinal power gain",
      x = "Binary power",
      y = "Ordinal power",
      color = "Subphenotype",
      shape = "Accuracy"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

p_a <- make_centered_bias_plot("PO", "A. PO: centered bias at N = 20,000")
p_b <- make_power_scatter("PO", paste0("B. PO: power comparison at N = ", scales::comma(power_n)))
p_c <- make_centered_bias_plot("nonPO", "C. nonPO: centered bias at N = 20,000")
p_d <- make_power_scatter("nonPO", paste0("D. nonPO: power comparison at N = ", scales::comma(power_n)))

combined <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(output_root, paste0("ordinal_nonPO_centered_bias_patchwork_", run_tag, ".pdf")),
  combined,
  width = 16,
  height = 12
)

message("Wrote centered-bias patchwork figure to ", output_root)
