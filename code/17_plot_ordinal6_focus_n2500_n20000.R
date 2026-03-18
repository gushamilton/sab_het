## =====================================================================
##  17_plot_ordinal6_focus_n2500_n20000.R
##  Focused plots:
##   - Power at N=2500
##   - Bias replicate plot at N=20000 (side-by-side Binary vs Ordinal)
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")
source("code/functions/02_metrics.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table
or_vector <- setNames(subphenotypes$or_arrest_shrunk, subphenotypes$subphenotype)
alpha <- params$alpha_primary

output_root <- "results/supp/ordinal6_sameeffect"
run_tag <- Sys.getenv("RUN_TAG", unset = "n2500_n20000")
tag_suffix <- if (nzchar(run_tag)) paste0("_", run_tag) else ""

raw_path <- file.path(output_root, "data", paste0("ordinal6_sameeffect_raw", tag_suffix, ".tsv.gz"))
if (!file.exists(raw_path)) {
  stop("Missing raw file: ", raw_path)
}

raw_results <- read_tsv(raw_path, show_col_types = FALSE) %>%
  mutate(
    group = factor(group, levels = subphenotypes$subphenotype),
    model_type = factor(model_type, levels = c("Binary", "Ordinal"))
  )

metrics <- raw_results %>%
  group_by(sample_size, accuracy, group, model_type) %>%
  group_modify(~ {
    g <- as.character(.y$group)
    truth <- log(or_vector[[g]])
    is_null <- abs(truth) < 1e-6
    out <- summarise_metrics(.x, log_or_true = truth, alpha = alpha, is_null = is_null)
    tibble(
      type1 = out$type1,
      power = out$power,
      type_s = out$type_s,
      median_signed_bias_log_or = median(.x$log_or_hat - .x$log_or_true, na.rm = TRUE)
    )
  }) %>%
  ungroup()

write_tsv(metrics, file.path(output_root, "tables", paste0("ordinal6_sameeffect_focus_metrics", tag_suffix, ".tsv")))

power_2500 <- metrics %>%
  filter(sample_size == 2500, group != "A") %>%
  select(accuracy, group, model_type, power) %>%
  pivot_wider(names_from = model_type, values_from = power) %>%
  ggplot(aes(x = Binary, y = Ordinal, color = group, shape = factor(accuracy))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3, alpha = 0.95) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "N=2,500: Power comparison",
    subtitle = "Ordinal power vs Binary power",
    x = "Binary power",
    y = "Ordinal power",
    color = "Subphenotype",
    shape = "Accuracy"
  ) +
  theme_minimal()

ggsave(
  file.path(output_root, "plots", paste0("ordinal6_sameeffect_power_scatter_n2500", tag_suffix, ".pdf")),
  power_2500,
  width = 8,
  height = 7
)

truth_df <- subphenotypes %>%
  transmute(
    group = factor(subphenotype, levels = subphenotypes$subphenotype),
    log_or_true = log(or_arrest_shrunk)
  ) %>%
  filter(group != "A")

bias_20000 <- raw_results %>%
  filter(sample_size == 20000, group != "A") %>%
  group_by(accuracy, group, model_type) %>%
  mutate(
    lo = quantile(log_or_hat, 0.01, na.rm = TRUE),
    hi = quantile(log_or_hat, 0.99, na.rm = TRUE),
    keep = log_or_hat >= lo & log_or_hat <= hi
  ) %>%
  ungroup() %>%
  filter(keep)

bias_plot <- bias_20000 %>%
  ggplot(aes(x = group, y = log_or_hat, color = model_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.6),
    alpha = 0.2,
    size = 0.75,
    stroke = 0
  ) +
  geom_errorbar(
    data = truth_df,
    aes(x = group, ymin = log_or_true, ymax = log_or_true),
    inherit.aes = FALSE,
    width = 0.55,
    linewidth = 1.0,
    color = "black"
  ) +
  facet_wrap(~ accuracy, ncol = 2) +
  labs(
    title = "N=20,000: replicate log(OR) estimates (trimmed 1st-99th percentile)",
    subtitle = "Black bars: true log(OR); Binary/Ordinal shown side-by-side",
    x = "Subphenotype",
    y = "Estimated log(OR)",
    color = "Model"
  ) +
  theme_minimal()

ggsave(
  file.path(output_root, "plots", paste0("ordinal6_sameeffect_bias_replicates_n20000_trimmed", tag_suffix, ".pdf")),
  bias_plot,
  width = 11,
  height = 7
)

message("Wrote focused N=2500 power and N=20000 bias plots.")
