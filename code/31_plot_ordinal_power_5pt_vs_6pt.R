## =====================================================================
##  31_plot_ordinal_power_5pt_vs_6pt.R -- Supplementary sensitivity figure
## =====================================================================

library(tidyverse)

run_tag_5pt <- Sys.getenv("RUN_TAG_5PT", unset = "paper_final_5pt_20260324")
run_tag_6pt <- Sys.getenv("RUN_TAG_6PT", unset = "ordinal6_figure5_20260325")

metrics_5pt_path <- file.path(
  "results", "supp", "ordinal_nonPO_comparison", "tables",
  paste0("ordinal_nonPO_comparison_metrics_", run_tag_5pt, ".tsv")
)
metrics_6pt_path <- file.path(
  "results", "supp", "ordinal_nonPO_comparison", "tables",
  paste0("ordinal_nonPO_comparison_metrics_", run_tag_6pt, ".tsv")
)

out_dir <- file.path("results", "paper", "final", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

metrics_5pt <- readr::read_tsv(metrics_5pt_path, show_col_types = FALSE) %>%
  filter(model_type == "Ordinal (polr)", group != "A") %>%
  rename(power_5pt = power)

metrics_6pt <- readr::read_tsv(metrics_6pt_path, show_col_types = FALSE) %>%
  filter(model_type == "Ordinal (polr)", group != "A") %>%
  rename(power_6pt = power)

plot_df <- metrics_5pt %>%
  select(sample_size, accuracy, dgm, group, power_5pt) %>%
  inner_join(
    metrics_6pt %>% select(sample_size, accuracy, dgm, group, power_6pt),
    by = c("sample_size", "accuracy", "dgm", "group")
  ) %>%
  mutate(
    dgm = factor(dgm, levels = c("PO", "nonPO")),
    accuracy = factor(accuracy, levels = c(0.7, 1.0), labels = c("70%", "100%")),
    group = factor(group, levels = c("B", "C", "D", "E"))
  )

p <- ggplot(plot_df, aes(x = power_5pt, y = power_6pt, color = group, shape = accuracy, size = sample_size)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(alpha = 0.9) +
  facet_wrap(~ dgm, nrow = 1) +
  scale_color_manual(values = c("B" = "#1b9e77", "C" = "#d95f02", "D" = "#7570b3", "E" = "#e7298a")) +
  scale_shape_manual(values = c("70%" = 16, "100%" = 17)) +
  scale_size_continuous(
    breaks = c(500, 1000, 1500, 2000, 3000, 5000, 10000, 20000),
    range = c(1.8, 4.8),
    labels = scales::comma
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Sensitivity of ordinal power to 5-point versus 6-point outcome scales",
    subtitle = "Ordinal model only; points above the diagonal favour the 6-point scale",
    x = "Power under 5-point ordinal scale",
    y = "Power under 6-point ordinal scale",
    color = "Subgroup",
    shape = "Accuracy",
    size = "Total N"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "FigureS1_ordinal_power_5pt_vs_6pt.pdf"), p, width = 11, height = 5.5)

message("Wrote supplementary 5-point vs 6-point power figure to ", out_dir)
