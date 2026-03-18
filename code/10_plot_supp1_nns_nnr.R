## =====================================================================
##  10_plot_supp1_nns_nnr.R -- Recreate enrichment NNS/NNR comparison plot
## =====================================================================

library(tidyverse)
library(patchwork)

input_path <- "results/supp/tables/supp1_enrichment_nns_nnr.tsv"
output_dir <- "results/supp/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_path)) {
  stop("Missing input: ", input_path)
}

df <- readr::read_tsv(input_path, show_col_types = FALSE) %>%
  mutate(
    target_group = factor(target_group, levels = c("B", "C", "D", "E")),
    test_type = factor(
      test_type,
      levels = c(
        "Perfect",
        "Near-perfect",
        "High/High",
        "Balanced",
        "High Sens/Low Spec",
        "Low Sens/High Spec"
      )
    )
  )

base_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p_nns <- df %>%
  ggplot(aes(x = test_type, y = nns, color = target_group, group = target_group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_y_log10(labels = scales::label_comma()) +
  labs(
    title = "A",
    x = "Test Type",
    y = "NNS (log scale)",
    color = "Target Group"
  ) +
  base_theme

p_nnr <- df %>%
  ggplot(aes(x = test_type, y = n_required, color = target_group, group = target_group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_y_log10(labels = scales::label_comma()) +
  labs(
    title = "B",
    x = "Test Type",
    y = "NNR (log scale)",
    color = "Target Group"
  ) +
  base_theme

combined <- p_nns / p_nnr + plot_layout(guides = "collect")

ggsave(
  file.path(output_dir, "figureS1_nns_nnr_comparison_refactor.pdf"),
  combined,
  width = 10,
  height = 10
)

