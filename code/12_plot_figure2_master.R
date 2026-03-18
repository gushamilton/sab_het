## =====================================================================
##  12_plot_figure2_master.R -- Core Figure 2 (enrichment master)
## =====================================================================

pacman::p_load(tidyverse, patchwork, ggdist, scales)

theme_set(ggdist::theme_ggdist())

input_path <- "results/supp/tables/supp1_enrichment_nns_nnr.tsv"
output_root <- "results/core_figures"
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_path)) {
  stop("Missing input: ", input_path)
}

palette_groups <- c(
  B = "#A3A500",
  C = "#00BF7D",
  D = "#00B0F6",
  E = "#E76BF3"
)

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
      ),
      labels = c(
        "Perfect (100%)",
        "Near-perfect",
        "High/High",
        "Balanced",
        "High Sens/Low Spec",
        "Low Sens/High Spec"
      )
    )
  )

base_theme <- theme(
  axis.text.x = element_text(angle = 40, hjust = 1),
  panel.grid.minor = element_blank()
)

panel_a <- df %>%
  ggplot(aes(x = test_type, y = nns, color = target_group, group = target_group)) +
  geom_hline(yintercept = c(1e3, 1e4, 1e5), linetype = "dashed", color = "grey70") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_manual(values = palette_groups) +
  scale_y_log10(labels = scales::comma) +
  labs(
    x = NULL,
    y = "NNS (log scale)",
    color = "Target Group"
  ) +
  base_theme

panel_b <- df %>%
  ggplot(aes(x = test_type, y = n_required, color = target_group, group = target_group)) +
  geom_hline(yintercept = c(1e3, 1e4, 1e5), linetype = "dashed", color = "grey70") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_manual(values = palette_groups) +
  scale_y_log10(labels = scales::comma) +
  labs(
    x = NULL,
    y = "NNR (log scale)",
    color = "Target Group"
  ) +
  base_theme

# Panel C: expected bias by test type and subphenotype.
# This is deterministic from the enrichment mixing model (N-independent).
panel_c <- df %>%
  ggplot(aes(x = test_type, y = bias_log_or, color = target_group, group = target_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#F8766D") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_manual(values = palette_groups) +
  labs(
    x = NULL,
    y = "Bias (log OR)",
    color = "Target Group"
  ) +
  base_theme

figure_2 <- (panel_a / panel_b / panel_c) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(output_root, "plots", "figure2_enrichment_master.pdf"),
  figure_2,
  width = 10,
  height = 13
)

# Keep a legacy alias for the old filename request while preserving old outputs.
ggsave(
  file.path(output_root, "plots", "figure2_02_nns_nnr_patchwork.pdf"),
  figure_2,
  width = 10,
  height = 13
)
