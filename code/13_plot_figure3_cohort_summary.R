## =====================================================================
##  13_plot_figure3_cohort_summary.R -- Core Figure 3 cohort summary
## =====================================================================

library(tidyverse)
library(patchwork)

input_prev_mort <- "results/supp/tables/supp2_cohort_prevalence_mortality.tsv"
input_n_required <- "results/supp/tables/supp2_required_n_by_cohort.tsv"
output_root <- "results/core_figures"
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_prev_mort)) stop("Missing input: ", input_prev_mort)
if (!file.exists(input_n_required)) stop("Missing input: ", input_n_required)

cohort_df <- readr::read_tsv(input_prev_mort, show_col_types = FALSE)
nreq_df <- readr::read_tsv(input_n_required, show_col_types = FALSE)

cohort_order <- c("ARREST", "Edinburgh", "Edinburgh 2", "IDISA", "SABG-PCS")

cohort_sizes <- cohort_df %>%
  distinct(cohort, n_cohort) %>%
  mutate(
    cohort = factor(cohort, levels = cohort_order),
    cohort_label = paste0(cohort, "\n(n=", n_cohort, ")")
  ) %>%
  arrange(cohort)

label_map <- setNames(as.character(cohort_sizes$cohort_label), as.character(cohort_sizes$cohort))

cohort_df <- cohort_df %>%
  mutate(
    subphenotype = factor(subphenotype, levels = c("A", "B", "C", "D", "E")),
    cohort = factor(cohort, levels = cohort_order),
    cohort_label = factor(label_map[as.character(cohort)], levels = cohort_sizes$cohort_label)
  )

nreq_df <- nreq_df %>%
  mutate(
    subphenotype = factor(subphenotype, levels = c("A", "B", "C", "D", "E")),
    cohort = factor(cohort, levels = rev(cohort_order))
  )

theme_fig <- theme_gray(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

panel_a <- cohort_df %>%
  ggplot(aes(x = subphenotype, y = prevalence)) +
  geom_boxplot(fill = "grey80", color = "grey40", width = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = cohort_label), width = 0.12, height = 0, size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Subphenotype",
    y = "Prevalence",
    color = "Cohort"
  ) +
  theme_fig +
  theme(legend.position = "bottom")

panel_b <- cohort_df %>%
  ggplot(aes(x = subphenotype, y = mortality)) +
  geom_boxplot(fill = "grey80", color = "grey40", width = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = cohort_label), width = 0.12, height = 0, size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Subphenotype",
    y = "Mortality",
    color = "Cohort"
  ) +
  theme_fig +
  theme(legend.position = "bottom")

panel_c_data <- nreq_df %>%
  filter(subphenotype %in% c("B", "C", "D", "E")) %>%
  mutate(
    n_required_plot = n_required,
    n_required_label = scales::comma(round(n_required_plot))
  )

panel_c <- panel_c_data %>%
  ggplot(aes(x = subphenotype, y = cohort, fill = n_required_plot)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = n_required_label), size = 3) +
  scale_fill_gradientn(
    colours = c("#3b0f70", "#8c1d99", "#de4968", "#fe9f6d", "#fcfdbf"),
    trans = "log10",
    breaks = c(50, 100, 300, 1000, 10000),
    labels = scales::comma
  ) +
  labs(
    x = "Subphenotype",
    y = "Cohort",
    fill = "Required N in\nTarget Subgroup\n(80% power)"
  ) +
  theme_fig +
  theme(legend.position = "right")

# Panel D: convert subgroup-specific required N to total trial sample size
# needed in an all-comers cohort, using cohort-specific subgroup prevalence.
panel_d_data <- panel_c_data %>%
  left_join(
    cohort_df %>% select(cohort, subphenotype, prevalence),
    by = c("cohort", "subphenotype")
  ) %>%
  mutate(
    n_total_required = n_required_plot / prevalence,
    n_total_label = scales::comma(round(n_total_required))
  )

panel_d <- panel_d_data %>%
  ggplot(aes(x = subphenotype, y = cohort, fill = n_total_required)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = n_total_label), size = 3) +
  scale_fill_gradientn(
    colours = c("#3b0f70", "#8c1d99", "#de4968", "#fe9f6d", "#fcfdbf"),
    trans = "log10",
    breaks = c(500, 1000, 5000, 10000, 50000, 100000),
    labels = scales::comma
  ) +
  labs(
    x = "Subphenotype",
    y = "Cohort",
    fill = "Total Trial N\nNeeded\n(80% power)"
  ) +
  theme_fig +
  theme(legend.position = "right")

left_block <- (panel_a / panel_b) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

right_block <- panel_c / panel_d

figure_3 <- (left_block | right_block) +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1.05, 1.15))

ggsave(
  file.path(output_root, "plots", "figure3.pdf"),
  figure_3,
  width = 15,
  height = 9
)

ggsave(
  file.path(output_root, "plots", "figure3.png"),
  figure_3,
  width = 15,
  height = 9,
  dpi = 200
)
