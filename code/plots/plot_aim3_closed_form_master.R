# Aim 3 closed-form master plot: NNS and NNR stacked

pacman::p_load(tidyverse, scales, patchwork, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_file <- file.path(results_root, "tables", "aim3_closed_form_summary.tsv")
output_file <- file.path(results_root, "plots", "aim3_closed_form_master.pdf")

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

closed <- readr::read_tsv(input_file, show_col_types = FALSE) %>%
  filter(target_group %in% c("B", "C", "D", "E"))

nns_plot <- closed %>%
  ggplot(aes(x = test_type, y = nns_closed_form, color = target_group, group = target_group)) +
  geom_hline(yintercept = c(1e3, 1e4, 1e5), linetype = "dashed", color = "grey70") +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Test Type", y = "NNS (log scale)", color = "Target Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

nnr_plot <- closed %>%
  ggplot(aes(x = test_type, y = nnr_closed_form, color = target_group, group = target_group)) +
  geom_hline(yintercept = c(1e3, 1e4, 1e5), linetype = "dashed", color = "grey70") +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Test Type", y = "NNR (log scale)", color = "Target Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

final_plot <- (nns_plot / nnr_plot) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(output_file, final_plot, width = 10, height = 10)
message("Saved: ", output_file)
