# Plot Aim 2 bias vs accuracy (single panel, log y-axis)

pacman::p_load(tidyverse, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "supp")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_file <- file.path(results_root, "tables", "aim2_accuracy_summary.tsv")
output_file <- file.path(results_root, "plots", "aim2_bias_vs_accuracy_log.pdf")

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

aim2 <- readr::read_tsv(input_file, show_col_types = FALSE)

bias_plot <- aim2 %>%
  filter(!group %in% c("Overall", "A")) %>%
  mutate(bias_mag = abs(bias)) %>%
  ggplot(aes(x = accuracy, y = bias_mag, color = group)) +
  geom_line() +
  scale_y_log10() +
  labs(title = "Aim 2: |Bias| vs Accuracy (log y)", x = "Accuracy", y = "|Bias| (log OR, log10)", color = "Group")

ggsave(output_file, bias_plot, width = 12, height = 6)
message("Saved: ", output_file)
