# Plot Aim 2 wrong-direction rate vs accuracy

pacman::p_load(tidyverse, scales, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "supp")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_file <- file.path(results_root, "tables", "aim2_accuracy_summary.tsv")
output_file <- file.path(results_root, "plots", "aim2_wrongdir_vs_accuracy.pdf")

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

aim2 <- readr::read_tsv(input_file, show_col_types = FALSE)

plot_title <- if (scenario_set == "supp") "Aim 2: Wrong Direction % vs Accuracy (Sensitivity)" else NULL

wrong_plot <- aim2 %>%
  filter(!group %in% c("Overall", "A")) %>%
  ggplot(aes(x = accuracy, y = wrong_dir, color = group)) +
  geom_line() +
  scale_x_continuous(breaks = sort(unique(aim2$accuracy))) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Accuracy", y = "Wrong Direction %", title = plot_title)

ggsave(output_file, wrong_plot, width = 10, height = 6)
message("Saved: ", output_file)
