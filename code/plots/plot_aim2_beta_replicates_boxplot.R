# Plot Aim 2: replicate beta estimates by subphenotype and accuracy

pacman::p_load(tidyverse, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_file <- file.path(results_root, "tables", "aim2_raw_results.tsv.gz")
output_file <- file.path(results_root, "plots", "aim2_beta_replicates_boxplot.pdf")

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

aim2 <- readr::read_tsv(input_file, show_col_types = FALSE) %>%
  filter(group != "Overall", accuracy %in% c(0.7, 1)) %>%
  mutate(
    group = factor(group, levels = c("A", "B", "C", "D", "E")),
    accuracy = factor(accuracy, levels = c(0.7, 1))
  )

true_beta_df <- aim2 %>%
  distinct(scenario_name, group, accuracy, true_beta)

plot <- ggplot(aim2, aes(x = group, y = beta, color = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.15, size = 0.6) +
  geom_segment(
    data = true_beta_df,
    aes(
      x = as.numeric(group) - 0.3,
      xend = as.numeric(group) + 0.3,
      y = true_beta,
      yend = true_beta
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.7
  ) +
  facet_grid(scenario_name ~ accuracy) +
  labs(
    title = "Aim 2: Replicate Beta Estimates by Subphenotype and Accuracy",
    x = "Subphenotype",
    y = "Estimated Beta (log OR)",
    color = "Subphenotype"
  ) +
  theme(legend.position = "none")

ggsave(output_file, plot, width = 12, height = 7)
message("Saved: ", output_file)
