# Plot Aim 3 replicate betas by subgroup, faceted by test type

pacman::p_load(tidyverse, scales, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
tables_root <- file.path(results_root, "tables")
input_files <- list.files(tables_root, pattern = "^aim3_detailed_.*\\.tsv$", full.names = TRUE)
output_file <- file.path(results_root, "plots", "aim3_replicate_betas_by_testtype.pdf")

if (length(input_files) == 0) {
  stop("No aim3_detailed_*.tsv files found in ", tables_root)
}

raw <- purrr::map_dfr(input_files, readr::read_tsv, show_col_types = FALSE)

test_type_map <- readr::read_tsv(file.path(results_root, "tables", "aim3_sens_spec_summary.tsv"), show_col_types = FALSE) %>%
  distinct(scenario_name, target_group, sensitivity, specificity, test_type)

plot_data <- raw %>%
  left_join(test_type_map, by = c("scenario_name", "target_group", "sensitivity", "specificity")) %>%
  mutate(
    target_group = factor(target_group, levels = c("B", "C", "D", "E")),
    test_type = factor(test_type),
    true_beta = as.numeric(true_beta),
    empirical_beta = as.numeric(empirical_beta)
  )

true_lines <- plot_data %>%
  distinct(target_group, test_type, true_beta)

plot <- ggplot(plot_data, aes(x = target_group, y = empirical_beta, color = target_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.2, size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_segment(
    data = true_lines,
    aes(
      x = as.numeric(target_group) - 0.3,
      xend = as.numeric(target_group) + 0.3,
      y = true_beta,
      yend = true_beta
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.7
  ) +
  facet_wrap(~ test_type, ncol = 3) +
  labs(
    x = "Subphenotype",
    y = "Estimated Beta (log OR)",
    color = "Subphenotype"
  )

ggsave(output_file, plot, width = 12, height = 8)
message("Saved: ", output_file)
