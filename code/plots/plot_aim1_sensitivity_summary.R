# Plot Aim 1 sensitivity summary: power vs sample size by cohort

pacman::p_load(tidyverse, scales, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_file <- file.path(results_root, "tables", "aim1_sensitivity_summary.tsv")
output_file <- file.path(results_root, "plots", "aim1_sensitivity_power_by_cohort.pdf")

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

aim1 <- readr::read_tsv(input_file, show_col_types = FALSE)

plot <- aim1 %>%
  filter(group != "A") %>%
  ggplot(aes(x = n_total, y = power, color = group)) +
  geom_line() +
  facet_wrap(~ cohort) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_x_log10(limits = c(100, NA), breaks = scales::pretty_breaks(6), labels = scales::comma) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Aim 1 Sensitivity: Power vs Sample Size by Cohort",
    x = "Total Trial Sample Size (log10)",
    y = "Power",
    color = "Group"
  )

ggsave(output_file, plot, width = 11, height = 7)
message("Saved: ", output_file)
