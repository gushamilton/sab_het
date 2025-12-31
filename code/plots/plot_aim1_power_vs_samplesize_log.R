# Plot Aim 1 power vs sample size with log-scaled x-axis

pacman::p_load(tidyverse, scales, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "supp")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_file <- file.path(results_root, "tables", "aim1_power_summary.tsv")
output_file <- file.path(results_root, "plots", "aim1_power_vs_samplesize_log.pdf")

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

aim1 <- readr::read_tsv(input_file, show_col_types = FALSE)

breaks <- c(100, 200, 500, 1000, 2000, 5000, 10000, 20000)
breaks <- breaks[breaks <= max(aim1$n_total)]

plot_title <- if (scenario_set == "supp") {
  "Aim 1: Power vs Sample Size (log scale) - Sensitivity"
} else {
  "Aim 1: Power vs Sample Size (log scale)"
}

gg <- ggplot(aim1, aes(x = n_total, y = power, color = group)) +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_x_log10(limits = c(100, NA), breaks = breaks, labels = scales::comma) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = plot_title, x = "Total Trial Sample Size (log10)", y = "Power")

ggsave(output_file, gg, width = 10, height = 6)
message("Saved: ", output_file)
