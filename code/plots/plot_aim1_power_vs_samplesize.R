# Plot Aim 1 power vs sample size (ARREST + Conservative in supp; shrunk in main)

pacman::p_load(tidyverse, scales, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "supp")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_file <- file.path(results_root, "tables", "aim1_power_summary.tsv")
output_file <- file.path(results_root, "plots", "aim1_power_vs_samplesize.pdf")

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

aim1 <- readr::read_tsv(input_file, show_col_types = FALSE)

gg <- ggplot(aim1, aes(x = n_total, y = power, color = group)) +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_x_continuous(limits = c(100, NA), breaks = scales::pretty_breaks(6), labels = scales::comma) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Aim 1: Power vs Sample Size", x = "Total Trial Sample Size", y = "Power")

ggsave(output_file, gg, width = 10, height = 6)
message("Saved: ", output_file)
