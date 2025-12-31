# Plot Aim 3 NNS comparison (simulation vs closed form)

pacman::p_load(tidyverse, scales, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "supp")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_sim <- file.path(results_root, "tables", "aim3_sens_spec_summary.tsv")
input_closed <- file.path(results_root, "tables", "aim3_closed_form_summary.tsv")
output_file <- file.path(results_root, "plots", "aim3_nns_comparison.pdf")

if (!file.exists(input_sim)) {
  stop("Missing input: ", input_sim)
}
if (!file.exists(input_closed)) {
  stop("Missing input: ", input_closed)
}

max_nns <- 1e6

sim <- readr::read_tsv(input_sim, show_col_types = FALSE) %>%
  transmute(
    scenario_name,
    test_type = test_type,
    target_group,
    nns = nns_needed,
    method = "Simulation",
    missing = is.na(nns_needed)
  )

closed <- readr::read_tsv(input_closed, show_col_types = FALSE) %>%
  transmute(
    scenario_name,
    test_type = test_type,
    target_group,
    nns = nns_closed_form,
    method = "Closed Form",
    missing = FALSE
  )

combined <- bind_rows(sim, closed)

nns_plot <- combined %>%
  filter(target_group %in% c("B", "C", "D", "E")) %>%
  mutate(nns_plot = ifelse(missing, max_nns, nns)) %>%
  filter(is.finite(nns_plot)) %>%
  ggplot(aes(x = test_type, y = nns_plot, color = method, group = method)) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(scenario_name ~ target_group, scales = "free_y") +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "Aim 3: NNS Comparison (Simulation vs Closed Form)",
    x = "Test Type",
    y = "Number Needed to Screen (log scale)",
    color = "Method"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(output_file, nns_plot, width = 12, height = 8)
message("Saved: ", output_file)
