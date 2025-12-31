# -----------------------------------------------------------------------------
# aim3_collect.R - Collect and combine parallel AIM 3 results
# 
# This script collects and combines results from the parallel execution of
# Aim 3 enrichment trial scenarios. It reads individual scenario results
# from the main results/ directory and combines them into final summary files.
#
# This script is part of the parallel execution system for Aim 3 analysis.
# -----------------------------------------------------------------------------

# --- 1. SETUP ---
message("--- Loading Packages ---")
pacman::p_load(tidyverse, broom, furrr, scales)

# +++ NEW: load shared constants +++
source("../code/common_parameters.R")
params <- get_common_params()
attach(params, warn.conflicts = FALSE)

# Output paths
scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}
output_root <- file.path("..", "results", scenario_set)

# --- 2. DEFINE SCENARIO PARAMETERS (same as parallel script) ---
or_arrest_raw    <- or_arrest
or_arrest_shrunk <- or_arrest_shrunk_k05
or_conservative  <- or_conservative

scenario_definitions <- if (scenario_set == "main") {
  tribble(
    ~scenario_name,   ~or_vector,
    "ARREST_shrunk",  or_arrest_shrunk
  )
} else {
  tribble(
    ~scenario_name,   ~or_vector,
    "ARREST_raw",     or_arrest_raw,
    "Conservative",   or_conservative
  )
}

# --- 4. DEFINE ALL SCENARIOS ---

all_aim3_scenarios <- expand_grid(
  scenario_definitions,
  sens_spec_scenarios,
  target_group = target_groups_aim3
) %>% mutate(seed = 1:n() + 4000)

total_scenarios <- nrow(all_aim3_scenarios)
message(paste("--- Collecting results for", total_scenarios, "scenarios ---"))

# --- 3. COLLECT SUMMARY RESULTS ---
message("--- Collecting summary results ---")

summary_files <- list.files(output_root, pattern = "aim3_summary_scenario_.*\\.tsv", full.names = TRUE)
message(paste("Found", length(summary_files), "summary files"))

if (length(summary_files) == 0) {
  stop("No summary files found in ../results/ directory")
}

# Read and combine all summary results
summary_results <- map_dfr(summary_files, read_tsv, .id = "file_id") %>%
  mutate(
    scenario_id = as.integer(str_extract(file_id, "\\d+")),
    file_id = NULL
  ) %>%
  arrange(scenario_id)

# Add scenario metadata
summary_results_final <- all_aim3_scenarios %>%
  select(scenario_name, test_type, sensitivity, specificity, target_group) %>%
  mutate(scenario_id = 1:n()) %>%
  right_join(summary_results) %>%
  select(-scenario_id) %>%
  arrange(scenario_name, target_group, sensitivity, specificity)

# --- 4. COLLECT DETAILED RESULTS ---
message("--- Collecting detailed results ---")

detailed_files <- list.files(output_root, pattern = "aim3_detailed_scenario_.*\\.tsv", full.names = TRUE)
message(paste("Found", length(detailed_files), "detailed files"))

if (length(detailed_files) > 0) {
  # Read and combine all detailed results
  detailed_results <- map_dfr(detailed_files, read_tsv, .id = "file_id") %>%
    mutate(
      scenario_id = as.integer(str_extract(file_id, "\\d+")),
      file_id = NULL
    ) %>%
    arrange(scenario_id, sim_id)
  
  # Add scenario metadata
  detailed_results_final <- all_aim3_scenarios %>%
    select(scenario_name, test_type, sensitivity, specificity, target_group) %>%
    mutate(scenario_id = 1:n()) %>%
    right_join(detailed_results, by = "scenario_id") %>%
    select(-scenario_id) %>%
    arrange(scenario_name, target_group, sensitivity, specificity, sim_id)
} else {
  detailed_results_final <- NULL
  message("No detailed files found")
}

# --- 5. SAVE COMBINED RESULTS ---
message("--- Saving combined results ---")

# Create main results directory if it doesn't exist
dir.create("../results", showWarnings = FALSE)
dir.create("../results/tables", showWarnings = FALSE)

# Save summary results
write_tsv(summary_results_final, "../results/tables/aim3_sens_spec_summary.tsv")
message("Summary results saved to: ../results/tables/aim3_sens_spec_summary.tsv")

# Save detailed results if available
if (!is.null(detailed_results_final)) {
  write_tsv(detailed_results_final, "../results/tables/aim3_detailed_all_scenarios.tsv")
  message("Detailed results saved to: ../results/tables/aim3_detailed_all_scenarios.tsv")
}

# --- 6. CREATE PLOT ---
message("--- Creating plot ---")

if (nrow(summary_results_final) > 0) {
  plot_aim3 <- ggplot(summary_results_final, aes(x = test_type, y = nns_needed, color = target_group)) +
    geom_point(size=3) +
    facet_grid(scenario_name ~ target_group, scales = "free_y") +
    scale_y_log10(labels = scales::comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title="Aim 3: NNS by Scenario, Test Type, and Target Group", x="Test Type", y="Number Needed to Screen (log scale)")
  
  dir.create("../results/plots", showWarnings = FALSE)
  dir.create("../results/objects", showWarnings = FALSE)
  
  ggsave("../results/plots/aim3_nns_summary.pdf", plot_aim3, width = 12, height = 8)
  saveRDS(plot_aim3, "../results/objects/aim3_plot.rds")
  message("Plot saved to: ../results/plots/aim3_nns_summary.pdf")
}

# --- 7. SUMMARY STATISTICS ---
message("--- Summary Statistics ---")
message(paste("Total scenarios processed:", nrow(summary_results_final)))
message(paste("Scenarios with valid NNS:", sum(!is.na(summary_results_final$nns_needed))))
message(paste("Scenarios with missing NNS:", sum(is.na(summary_results_final$nns_needed))))

if (!is.null(detailed_results_final)) {
  message(paste("Total individual simulations:", nrow(detailed_results_final)))
  message(paste("Average bias across all scenarios:", round(mean(detailed_results_final$bias, na.rm = TRUE), 4)))
  message(paste("Average MSE across all scenarios:", round(mean(detailed_results_final$mse, na.rm = TRUE), 4)))
  message(paste("Proportion with wrong direction:", round(mean(detailed_results_final$wrong_direction, na.rm = TRUE), 4)))
}

message("--- Collection complete ---") 
