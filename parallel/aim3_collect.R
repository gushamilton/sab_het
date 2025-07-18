# -----------------------------------------------------------------------------
# aim3_collect.R - Collect and combine parallel AIM 3 results
# 
# This script collects and combines results from the parallel execution of
# Aim 3 enrichment trial scenarios. It reads individual scenario results
# from parallel/results/ and combines them into final summary files.
#
# This script is part of the parallel execution system for Aim 3 analysis.
# -----------------------------------------------------------------------------

# --- 1. SETUP ---
message("--- Loading Packages ---")
pacman::p_load(tidyverse, broom, furrr, scales)

# --- 2. DEFINE SCENARIO PARAMETERS (same as parallel script) ---
# --- Define the two main OR scenarios ---
or_arrest <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
or_conservative <- c(A = 1.0, B = 2.0, C = 0.7, D = 1.2, E = 0.8)

scenario_definitions <- tribble(
  ~scenario_name, ~or_vector,
  "ARREST",        or_arrest,
  "Conservative",  or_conservative
)

# --- Shared Parameters ---
# Updated to use overall prevalence (both treatment and control arms combined)
# From Swets et al. CID Supplementary Table 2:
# Group A: 60+55=115, Group B: 52+55=107, Group C: 138+138=276, Group D: 69+52=121, Group E: 69+70=139
# Total: 115+107+276+121+139 = 758
freq_arrest        <- c(A = 115/758, B = 107/758, C = 276/758, D = 121/758, E = 139/758)
# Updated mortality data based on paper (overall mortality across both arms)
# From Swets et al. CID Supplementary Table 2:
# Group A: (13+12)/(60+55)=25/115=21.7%, Group B: (0+8)/(52+55)=8/107=7.5%, 
# Group C: (29+24)/(138+138)=53/276=19.2%, Group D: (11+11)/(69+52)=22/121=18.2%, 
# Group E: (3+1)/(69+70)=4/139=2.9%
p0_arrest_raw      <- c(A = 25/115, B = 8/107, C = 53/276, D = 22/121, E = 4/139)
p0_B_corrected     <- 0.5 / (107 + 0.5)
p0_arrest_adjusted <- p0_arrest_raw
p0_arrest_adjusted["B"] <- p0_B_corrected

# --- 4. DEFINE ALL SCENARIOS ---
sens_spec_scenarios <- tribble(
  ~test_type,                  ~sensitivity, ~specificity,
  "Perfect (100%)",            1.00,         1.00,
  "Near-Perfect",              0.99,         0.99,
  "High Sens/High Spec",       0.95,         0.95,
  "High Sens/Low Spec",        0.95,         0.70,
  "Low Sens/High Spec",        0.70,         0.95,
  "Balanced/Moderate",         0.80,         0.80
)
target_groups_aim3 <- c("B", "C", "D", "E")

all_aim3_scenarios <- expand_grid(
  scenario_definitions,
  sens_spec_scenarios,
  target_group = target_groups_aim3
) %>% mutate(seed = 1:n() + 4000)

total_scenarios <- nrow(all_aim3_scenarios)
message(paste("--- Collecting results for", total_scenarios, "scenarios ---"))

# --- 3. COLLECT SUMMARY RESULTS ---
message("--- Collecting summary results ---")

summary_files <- list.files("parallel/results", pattern = "aim3_summary_scenario_.*\\.tsv", full.names = TRUE)
message(paste("Found", length(summary_files), "summary files"))

if (length(summary_files) == 0) {
  stop("No summary files found in results/ directory")
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

detailed_files <- list.files("parallel/results", pattern = "aim3_detailed_scenario_.*\\.tsv", full.names = TRUE)
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
dir.create("parallel/results", showWarnings = FALSE)
dir.create("parallel/results/tables", showWarnings = FALSE)

# Save summary results
write_tsv(summary_results_final, "parallel/results/tables/aim3_sens_spec_summary.tsv")
message("Summary results saved to: parallel/results/tables/aim3_sens_spec_summary.tsv")

# Save detailed results if available
if (!is.null(detailed_results_final)) {
  write_tsv(detailed_results_final, "parallel/results/tables/aim3_detailed_all_scenarios.tsv")
  message("Detailed results saved to: parallel/results/tables/aim3_detailed_all_scenarios.tsv")
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
  
  dir.create("parallel/results/plots", showWarnings = FALSE)
  dir.create("parallel/results/objects", showWarnings = FALSE)
  
  ggsave("parallel/results/plots/aim3_nns_summary.pdf", plot_aim3, width = 12, height = 8)
  saveRDS(plot_aim3, "parallel/results/objects/aim3_plot.rds")
  message("Plot saved to: parallel/results/plots/aim3_nns_summary.pdf")
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
