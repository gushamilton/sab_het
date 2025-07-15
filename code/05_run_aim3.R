# -----------------------------------------------------------------------------
# 05_run_aim3.R - AIM 3: Enrichment Trial (Sensitivity/Specificity Scenarios)
# -----------------------------------------------------------------------------

# --- 1. SETUP ---
message("--- Loading Packages ---")
pacman::p_load(tidyverse, broom, furrr, scales)

theme_set(theme_minimal() + theme(legend.position = "bottom"))
source("code/functions/functions.R")
set.seed(20240423)

# Create results directories
dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)
dir.create("results/objects", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)

# --- 2. DEFINE SIMULATION PARAMETERS ---
message("--- Defining Simulation Parameters ---")

# --- Repetitions and Parallel Cores ---
n_reps_aim3 <- 3000  # increased from 300 to reduce Monte Carlo error in Aim 3 power estimates
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores() - 1))
if (is.na(n_cores) || n_cores < 1) n_cores <- 1
plan(multisession, workers = n_cores)
message(paste("--- Parallel processing enabled on", n_cores, "cores ---"))

# --- Define the two main OR scenarios ---
or_arrest <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
or_conservative <- c(A = 1.0, B = 2.0, C = 0.7, D = 1.2, E = 0.8)

scenario_definitions <- tribble(
  ~scenario_name, ~or_vector,
  "ARREST",        or_arrest,
  "Conservative",  or_conservative
)

# --- Shared Parameters ---
freq_arrest        <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)
p0_arrest_raw      <- c(A = 7/60, B = 0/52, C = 11/138, D = 11/69, E = 1/69)
p0_B_corrected     <- 0.5 / (52 + 0.5)
p0_arrest_adjusted <- p0_arrest_raw
p0_arrest_adjusted["B"] <- p0_B_corrected
alpha_bonferroni   <- 0.05 / 5
alpha_overall      <- 0.05

# --- 5. RUN AIM 3: Enrichment Trial (Sensitivity/Specificity Scenarios) ---
message("\n--- STARTING AIM 3: Enrichment Trial (Sens/Spec Scenarios) ---")
start_time_aim3 <- Sys.time()

sens_spec_scenarios <- tribble(
  ~test_type,                  ~sensitivity, ~specificity,
  "Perfect (100%)",            1.00,         1.00,  # Hypothetical perfect classifier
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

reset_global_counter(nrow(all_aim3_scenarios))

aim3_results <- pmap_dfr(list(
    scenario_name = all_aim3_scenarios$scenario_name,
    or_vector = all_aim3_scenarios$or_vector,
    target_group = all_aim3_scenarios$target_group,
    sensitivity = all_aim3_scenarios$sensitivity,
    specificity = all_aim3_scenarios$specificity,
    seed = all_aim3_scenarios$seed
  ), function(scenario_name, or_vector, target_group, sensitivity, specificity, seed) {
    
  set_current_scenario(scenario_name)
  increment_global_counter()
  
  find_nns_for_scenario_sens_spec(
    target_group = target_group, sensitivity = sensitivity, specificity = specificity, seed = seed,
    or_vector = or_vector, p0_vector = p0_arrest_adjusted, freq_vector = freq_arrest,
    n_reps_per_calc = n_reps_aim3
  ) %>%
    dplyr::mutate(scenario_name = scenario_name) # <--- add scenario name to results
})

aim3_results_final <- all_aim3_scenarios %>%
  dplyr::select(scenario_name, test_type, sensitivity, specificity, target_group) %>%
  dplyr::left_join(aim3_results, by = c("scenario_name", "sensitivity", "specificity", "target_group")) # <--- join on scenario_name too

plot_aim3 <- ggplot(aim3_results_final, aes(x = test_type, y = nns_needed, color = target_group)) +
    geom_point(size=3) +
    facet_grid(scenario_name ~ target_group, scales = "free_y") +
    scale_y_log10(labels = scales::comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title="Aim 3: NNS by Scenario, Test Type, and Target Group", x="Test Type", y="Number Needed to Screen (log scale)")

# Ensure bias and mse columns are included in the final output
write_tsv(aim3_results_final, "results/tables/aim3_sens_spec_summary.tsv")
ggsave("results/plots/aim3_nns_summary.pdf", plot_aim3, width = 12, height = 8)
saveRDS(plot_aim3, "results/objects/aim3_plot.rds")

message(paste("--- AIM 3 COMPLETE --- (Duration:", round(difftime(Sys.time(), start_time_aim3, units = "mins"), 1), "minutes) ---"))
message("\n\n--- AIM 3 SIMULATION COMPLETE ---") 