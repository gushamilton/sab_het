# -----------------------------------------------------------------------------
# 05_run_aim3.R - AIM 3: Enrichment Trial (Sensitivity/Specificity Scenarios)
# 
# This script runs enrichment trial simulations with various sensitivity/specificity
# scenarios and saves both summary results and detailed individual simulation results.
# 
# Summary results include: NNS, bias, MSE, RMSE, sign flipping rates, etc.
# Detailed results include: individual simulation betas, p-values, bias, MSE for each rep.
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
n_reps_aim3 <- 200  # increased from 300 to reduce Monte Carlo error in Aim 3 power estimates
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

# Function to run one scenario and save detailed results
run_scenario_with_details <- function(scenario_name, or_vector, target_group, sensitivity, specificity, seed) {
  set_current_scenario(scenario_name)
  increment_global_counter()
  
  # Run the main function to get summary results
  summary_result <- find_nns_for_scenario_sens_spec(
    target_group = target_group, sensitivity = sensitivity, specificity = specificity, seed = seed,
    or_vector = or_vector, p0_vector = p0_arrest_adjusted, freq_vector = freq_arrest,
    n_reps_per_calc = n_reps_aim3
  ) %>%
    dplyr::mutate(scenario_name = scenario_name)
  
  # Also run the final scenario to get individual simulation results
  if (!is.na(summary_result$nns_needed)) {
    final_nns <- summary_result$nns_needed
    final_scenario_result <- run_enrichment_scenario_sens_spec(
      n_screened = final_nns, target_group = target_group, sensitivity = sensitivity, specificity = specificity,
      or_vector = or_vector, p0_vector = p0_arrest_adjusted, freq_vector = freq_arrest,
      n_reps = n_reps_aim3 * 2, base_seed = seed + 999
    )
    
    # Create detailed results dataframe
    true_beta <- log(or_vector[target_group])
    detailed_results <- tibble(
      scenario_name = scenario_name,
      target_group = target_group,
      sensitivity = sensitivity,
      specificity = specificity,
      sim_id = 1:length(final_scenario_result$beta_hats),
      true_beta = true_beta,
      empirical_beta = final_scenario_result$beta_hats,
      p_value = final_scenario_result$p_values,
      n_enrolled = final_scenario_result$mean_n_enrolled,
      bias = empirical_beta - true_beta,
      mse = (empirical_beta - true_beta)^2,
      significant = p_value < 0.05,
      wrong_direction = sign(empirical_beta) != sign(true_beta),
      significant_wrong = significant & wrong_direction
    )
    
    # Save detailed results to file
    detailed_filename <- sprintf("results/tables/aim3_detailed_%s_%s_sens%.2f_spec%.2f.tsv", 
                                scenario_name, target_group, sensitivity, specificity)
    write_tsv(detailed_results, detailed_filename)
  }
  
  return(summary_result)
}

aim3_results <- pmap_dfr(list(
    scenario_name = all_aim3_scenarios$scenario_name,
    or_vector = all_aim3_scenarios$or_vector,
    target_group = all_aim3_scenarios$target_group,
    sensitivity = all_aim3_scenarios$sensitivity,
    specificity = all_aim3_scenarios$specificity,
    seed = all_aim3_scenarios$seed
  ), run_scenario_with_details)

aim3_results_final <- all_aim3_scenarios %>%
  dplyr::select(scenario_name, test_type, sensitivity, specificity, target_group) %>%
  dplyr::left_join(aim3_results, by = c("scenario_name", "sensitivity", "specificity", "target_group")) # <--- join on scenario_name too

plot_aim3 <- ggplot(aim3_results_final, aes(x = test_type, y = nns_needed, color = target_group)) +
    geom_point(size=3) +
    facet_grid(scenario_name ~ target_group, scales = "free_y") +
    scale_y_log10(labels = scales::comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title="Aim 3: NNS by Scenario, Test Type, and Target Group", x="Test Type", y="Number Needed to Screen (log scale)")

# Save summary results with enhanced metrics
write_tsv(aim3_results_final, "results/tables/aim3_sens_spec_summary.tsv")
ggsave("results/plots/aim3_nns_summary.pdf", plot_aim3, width = 12, height = 8)
saveRDS(plot_aim3, "results/objects/aim3_plot.rds")

message(paste("--- AIM 3 COMPLETE --- (Duration:", round(difftime(Sys.time(), start_time_aim3, units = "mins"), 1), "minutes) ---"))
message("\n\n--- AIM 3 SIMULATION COMPLETE ---") 