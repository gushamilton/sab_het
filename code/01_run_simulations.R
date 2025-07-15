# -----------------------------------------------------------------------------
# MASTER SIMULATION SCRIPT
#
# Desc: This script runs the full set of simulations for the HTE paper.
#       It is designed to be run as a standalone script from the terminal.
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
n_reps_aim1_2 <- 1000
n_reps_aim3 <- 2000  # increased from 300 to reduce Monte Carlo error in Aim 3 power estimates
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
# --- 3. RUN AIM 1: Power vs. Sample Size (for both scenarios) ---
message("\n--- STARTING AIM 1: Power vs. Sample Size ---")
start_time_aim1 <- Sys.time()

results_aim1 <- map_dfr(1:nrow(scenario_definitions), ~{
  scenario <- scenario_definitions$scenario_name[.x]
  or_vec <- scenario_definitions$or_vector[[.x]]
  message(paste("  Running Aim 1 for scenario:", scenario))
  
  sample_sizes <- c(500, 1000, 2000, 3000, 5000, 7500, 10000, 15000, 20000)
  
  future_map_dfr(sample_sizes, function(ss) {
    replicate_sims(
      or_vector  = or_vec,
      freq_vector = freq_arrest,
      p0_vector  = p0_arrest_adjusted,
      n          = ss,
      accuracy   = 1.0,
      n_reps     = n_reps_aim1_2
    )
  }, .options = furrr_options(seed = TRUE), .id = "size_id") %>%
    mutate(n_total = sample_sizes[as.numeric(size_id)], scenario_name = scenario)
})

power_aim1 <- results_aim1 %>%
  filter(group != "Overall") %>%
  group_by(scenario_name, n_total, group) %>%
  summarise(power = mean(pval < alpha_bonferroni, na.rm = TRUE), .groups = "drop")

plot_aim1 <- ggplot(power_aim1, aes(x = n_total, y = power, color = group)) +
  geom_line() + geom_point() +
  facet_wrap(~ scenario_name) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_x_continuous(breaks = scales::pretty_breaks(6), labels = scales::comma) +
  labs(title="Aim 1: Power vs. Sample Size", x="Total Trial Sample Size", y="Power")

write_tsv(power_aim1, "results/tables/aim1_power_summary.tsv")
write_tsv(results_aim1, "results/tables/aim1_raw_results.tsv.gz")
ggsave("results/plots/aim1_power_vs_samplesize.pdf", plot_aim1, width = 10, height = 6)
saveRDS(plot_aim1, "results/objects/aim1_power_plot.rds")

message(paste("--- AIM 1 COMPLETE --- (Duration:", round(difftime(Sys.time(), start_time_aim1, units = "mins"), 1), "minutes) ---"))

# --- 4. RUN AIM 2: Impact of Accuracy (for both scenarios) ---
message("\n--- STARTING AIM 2: Impact of Accuracy ---")
start_time_aim2 <- Sys.time()

results_aim2 <- map_dfr(1:nrow(scenario_definitions), ~{
  scenario <- scenario_definitions$scenario_name[.x]
  or_vec <- scenario_definitions$or_vector[[.x]]
  message(paste("  Running Aim 2 for scenario:", scenario))
  
  accuracy_levels <- c(1.0, 0.99, 0.95, 0.9, 0.8, 0.7)
  n_fixed <- 20000
  
  future_map_dfr(accuracy_levels, function(acc) {
    replicate_sims(
      or_vector  = or_vec,
      freq_vector = freq_arrest,
      p0_vector  = p0_arrest_adjusted,
      n          = n_fixed,
      accuracy   = acc,
      n_reps     = n_reps_aim1_2
    )
  }, .options = furrr_options(seed = TRUE), .id = "acc_id") %>%
    mutate(accuracy = accuracy_levels[as.numeric(acc_id)], scenario_name = scenario)
})

summary_aim2 <- results_aim2 %>%
  filter(group != "Overall") %>%
  group_by(scenario_name, accuracy, group) %>%
  summarise(
    power = mean(pval < alpha_bonferroni, na.rm = TRUE),
    bias = mean(log(or) - log(true_or), na.rm = TRUE),
    mse = mean((log(or) - log(true_or))^2, na.rm = TRUE),
    wrong_dir = mean(sign(or - 1) != sign(true_or - 1), na.rm = TRUE) * 100,
    .groups = "drop"
  )

plot_aim2_power <- summary_aim2 %>%
  filter(group != "Overall") %>%
  ggplot(aes(x = accuracy, y = power, color = group)) +
  geom_line() + geom_point() + facet_wrap(~ scenario_name) + labs(title="Aim 2: Power vs. Accuracy")

plot_aim2_bias_a <- ggplot(summary_aim2 %>% filter(group != "Overall" & scenario_name == "ARREST"), aes(x = accuracy, y = bias, color = group)) + geom_line() + geom_point() + facet_wrap(scenario_name ~ group, scales="free_y") + labs(title="Aim 2: Bias vs. Accuracy") + 
  geom_hline(yintercept = 0, linetype = "dashed")

plot_aim2_bias_c <- ggplot(summary_aim2 %>% filter(group != "Overall" & scenario_name == "Conservative"), aes(x = accuracy, y = bias, color = group)) + geom_line() + geom_point() + facet_wrap(scenario_name ~ group, scales="free_y") + labs(title="Aim 2: Bias vs. Accuracy") + 
  geom_hline(yintercept = 0, linetype = "dashed")

plot_aim2_wrong_dir <- summary_aim2 %>%
  ggplot(aes(x = accuracy, y = wrong_dir, color = group)) +
  geom_line() + geom_point() + facet_wrap(~ scenario_name) +
  labs(title="Aim 2: Wrong Direction % vs. Accuracy", x="Accuracy", y="Wrong Direction %")

write_tsv(summary_aim2, "results/tables/aim2_summary.tsv")
write_tsv(results_aim2, "results/tables/aim2_raw_results.tsv.gz")
ggsave("results/plots/aim2_power_vs_accuracy.pdf", plot_aim2_power, width = 10, height = 6)
ggsave("results/plots/aim2_bias_vs_accuracy_arrest.pdf", plot_aim2_bias_a, width = 10, height = 6)
ggsave("results/plots/aim2_bias_vs_accuracy_conservative.pdf", plot_aim2_bias_c, width = 10, height = 6)
ggsave("results/plots/aim2_wrong_dir_vs_accuracy.pdf", plot_aim2_wrong_dir, width = 10, height = 6)
saveRDS(plot_aim2_power, "results/objects/aim2_power_plot.rds")
saveRDS(plot_aim2_bias_a, "results/objects/aim2_bias_arrest_plot.rds")
saveRDS(plot_aim2_bias_c, "results/objects/aim2_bias_conservative_plot.rds")
saveRDS(plot_aim2_wrong_dir, "results/objects/aim2_wrong_dir_plot.rds")

message(paste("--- AIM 2 COMPLETE --- (Duration:", round(difftime(Sys.time(), start_time_aim2, units = "mins"), 1), "minutes) ---"))


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
      dplyr::mutate(scenario_name = scenario_name)
  })

  aim3_results_final <- all_aim3_scenarios %>%
    dplyr::select(scenario_name, test_type, sensitivity, specificity, target_group) %>%
    dplyr::left_join(aim3_results, by = c("scenario_name", "sensitivity", "specificity", "target_group"))

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
  message("\n\n--- ALL SIMULATIONS COMPLETE ---") 