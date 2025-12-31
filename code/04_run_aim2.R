# -----------------------------------------------------------------------------
# AIM 2 SIMULATION SCRIPT: Impact of Accuracy
#
# Desc: This script runs simulations to assess how diagnostic accuracy affects
#       statistical power, bias, MSE, and direction of effect estimates.
# -----------------------------------------------------------------------------

# --- 1. SETUP ---
message("--- Loading Packages ---")
pacman::p_load(tidyverse, broom, furrr, scales)

theme_set(theme_minimal() + theme(legend.position = "bottom"))
source("code/functions/functions.R")
source("code/common_parameters.R")
set.seed(20240423)

# Create results directories
scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}
output_root <- file.path("results", scenario_set)
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "objects"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)

# --- 2. DEFINE SIMULATION PARAMETERS ---
message("--- Defining Simulation Parameters ---")

# --- Repetitions and Parallel Cores ---
n_reps_aim2 <- as.integer(Sys.getenv("N_REPS_AIM2", unset = "1000"))
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores() - 1))
if (is.na(n_cores) || n_cores < 1) n_cores <- 1
plan(multisession, workers = n_cores)
message(paste("--- Parallel processing enabled on", n_cores, "cores ---"))

params <- get_common_params()

# --- Define OR scenarios ---
or_arrest_raw    <- params$or_arrest
or_arrest_shrunk <- params$or_arrest_shrunk_k05
or_conservative  <- params$or_conservative

scenario_definitions <- if (scenario_set == "main") {
  tribble(
    ~scenario_name,   ~or_vector,
    "ARREST_raw",     or_arrest_raw
  )
} else {
  tribble(
    ~scenario_name,   ~or_vector,
    "ARREST_shrunk",  or_arrest_shrunk
  )
}

# --- Shared Parameters (overall mortality across both arms) ---
freq_arrest        <- params$freq_arrest
p0_arrest_adjusted <- params$p0_overall
alpha_bonferroni   <- 0.05 / 5
alpha_overall      <- 0.05

# --- 3. RUN AIM 2: Impact of Accuracy (for both scenarios) ---
message("\n--- STARTING AIM 2: Impact of Accuracy ---")
start_time_aim2 <- Sys.time()

# --- Pre-calculate the true marginal OR for each scenario ---
true_or_n <- as.integer(Sys.getenv("TRUE_OR_N", unset = "10000000"))
message(paste("--- Pre-calculating true marginal ORs with large simulation (N=", format(true_or_n, big.mark = ","), ") ---", sep = ""))
true_marginal_ors <- scenario_definitions %>%
  mutate(
    true_overall_beta = map_dbl(or_vector, function(or_vec) {
      large_sim_data <- simulate_trial_data(
        or_vector = or_vec,
        freq_vector = freq_arrest,
        p0_vector = p0_arrest_adjusted,
        n = true_or_n,
        seed = 20240729
      )
      
      model <- glm(success ~ treatment, data = large_sim_data, family = binomial())
      coef(model)["treatment"]
    })
  ) %>%
  select(scenario_name, true_overall_beta)

print(true_marginal_ors)

results_aim2 <- map_dfr(1:nrow(scenario_definitions), ~{
  scenario <- scenario_definitions$scenario_name[.x]
  or_vec <- scenario_definitions$or_vector[[.x]]
  
  # Get the pre-calculated true overall beta for the current scenario
  true_beta_overall_scenario <- true_marginal_ors %>%
    filter(scenario_name == scenario) %>%
    pull(true_overall_beta)
    
  message(paste("  Running Aim 2 for scenario:", scenario))
  
  accuracy_levels <- seq(0.7, 1.0, by = 0.05)
  n_fixed <- 10000
  
  future_map_dfr(accuracy_levels, function(acc) {
    replicate_sims(
      or_vector  = or_vec,
      freq_vector = freq_arrest,
      p0_vector  = p0_arrest_adjusted,
      n          = n_fixed,
      accuracy   = acc,
      n_reps     = n_reps_aim2,
      true_overall_beta = true_beta_overall_scenario # Pass the correct value
    )
  }, .options = furrr_options(seed = TRUE), .id = "acc_id") %>%
    mutate(accuracy = accuracy_levels[as.numeric(acc_id)], scenario_name = scenario)
})

summary_aim2 <- results_aim2 %>%
  # Separate summary for overall group because it uses a different alpha
  # and has a different 'true' effect to compare against
  bind_rows(
    results_aim2 %>%
      filter(group == "Overall") %>%
      group_by(scenario_name, accuracy, group) %>%
      summarise(
        power = mean(pval < alpha_overall, na.rm = TRUE),
        bias = mean(beta - true_beta, na.rm = TRUE),
        mse = mean((beta - true_beta)^2, na.rm = TRUE),
        wrong_dir = mean(sign(or - 1) != sign(true_or - 1), na.rm = TRUE) * 100,
        .groups = "drop"
      )
  ) %>%
  filter(group != "Overall") %>% # Now remove the original overall rows
  group_by(scenario_name, accuracy, group) %>%
  summarise(
    power = mean(pval < alpha_bonferroni, na.rm = TRUE),
    bias = mean(beta - true_beta, na.rm = TRUE),
    mse = mean((beta - true_beta)^2, na.rm = TRUE),
    wrong_dir = mean(sign(or - 1) != sign(true_or - 1), na.rm = TRUE) * 100,
    .groups = "drop"
  )

plot_aim2_power <- summary_aim2 %>%
  filter(group != "Overall") %>%
  ggplot(aes(x = accuracy, y = power, color = group)) +
  geom_line() + geom_point() + facet_wrap(~ scenario_name) + 
  labs(title="Aim 2: Power vs. Accuracy", x="Accuracy", y="Power")

primary_scenario_name <- if (scenario_set == "main") "ARREST_raw" else "ARREST_shrunk"

plot_aim2_bias_arrest <- ggplot(summary_aim2 %>% filter(group != "Overall" & scenario_name == primary_scenario_name), 
                               aes(x = accuracy, y = bias, color = group)) + 
  geom_line() + geom_point() + facet_wrap(scenario_name ~ group, scales="free_y") + 
  labs(title="Aim 2: Bias vs. Accuracy (ARREST)", x="Accuracy", y="Bias") + 
  geom_hline(yintercept = 0, linetype = "dashed")

plot_aim2_bias_conservative <- NULL

plot_aim2_wrong_dir <- summary_aim2 %>%
  ggplot(aes(x = accuracy, y = wrong_dir, color = group)) +
  geom_line() + geom_point() + facet_wrap(~ scenario_name) +
  labs(title="Aim 2: Wrong Direction % vs. Accuracy", x="Accuracy", y="Wrong Direction %")

# --- 4. SAVE RESULTS ---
summary_aim2 <- summary_aim2 %>% mutate(scenario_set = scenario_set)
results_aim2 <- results_aim2 %>% mutate(scenario_set = scenario_set)

write_tsv(summary_aim2, file.path(output_root, "tables/aim2_accuracy_summary.tsv"))
write_tsv(results_aim2, file.path(output_root, "tables/aim2_raw_results.tsv.gz"))
ggsave(file.path(output_root, "plots/aim2_power_vs_accuracy.pdf"), plot_aim2_power, width = 10, height = 6)
ggsave(file.path(output_root, "plots/aim2_bias_vs_accuracy_arrest.pdf"), plot_aim2_bias_arrest, width = 10, height = 6)
if (!is.null(plot_aim2_bias_conservative)) {
  ggsave(file.path(output_root, "plots/aim2_bias_vs_accuracy_conservative.pdf"), plot_aim2_bias_conservative, width = 10, height = 6)
}
ggsave(file.path(output_root, "plots/aim2_wrongdir_vs_accuracy.pdf"), plot_aim2_wrong_dir, width = 10, height = 6)
saveRDS(plot_aim2_power, file.path(output_root, "objects/aim2_power_plot.rds"))
saveRDS(plot_aim2_bias_arrest, file.path(output_root, "objects/aim2_bias_arrest_plot.rds"))
if (!is.null(plot_aim2_bias_conservative)) {
  saveRDS(plot_aim2_bias_conservative, file.path(output_root, "objects/aim2_bias_conservative_plot.rds"))
}
saveRDS(plot_aim2_wrong_dir, file.path(output_root, "objects/aim2_wrong_dir_plot.rds"))

message(paste("--- AIM 2 COMPLETE --- (Duration:", round(difftime(Sys.time(), start_time_aim2, units = "mins"), 1), "minutes) ---"))
message("\n--- AIM 2 SIMULATION COMPLETE ---") 
