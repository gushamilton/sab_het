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
set.seed(20240423)

# Create results directories
dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)
dir.create("results/objects", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)

# --- 2. DEFINE SIMULATION PARAMETERS ---
message("--- Defining Simulation Parameters ---")

# --- Repetitions and Parallel Cores ---
n_reps_aim2 <- 1000
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

# --- 3. RUN AIM 2: Impact of Accuracy (for both scenarios) ---
message("\n--- STARTING AIM 2: Impact of Accuracy ---")
start_time_aim2 <- Sys.time()

# --- Pre-calculate the true marginal OR for each scenario ---
message("--- Pre-calculating true marginal ORs with large simulation (N=10M) ---")
true_marginal_ors <- scenario_definitions %>%
  mutate(
    true_overall_beta = map_dbl(or_vector, function(or_vec) {
      large_sim_data <- simulate_trial_data(
        or_vector = or_vec,
        freq_vector = freq_arrest,
        p0_vector = p0_arrest_adjusted,
        n = 10e6, # 10 million
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
  
  accuracy_levels <- c(1.0, 0.99, 0.95, 0.9, 0.8, 0.7)
  n_fixed <- 20000
  
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

plot_aim2_bias_arrest <- ggplot(summary_aim2 %>% filter(group != "Overall" & scenario_name == "ARREST"), 
                                aes(x = accuracy, y = bias, color = group)) + 
  geom_line() + geom_point() + facet_wrap(scenario_name ~ group, scales="free_y") + 
  labs(title="Aim 2: Bias vs. Accuracy (ARREST)", x="Accuracy", y="Bias") + 
  geom_hline(yintercept = 0, linetype = "dashed")

plot_aim2_bias_conservative <- ggplot(summary_aim2 %>% filter(group != "Overall" & scenario_name == "Conservative"), 
                                     aes(x = accuracy, y = bias, color = group)) + 
  geom_line() + geom_point() + facet_wrap(scenario_name ~ group, scales="free_y") + 
  labs(title="Aim 2: Bias vs. Accuracy (Conservative)", x="Accuracy", y="Bias") + 
  geom_hline(yintercept = 0, linetype = "dashed")

plot_aim2_wrong_dir <- summary_aim2 %>%
  ggplot(aes(x = accuracy, y = wrong_dir, color = group)) +
  geom_line() + geom_point() + facet_wrap(~ scenario_name) +
  labs(title="Aim 2: Wrong Direction % vs. Accuracy", x="Accuracy", y="Wrong Direction %")

# --- 4. SAVE RESULTS ---
write_tsv(summary_aim2, "results/tables/aim2_accuracy_summary.tsv")
write_tsv(results_aim2, "results/tables/aim2_raw_results.tsv.gz")
ggsave("results/plots/aim2_power_vs_accuracy.pdf", plot_aim2_power, width = 10, height = 6)
ggsave("results/plots/aim2_bias_vs_accuracy_arrest.pdf", plot_aim2_bias_arrest, width = 10, height = 6)
ggsave("results/plots/aim2_bias_vs_accuracy_conservative.pdf", plot_aim2_bias_conservative, width = 10, height = 6)
ggsave("results/plots/aim2_wrongdir_vs_accuracy.pdf", plot_aim2_wrong_dir, width = 10, height = 6)
saveRDS(plot_aim2_power, "results/objects/aim2_power_plot.rds")
saveRDS(plot_aim2_bias_arrest, "results/objects/aim2_bias_arrest_plot.rds")
saveRDS(plot_aim2_bias_conservative, "results/objects/aim2_bias_conservative_plot.rds")
saveRDS(plot_aim2_wrong_dir, "results/objects/aim2_wrong_dir_plot.rds")

message(paste("--- AIM 2 COMPLETE --- (Duration:", round(difftime(Sys.time(), start_time_aim2, units = "mins"), 1), "minutes) ---"))
message("\n--- AIM 2 SIMULATION COMPLETE ---") 