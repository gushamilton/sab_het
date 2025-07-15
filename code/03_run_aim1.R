# -----------------------------------------------------------------------------
# 03_run_aim1.R - AIM 1: Power vs. Sample Size
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