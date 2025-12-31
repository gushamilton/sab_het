# -----------------------------------------------------------------------------
# 03_run_aim1.R - AIM 1: Power vs. Sample Size
# -----------------------------------------------------------------------------

# --- 1. SETUP ---
message("--- Loading Packages ---")
pacman::p_load(tidyverse, broom, furrr, scales)

theme_set(theme_minimal() + theme(legend.position = "bottom"))
source("code/functions/functions.R")
source("code/functions/aim1_closed_form.R")
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
n_reps_aim1_2 <- as.integer(Sys.getenv("N_REPS_AIM1_2", unset = "1000"))
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

# --- 3. RUN AIM 1: Power vs. Sample Size (for both scenarios) ---
message("\n--- STARTING AIM 1: Power vs. Sample Size ---")
start_time_aim1 <- Sys.time()

sample_sizes <- seq(100, 20000, by = 100)
aim1_mode <- Sys.getenv("AIM1_MODE", unset = "closed_form")

if (aim1_mode == "closed_form") {
  power_aim1 <- expand_grid(
    scenario_definitions,
    n_total = sample_sizes,
    group = names(params$or_arrest)
  ) %>%
    rowwise() %>%
    mutate(
      power = calc_power_closed_form(
        n_total = n_total,
        p0 = p0_arrest_adjusted[group],
        or = or_vector[group],
        prevalence = freq_arrest[group],
        alpha = alpha_bonferroni
      )
    ) %>%
    ungroup() %>%
    select(scenario_name, n_total, group, power)
  results_aim1 <- NULL
} else {
  results_aim1 <- map_dfr(1:nrow(scenario_definitions), ~{
    scenario <- scenario_definitions$scenario_name[.x]
    or_vec <- scenario_definitions$or_vector[[.x]]
    message(paste("  Running Aim 1 for scenario:", scenario))
    
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
}

plot_aim1 <- ggplot(power_aim1, aes(x = n_total, y = power, color = group)) +
  geom_line() + geom_point() +
  facet_wrap(~ scenario_name) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_x_continuous(breaks = scales::pretty_breaks(6), labels = scales::comma) +
  labs(title="Aim 1: Power vs. Sample Size", x="Total Trial Sample Size", y="Power")

power_aim1 <- power_aim1 %>% mutate(scenario_set = scenario_set)
write_tsv(power_aim1, file.path(output_root, "tables/aim1_power_summary.tsv"))

if (!is.null(results_aim1)) {
  results_aim1 <- results_aim1 %>% mutate(scenario_set = scenario_set)
  write_tsv(results_aim1, file.path(output_root, "tables/aim1_raw_results.tsv.gz"))
}
ggsave(file.path(output_root, "plots/aim1_power_vs_samplesize.pdf"), plot_aim1, width = 10, height = 6)
saveRDS(plot_aim1, file.path(output_root, "objects/aim1_power_plot.rds"))

if (aim1_mode == "closed_form" && Sys.getenv("AIM1_VALIDATE", unset = "0") == "1") {
  validate_sizes <- as.integer(strsplit(Sys.getenv("AIM1_VALIDATE_SIZES", unset = "1000,5000,10000,20000"), ",")[[1]])
  validate_sizes <- validate_sizes[!is.na(validate_sizes)]
  validate_reps <- as.integer(Sys.getenv("N_REPS_AIM1_VALIDATE", unset = "200"))
  
  sim_validate <- map_dfr(1:nrow(scenario_definitions), ~{
    scenario <- scenario_definitions$scenario_name[.x]
    or_vec <- scenario_definitions$or_vector[[.x]]
    future_map_dfr(validate_sizes, function(ss) {
      replicate_sims(
        or_vector  = or_vec,
        freq_vector = freq_arrest,
        p0_vector  = p0_arrest_adjusted,
        n          = ss,
        accuracy   = 1.0,
        n_reps     = validate_reps
      )
    }, .options = furrr_options(seed = TRUE), .id = "size_id") %>%
      mutate(n_total = validate_sizes[as.numeric(size_id)], scenario_name = scenario)
  }) %>%
    filter(group != "Overall") %>%
    group_by(scenario_name, n_total, group) %>%
    summarise(power_sim = mean(pval < alpha_bonferroni, na.rm = TRUE), .groups = "drop")
  
  cf_validate <- power_aim1 %>% filter(n_total %in% validate_sizes)
  validate_cmp <- sim_validate %>% left_join(cf_validate, by = c("scenario_name", "n_total", "group")) %>%
    mutate(diff = power_sim - power)
  
  write_tsv(sim_validate, file.path(output_root, "tables/aim1_power_validation_sim.tsv"))
  write_tsv(validate_cmp, file.path(output_root, "tables/aim1_power_validation_compare.tsv"))
}

message(paste("--- AIM 1 COMPLETE --- (Duration:", round(difftime(Sys.time(), start_time_aim1, units = "mins"), 1), "minutes) ---")) 
