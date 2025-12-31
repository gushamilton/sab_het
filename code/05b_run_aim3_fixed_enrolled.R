# -----------------------------------------------------------------------------
# AIM 3 FIXED-ENROLLED REPLICATES: Generate beta replicates at target enrolled N
# -----------------------------------------------------------------------------

message("--- Loading Packages ---")
pacman::p_load(tidyverse, furrr)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))
source("code/functions/functions.R")
source("code/common_parameters.R")
set.seed(20240423)

# Create results directories
scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}
output_root <- file.path("results", scenario_set)
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)

# --- Parameters ---
n_reps_fixed <- as.integer(Sys.getenv("N_REPS_AIM3_FIXED", unset = "500"))
target_enrolled <- as.integer(Sys.getenv("AIM3_TARGET_ENROLLED", unset = "5000"))
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores() - 1))
if (is.na(n_cores) || n_cores < 1) n_cores <- 1
plan(multisession, workers = n_cores)
message(paste("--- Parallel processing enabled on", n_cores, "cores ---"))

params <- get_common_params()

# --- Define OR scenarios ---
or_arrest_raw    <- params$or_arrest
or_arrest_shrunk <- params$or_arrest_shrunk_k05

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

freq_arrest        <- params$freq_arrest
p0_arrest_adjusted <- params$p0_overall
sens_spec_scenarios <- params$sens_spec_scenarios
target_groups_aim3  <- params$target_groups_aim3

all_scenarios <- expand_grid(
  scenario_definitions,
  sens_spec_scenarios,
  target_group = target_groups_aim3
) %>% mutate(seed = 1:n() + 9000)

run_fixed_scenario <- function(scenario_name, or_vector, target_group, sensitivity, specificity, seed) {
  prevalence <- freq_arrest[[target_group]]
  enrol_rate <- prevalence * sensitivity + (1 - prevalence) * (1 - specificity)
  n_screened <- ceiling(target_enrolled / enrol_rate)

  res <- run_enrichment_scenario_sens_spec(
    n_screened = n_screened,
    target_group = target_group,
    sensitivity = sensitivity,
    specificity = specificity,
    or_vector = or_vector,
    p0_vector = p0_arrest_adjusted,
    freq_vector = freq_arrest,
    n_reps = n_reps_fixed,
    base_seed = seed
  )

  true_beta <- log(or_vector[target_group])
  tibble(
    scenario_name = scenario_name,
    target_group = target_group,
    sensitivity = sensitivity,
    specificity = specificity,
    test_type = sens_spec_scenarios %>%
      filter(sensitivity == !!sensitivity, specificity == !!specificity) %>%
      pull(test_type) %>%
      first(),
    n_screened = n_screened,
    target_enrolled = target_enrolled,
    sim_id = seq_along(res$beta_hats),
    true_beta = true_beta,
    empirical_beta = res$beta_hats,
    p_value = res$p_values,
    n_enrolled = res$mean_n_enrolled
  )
}

message("--- STARTING AIM 3 FIXED-ENROLLED REPLICATES ---")
start_time <- Sys.time()

replicates <- pmap_dfr(
  list(
    scenario_name = all_scenarios$scenario_name,
    or_vector = all_scenarios$or_vector,
    target_group = all_scenarios$target_group,
    sensitivity = all_scenarios$sensitivity,
    specificity = all_scenarios$specificity,
    seed = all_scenarios$seed
  ),
  run_fixed_scenario
)

replicates <- replicates %>% mutate(scenario_set = scenario_set)
output_file <- file.path(output_root, "tables", "aim3_fixed_enrolled_replicates.tsv.gz")
write_tsv(replicates, output_file)

message("Saved: ", output_file)
message(paste("--- AIM 3 FIXED-ENROLLED COMPLETE --- (Duration:",
              round(difftime(Sys.time(), start_time, units = "mins"), 1), "minutes) ---"))
