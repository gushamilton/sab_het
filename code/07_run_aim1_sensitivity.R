# -----------------------------------------------------------------------------
# 07_run_aim1_sensitivity.R - AIM 1 sensitivity: baseline risk + prevalence
# -----------------------------------------------------------------------------

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
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)

# --- 2. DEFINE SIMULATION PARAMETERS ---
message("--- Defining Simulation Parameters ---")

# --- Repetitions and Parallel Cores ---
n_reps <- as.integer(Sys.getenv("N_REPS_AIM1_SENS", unset = "500"))
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores() - 1))
if (is.na(n_cores) || n_cores < 1) n_cores <- 1
plan(multisession, workers = n_cores)
message(paste("--- Parallel processing enabled on", n_cores, "cores ---"))

params <- get_common_params()

# --- Define OR scenario ---
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

# --- Cohort-specific baseline risks and prevalences ---
cohort_input <- "code/tables/mortality_by_cohort_with_n.tsv"
if (!file.exists(cohort_input)) {
  stop("Missing input: ", cohort_input)
}

cohort_data <- readr::read_tsv(cohort_input, show_col_types = FALSE)
required_cols <- c("cohort", "n_cohort", "subphenotype", "n_subpheno", "mortality")
missing_cols <- setdiff(required_cols, names(cohort_data))
if (length(missing_cols) > 0) {
  stop("Cohort input missing columns: ", paste(missing_cols, collapse = ", "))
}

cohort_params <- cohort_data %>%
  group_by(cohort, n_cohort) %>%
  summarise(
    p0_vector = list(setNames(mortality, subphenotype)),
    freq_vector = list(setNames(n_subpheno / sum(n_subpheno), subphenotype)),
    .groups = "drop"
  )

cohort_params_long <- cohort_data %>%
  mutate(
    group = subphenotype,
    p0 = mortality,
    prevalence = n_subpheno / n_cohort
  ) %>%
  select(cohort, n_cohort, group, p0, prevalence)

alpha_bonferroni <- 0.05 / 5
sample_sizes <- seq(100, 20000, by = 100)
aim1_mode <- Sys.getenv("AIM1_MODE", unset = "closed_form")

message("\n--- STARTING AIM 1 SENSITIVITY ---")
start_time <- Sys.time()

all_combos <- expand_grid(
  scenario_definitions,
  cohort_params
)

if (aim1_mode == "closed_form") {
  scenario_long <- scenario_definitions %>%
    mutate(or_list = map(or_vector, ~tibble(group = names(.x), or = as.numeric(.x)))) %>%
    select(-or_vector) %>%
    unnest(or_list)

  summary <- scenario_long %>%
    inner_join(cohort_params_long, by = "group") %>%
    crossing(n_total = sample_sizes) %>%
    mutate(
      power = calc_power_closed_form(
        n_total = n_total,
        p0 = p0,
        or = or,
        prevalence = prevalence,
        alpha = alpha_bonferroni
      )
    ) %>%
    select(scenario_name, cohort, n_cohort, n_total, group, power)
  results <- NULL
} else {
  results <- pmap_dfr(list(
      scenario_name = all_combos$scenario_name,
      or_vector = all_combos$or_vector,
      cohort = all_combos$cohort,
      n_cohort = all_combos$n_cohort,
      p0_vector = all_combos$p0_vector,
      freq_vector = all_combos$freq_vector
    ), function(scenario_name, or_vector, cohort, n_cohort, p0_vector, freq_vector) {
      message(paste("  Running:", scenario_name, "cohort=", cohort))
      future_map_dfr(sample_sizes, function(ss) {
        replicate_sims(
          or_vector  = or_vector,
          freq_vector = freq_vector,
          p0_vector  = p0_vector,
          n          = ss,
          accuracy   = 1.0,
          n_reps     = n_reps
        )
      }, .options = furrr_options(seed = TRUE), .id = "size_id") %>%
        mutate(
          n_total = sample_sizes[as.numeric(size_id)],
          scenario_name = scenario_name,
          cohort = cohort,
          n_cohort = n_cohort
        )
    })
  
  summary <- results %>%
    filter(group != "Overall") %>%
    group_by(scenario_name, cohort, n_cohort, n_total, group) %>%
    summarise(power = mean(pval < alpha_bonferroni, na.rm = TRUE), .groups = "drop")
}

summary <- summary %>% mutate(scenario_set = scenario_set)
write_tsv(summary, file.path(output_root, "tables/aim1_sensitivity_summary.tsv"))

if (!is.null(results)) {
  results <- results %>% mutate(scenario_set = scenario_set)
  write_tsv(results, file.path(output_root, "tables/aim1_sensitivity_raw.tsv.gz"))
}

message(paste("--- AIM 1 SENSITIVITY COMPLETE --- (Duration:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "minutes) ---"))
