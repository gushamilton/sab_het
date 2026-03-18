## =====================================================================
##  05_run_cohort_variation.R -- Supplement 2 cohort sensitivity
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")
source("code/functions/01_simulation_helpers.R")
source("code/functions/02_metrics.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table

or_vector <- setNames(subphenotypes$or_arrest_shrunk, subphenotypes$subphenotype)

output_root <- file.path("results", "supp")
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "data"), showWarnings = FALSE, recursive = TRUE)

cohort_data <- readr::read_tsv(params$cohort_mortality_path, show_col_types = FALSE)

cohort_summary <- cohort_data %>%
  mutate(prevalence = n_subpheno / n_cohort) %>%
  select(cohort, subphenotype, n_subpheno, n_cohort, prevalence, mortality)

write_tsv(cohort_summary, file.path(output_root, "tables", "supp2_cohort_prevalence_mortality.tsv"))

calc_required_n <- function(or_value, p0, alpha = 0.05, power = 0.80) {
  if (abs(or_value - 1) < 1e-6 || p0 %in% c(0, 1)) return(Inf)
  p1 <- calculate_treated_risk(or_value, p0)
  if (p1 %in% c(0, 1)) return(Inf)

  z_alpha <- qnorm(1 - alpha / 2)
  z_beta <- qnorm(power)
  var_term <- (1 / (p0 * (1 - p0))) + (1 / (p1 * (1 - p1)))
  2 * (z_alpha + z_beta)^2 * var_term / (log(or_value)^2)
}

required_n_table <- cohort_summary %>%
  mutate(
    or_value = or_vector[subphenotype],
    n_required = map2_dbl(or_value, mortality, calc_required_n)
  ) %>%
  select(cohort, subphenotype, n_required)

write_tsv(required_n_table, file.path(output_root, "tables", "supp2_required_n_by_cohort.tsv"))

heatmap_plot <- required_n_table %>%
  mutate(n_required = if_else(is.finite(n_required), n_required, NA_real_)) %>%
  ggplot(aes(x = cohort, y = subphenotype, fill = n_required)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C", trans = "log10", na.value = "grey90") +
  labs(
    title = "Supplement 2: Required N for 80% power (perfect classification)",
    x = "Cohort",
    y = "Subphenotype",
    fill = "N required"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  file.path(output_root, "plots", "figureS2_required_n_heatmap.pdf"),
  heatmap_plot,
  width = 10,
  height = 6
)

run_cohort_sim <- function(prevalence, baseline, n, n_reps, seed_base = 500) {
  groups <- names(or_vector)
  map_dfr(seq_len(n_reps), function(rep_id) {
    sim_data <- simulate_binary_trial(
      or_vector = or_vector,
      prevalence = prevalence,
      baseline_mortality = baseline,
      n = n,
      seed = seed_base + rep_id
    )
    sim_data$observed_group <- sim_data$group

    map_dfr(groups, function(group_label) {
      fit <- fit_logistic_or(sim_data %>% filter(observed_group == group_label))
      fit %>% mutate(group = group_label, rep_id = rep_id)
    })
  })
}

n_reps <- as.integer(Sys.getenv("N_REPS_COHORT", unset = "1000"))
fixed_n <- as.integer(Sys.getenv("COHORT_FIXED_N", unset = "5000"))
save_cohort_raw <- tolower(Sys.getenv("SAVE_COHORT_RAW", unset = "0")) %in% c("1", "true", "yes")

cohort_names <- unique(cohort_summary$cohort)

cohort_results <- map(cohort_names, function(cohort_name) {
  cohort_slice <- cohort_summary %>% filter(cohort == cohort_name)
  prevalence <- setNames(cohort_slice$prevalence, cohort_slice$subphenotype)
  baseline <- setNames(cohort_slice$mortality, cohort_slice$subphenotype)

  sim_results <- run_cohort_sim(prevalence, baseline, fixed_n, n_reps)

  metrics <- sim_results %>%
    group_by(group) %>%
    group_modify(~ {
      log_or_true <- log(or_vector[[.y$group]])
      is_null <- abs(log_or_true) < 1e-6
      summarise_metrics(.x, log_or_true, params$alpha_primary, is_null = is_null)
    }) %>%
    ungroup() %>%
    mutate(cohort = cohort_name)

  list(
    metrics = metrics,
    raw = sim_results %>% mutate(cohort = cohort_name)
  )
})

cohort_power_results <- bind_rows(map(cohort_results, "metrics")) %>%
  mutate(sample_size = fixed_n)

write_tsv(cohort_power_results, file.path(output_root, "tables", "supp2_power_typeS_fixedN.tsv"))

if (save_cohort_raw) {
  raw_cohort <- bind_rows(map(cohort_results, "raw"))
  raw_path <- file.path(output_root, "data", sprintf("supp2_cohort_raw_fixedN%d.tsv.gz", fixed_n))
  write_tsv(raw_cohort, raw_path)
}
