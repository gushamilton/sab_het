## =====================================================================
##  03_run_main_binary.R -- Main binary outcome simulations
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")
source("code/functions/01_simulation_helpers.R")
source("code/functions/02_metrics.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table

or_vector <- setNames(subphenotypes$or_arrest_shrunk, subphenotypes$subphenotype)
prevalence <- setNames(subphenotypes$prevalence, subphenotypes$subphenotype)
baseline <- setNames(subphenotypes$baseline_mortality, subphenotypes$subphenotype)

print_effect_directions(or_vector)

sample_sizes <- params$sample_sizes
sample_override <- Sys.getenv("SAMPLE_SIZES", unset = "")
if (nzchar(sample_override)) {
  sample_sizes <- as.integer(strsplit(sample_override, ",")[[1]])
}
accuracy_grid <- params$accuracy_grid
accuracy_override <- Sys.getenv("ACCURACY_GRID", unset = "")
if (nzchar(accuracy_override)) {
  accuracy_grid <- as.numeric(strsplit(accuracy_override, ",")[[1]])
}
alpha_levels <- c(params$alpha_primary, params$alpha_bonferroni)

n_reps <- as.integer(Sys.getenv("N_REPS_MAIN", unset = "2000"))
seed_base <- as.integer(Sys.getenv("SEED_BASE", unset = "123"))

output_root <- file.path("results", "main")
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "data"), showWarnings = FALSE, recursive = TRUE)

run_scenario <- function(n, accuracy, n_reps, seed_base) {
  groups <- names(or_vector)

  map_dfr(seq_len(n_reps), function(rep_id) {
    sim_data <- simulate_binary_trial(
      or_vector = or_vector,
      prevalence = prevalence,
      baseline_mortality = baseline,
      n = n,
      seed = seed_base + rep_id
    )

    observed_group <- misclassify_groups(
      true_group = sim_data$group,
      prevalence = prevalence,
      accuracy = accuracy,
      seed = seed_base + rep_id + 10000
    )

    sim_data$observed_group <- observed_group

    map_dfr(groups, function(group_label) {
      fit <- fit_logistic_or(sim_data %>% filter(observed_group == group_label))
      fit %>% mutate(group = group_label, rep_id = rep_id)
    })
  })
}

scenario_grid <- expand_grid(
  sample_size = sample_sizes,
  accuracy = accuracy_grid
)

all_results <- pmap_dfr(scenario_grid, function(sample_size, accuracy) {
  message(sprintf("Running N=%d, accuracy=%.2f", sample_size, accuracy))
  scenario_results <- run_scenario(sample_size, accuracy, n_reps, seed_base)
  scenario_results %>%
    mutate(sample_size = sample_size, accuracy = accuracy)
})

write_tsv(all_results, file.path(output_root, "data", "main_binary_raw.tsv"))

metric_results <- all_results %>%
  group_by(sample_size, accuracy, group) %>%
  group_modify(~ {
    log_or_true <- log(or_vector[[.y$group]])
    is_null <- abs(log_or_true) < 1e-6

    map_dfr(alpha_levels, function(alpha_level) {
      summarise_metrics(.x, log_or_true = log_or_true, alpha = alpha_level, is_null = is_null) %>%
        mutate(alpha = alpha_level)
    })
  }) %>%
  ungroup() %>%
  mutate(
    scenario = "ARREST_shrunk",
    group_label = subphenotypes$label[match(group, subphenotypes$subphenotype)]
  )

write_tsv(metric_results, file.path(output_root, "data", "main_binary_metrics.tsv"))

table_summary <- metric_results %>%
  filter(sample_size %in% c(5000, 10000), alpha == params$alpha_primary) %>%
  select(sample_size, accuracy, group, group_label, power, type_s, type_m, type1)

write_tsv(
  table_summary,
  file.path(output_root, "tables", "main_binary_summary_N5000_N10000.tsv")
)

power_plot <- metric_results %>%
  filter(alpha == params$alpha_primary, group != "A") %>%
  ggplot(aes(x = sample_size, y = power, color = factor(accuracy))) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_wrap(~ group, scales = "free_y") +
  scale_x_continuous(breaks = sample_sizes) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Power by sample size and accuracy (binary outcome)",
    x = "Sample size",
    y = "Power",
    color = "Accuracy"
  ) +
  theme_minimal()

ggsave(
  file.path(output_root, "plots", "figure1_power_curves_binary.pdf"),
  power_plot,
  width = 11,
  height = 7
)

type_s_plot <- metric_results %>%
  filter(alpha == params$alpha_primary, group != "A", sample_size == 5000) %>%
  ggplot(aes(x = factor(accuracy), y = type_s, fill = group)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Type S error by subphenotype (N=5,000)",
    x = "Accuracy",
    y = "Type S error",
    fill = "Subphenotype"
  ) +
  theme_minimal()

ggsave(
  file.path(output_root, "plots", "figure2_type_s_by_accuracy.pdf"),
  type_s_plot,
  width = 10,
  height = 6
)
