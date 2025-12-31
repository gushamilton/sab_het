# Plot Aim 1 sensitivity: exact required N for 80% power by cohort and subgroup

pacman::p_load(tidyverse, scales, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

source("code/functions/aim1_closed_form.R")
source("code/common_parameters.R")

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
output_file <- file.path(results_root, "plots", "aim1_sensitivity_nrequired_by_cohort.pdf")

cohort_input <- "code/tables/mortality_by_cohort_with_n.tsv"
if (!file.exists(cohort_input)) {
  stop("Missing input: ", cohort_input)
}

params <- get_common_params()
or_vector <- if (scenario_set == "main") params$or_arrest else params$or_arrest_shrunk_k05

cohort_params <- readr::read_tsv(cohort_input, show_col_types = FALSE) %>%
  mutate(
    group = subphenotype,
    p0 = mortality,
    prevalence = n_subpheno / n_cohort
  ) %>%
  select(cohort, n_cohort, group, p0, prevalence)

n_max <- 200000

n_required <- cohort_params %>%
  mutate(
    n_required = purrr::pmap_dbl(
      list(p0, group, prevalence),
      ~ solve_n_required(
        target_power = 0.8,
        p0 = ..1,
        or = or_vector[..2],
        prevalence = ..3,
        alpha = 0.05 / 5,
        n_min = 500,
        n_max = n_max
      )
    )
  )

plot <- n_required %>%
  ggplot(aes(x = group, y = n_required, color = cohort, group = cohort)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "Aim 1 Sensitivity: Required N for 80% Power",
    x = "Subphenotype",
    y = "Required Total Sample Size (log10)",
    color = "Cohort"
  ) +
  theme(axis.text.x = element_text(angle = 0))

ggsave(output_file, plot, width = 9, height = 6)
message("Saved: ", output_file)
