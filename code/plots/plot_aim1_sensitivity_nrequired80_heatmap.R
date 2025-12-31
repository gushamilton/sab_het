# Heatmap of required N for 80% power by cohort and subphenotype

pacman::p_load(tidyverse, scales, viridisLite, ggdist)

theme_set(ggdist::theme_ggdist())

source("code/functions/aim1_closed_form.R")
source("code/common_parameters.R")

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
output_file <- file.path(results_root, "plots", "aim1_sensitivity_nrequired80_heatmap.pdf")

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

n_min <- 100
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
        n_min = n_min,
        n_max = n_max
      )
    ),
    n_required = as.integer(n_required),
    label = ifelse(is.na(n_required), paste0(">", format(n_max, big.mark = ",")), format(n_required, big.mark = ",")),
    cohort = factor(cohort, levels = sort(unique(cohort))),
    group = factor(group, levels = c("A", "B", "C", "D", "E"))
  ) %>%
  filter(group != "A") %>%
  mutate(group = factor(group, levels = c("B", "C", "D", "E")))

plot_title <- if (scenario_set == "supp") {
  "Aim 1 Sensitivity: Required N for 80% Power (Sensitivity)"
} else {
  "Aim 1 Sensitivity: Required N for 80% Power"
}

plot <- ggplot(n_required, aes(x = group, y = cohort, fill = n_required)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_viridis_c(
    option = "C",
    trans = "log10",
    labels = scales::comma,
    na.value = "grey90",
    name = "Required N\n(80% power)"
  ) +
  labs(
    title = plot_title,
    x = "Subphenotype",
    y = "Cohort"
  )

ggsave(output_file, plot, width = 8.5, height = 5.5)
message("Saved: ", output_file)
