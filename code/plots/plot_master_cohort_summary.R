# Master plot: prevalence, mortality, and N-required heatmap by cohort

pacman::p_load(tidyverse, scales, viridisLite, patchwork, ggdist)

theme_set(ggdist::theme_ggdist())

source("code/functions/aim1_closed_form.R")
source("code/common_parameters.R")

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

output_file <- file.path("results", scenario_set, "plots", "cohort_summary_master.pdf")
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

cohort_input <- "code/tables/mortality_by_cohort_with_n.tsv"
mortality_input <- "code/tables/mortality_by_cohort.tsv"

if (!file.exists(cohort_input)) {
  stop("Missing input: ", cohort_input)
}
if (!file.exists(mortality_input)) {
  stop("Missing input: ", mortality_input)
}

cohort_data <- readr::read_tsv(cohort_input, show_col_types = FALSE)
mortality <- readr::read_tsv(mortality_input, show_col_types = FALSE)

cohort_data <- cohort_data %>%
  mutate(
    subphenotype = factor(subphenotype, levels = c("A", "B", "C", "D", "E")),
    prevalence = n_subpheno / n_cohort,
    cohort_label = paste0(cohort, "\n(n=", n_cohort, ")")
  )

mortality <- mortality %>%
  mutate(
    subphenotype = factor(subphenotype, levels = c("A", "B", "C", "D", "E")),
    cohort_label = paste0(cohort, "\n(n=", n, ")")
  )

prevalence_plot <- ggplot(cohort_data, aes(x = subphenotype, y = prevalence, color = cohort_label)) +
  geom_boxplot(
    aes(x = subphenotype, y = prevalence),
    inherit.aes = FALSE,
    fill = "grey80",
    color = "grey60",
    alpha = 0.5,
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 2, alpha = 0.9) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Subphenotype",
    y = "Prevalence",
    color = "Cohort"
  ) +
  theme(legend.position = "none")

mortality_plot <- ggplot(mortality, aes(x = subphenotype, y = mortality, color = cohort_label)) +
  geom_boxplot(
    aes(x = subphenotype, y = mortality),
    inherit.aes = FALSE,
    fill = "grey80",
    color = "grey60",
    alpha = 0.5,
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 2, alpha = 0.9) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Subphenotype",
    y = "Mortality",
    color = "Cohort"
  )

params <- get_common_params()
or_vector <- if (scenario_set == "main") params$or_arrest else params$or_arrest_shrunk_k05

n_min <- 100
n_max <- 200000

baseline_mortality <- readr::read_tsv(mortality_input, show_col_types = FALSE) %>%
  select(cohort, subphenotype, mortality) %>%
  distinct()

n_required <- baseline_mortality %>%
  mutate(
    group = subphenotype,
    cohort = as.character(cohort)
  ) %>%
  left_join(
    cohort_data %>%
      select(cohort, subphenotype, prevalence) %>%
      distinct(),
    by = c("cohort", "subphenotype")
  ) %>%
  mutate(
    n_required = purrr::pmap_dbl(
      list(mortality, group, prevalence),
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
    label = dplyr::case_when(
      is.na(n_required) ~ ">100,000",
      n_required > 100000 ~ ">100,000",
      TRUE ~ format(n_required, big.mark = ",")
    ),
    cohort = factor(cohort, levels = sort(unique(cohort))),
    group = factor(group, levels = c("A", "B", "C", "D", "E"))
  )

heatmap_plot <- ggplot(n_required, aes(x = group, y = cohort, fill = n_required)) +
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
    x = "Subphenotype",
    y = "Cohort"
  )

final_plot <- (prevalence_plot / mortality_plot / heatmap_plot) +
  plot_annotation(tag_levels = "A")

ggsave(output_file, final_plot, width = 9, height = 13)
message("Saved: ", output_file)
