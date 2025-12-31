# Plot subphenotype prevalence by cohort (cohort colors)

pacman::p_load(tidyverse, scales, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

input_file <- "code/tables/mortality_by_cohort_with_n.tsv"
output_file <- "results/main/plots/subphenotype_prevalence_by_cohort.pdf"

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

cohort_data <- readr::read_tsv(input_file, show_col_types = FALSE)

required_cols <- c("cohort", "n_cohort", "subphenotype", "n_subpheno")
missing_cols <- setdiff(required_cols, names(cohort_data))
if (length(missing_cols) > 0) {
  stop("Input is missing columns: ", paste(missing_cols, collapse = ", "))
}

cohort_data <- cohort_data %>%
  mutate(
    subphenotype = factor(subphenotype, levels = c("A", "B", "C", "D", "E")),
    prevalence = n_subpheno / n_cohort,
    cohort_label = paste0(cohort, "\n(n=", n_cohort, ")")
  ) %>%
  filter(subphenotype != "A") %>%
  mutate(subphenotype = factor(subphenotype, levels = c("B", "C", "D", "E")))

plot <- ggplot(cohort_data, aes(x = subphenotype, y = prevalence, color = cohort_label)) +
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
    title = "Subphenotype Prevalence by Cohort",
    x = "Subphenotype",
    y = "Prevalence",
    color = "Cohort"
  )

ggsave(output_file, plot, width = 8, height = 5)
message("Saved: ", output_file)
