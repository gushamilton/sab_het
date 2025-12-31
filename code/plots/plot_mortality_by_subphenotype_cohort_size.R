# Plot mortality by subphenotype with cohort dots sized by subphenotype n

pacman::p_load(tidyverse, scales, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

input_file <- "code/tables/mortality_by_cohort_with_n.tsv"
output_file <- "results/main/plots/mortality_by_subphenotype_cohort_size.pdf"

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

mortality <- readr::read_tsv(input_file, show_col_types = FALSE)

# Expected columns: cohort, n_cohort, subphenotype, n_subpheno, mortality
missing_cols <- setdiff(c("cohort", "n_cohort", "subphenotype", "n_subpheno", "mortality"), names(mortality))
if (length(missing_cols) > 0) {
  stop("Input is missing columns: ", paste(missing_cols, collapse = ", "))
}

mortality <- mortality %>%
  mutate(
    subphenotype = factor(subphenotype, levels = c("A", "B", "C", "D", "E")),
    cohort = factor(cohort),
    cohort_label = paste0(cohort, "\n(n=", n_cohort, ")")
  )

gg <- ggplot(mortality, aes(x = subphenotype, y = mortality, color = cohort_label)) +
  geom_boxplot(
    aes(x = subphenotype, y = mortality),
    inherit.aes = FALSE,
    fill = "grey80",
    color = "grey60",
    alpha = 0.5,
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_point(
    aes(size = n_subpheno),
    position = position_jitter(width = 0.1, height = 0),
    alpha = 0.9
  ) +
  scale_size_continuous(range = c(2, 6)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Mortality by Subphenotype and Cohort (Point Size = Subphenotype N)",
    x = "Subphenotype",
    y = "Mortality",
    color = "Cohort",
    size = "Subphenotype N"
  )

ggsave(output_file, gg, width = 8, height = 5)
message("Saved: ", output_file)
