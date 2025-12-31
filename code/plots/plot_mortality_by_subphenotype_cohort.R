# Plot mortality by subphenotype with cohort dots

pacman::p_load(tidyverse, scales, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

input_file <- "code/tables/mortality_by_cohort.tsv"
output_file <- "results/main/plots/mortality_by_subphenotype_cohort.pdf"

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

mortality <- readr::read_tsv(input_file, show_col_types = FALSE)

# Expected columns: cohort, n, subphenotype, mortality
missing_cols <- setdiff(c("cohort", "n", "subphenotype", "mortality"), names(mortality))
if (length(missing_cols) > 0) {
  stop("Input is missing columns: ", paste(missing_cols, collapse = ", "))
}

mortality <- mortality %>%
  mutate(
    subphenotype = factor(subphenotype, levels = c("A", "B", "C", "D", "E")),
    cohort = factor(cohort),
    cohort_label = paste0(cohort, "\n(n=", n, ")")
  ) %>%
  filter(subphenotype != "A") %>%
  mutate(subphenotype = factor(subphenotype, levels = c("B", "C", "D", "E")))

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
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 2, alpha = 0.9) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Mortality by Subphenotype and Cohort",
    x = "Subphenotype",
    y = "Mortality",
    color = "Cohort"
  )

ggsave(output_file, gg, width = 8, height = 5)
message("Saved: ", output_file)
