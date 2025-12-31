# Plot baseline mortality by subphenotype (cohort colors)

pacman::p_load(tidyverse, scales, ggdist)

theme_set(ggdist::theme_ggdist() + theme(legend.position = "bottom"))

input_file <- "code/tables/mortality_by_cohort_with_n.tsv"
output_file <- "results/main/plots/baseline_mortality_by_cohort.pdf"

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

mortality <- readr::read_tsv(input_file, show_col_types = FALSE)

required_cols <- c("cohort", "subphenotype", "mortality")
missing_cols <- setdiff(required_cols, names(mortality))
if (length(missing_cols) > 0) {
  stop("Input is missing columns: ", paste(missing_cols, collapse = ", "))
}

mortality <- mortality %>%
  mutate(subphenotype = factor(subphenotype, levels = c("A", "B", "C", "D", "E")))

plot <- ggplot(mortality, aes(x = subphenotype, y = mortality, color = cohort)) +
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 2, alpha = 0.9) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Baseline Mortality by Subphenotype",
    x = "Subphenotype",
    y = "Baseline Mortality",
    color = "Cohort"
  )

ggsave(output_file, plot, width = 8, height = 5)
message("Saved: ", output_file)
