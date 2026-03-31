## =====================================================================
##  01_parameters.R -- Trial parameters and cohort inputs
## =====================================================================

library(tibble)

subphenotype_table <- tibble(
  subphenotype = c("A", "B", "C", "D", "E"),
  label = c(
    "Elderly/comorbid",
    "Nosocomial IV",
    "Metastatic",
    "CKD",
    "IDU"
  ),
  prevalence = c(0.155, 0.134, 0.356, 0.178, 0.178),
  baseline_mortality = c(0.217, 0.075, 0.192, 0.182, 0.029),
  or_arrest_raw = c(1.00, 18.8, 0.79, 1.42, 0.31),
  or_arrest_shrunk = c(1.00, 4.34, 0.89, 1.19, 0.56)
)

sample_sizes <- c(500, 1000, 1500, 2000, 3000, 5000, 10000, 20000)
ordinal_sample_sizes <- c(500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000, 7500, 10000, 15000, 20000)
accuracy_grid <- seq(0.70, 1.00, by = 0.05)
alpha_primary <- 0.05
alpha_bonferroni <- 0.01

# Five-level sensitivity baseline retained for the supplementary 5-point
# comparison; the main manuscript-facing ordinal analyses now default to a
# six-level setup defined within the analysis scripts.
ordinal_baseline <- tibble(
  level = 1:5,
  label = c(
    "Dead",
    "ICU/ventilated",
    "Still hospitalised",
    "Discharged with complications",
    "Discharged well"
  ),
  prevalence = c(0.15, 0.05, 0.10, 0.20, 0.50)
)

enrichment_tests <- tibble(
  test_type = c(
    "Perfect",
    "Near-perfect",
    "High/High",
    "Balanced",
    "High Sens/Low Spec",
    "Low Sens/High Spec"
  ),
  sensitivity = c(1.00, 0.99, 0.95, 0.80, 0.95, 0.70),
  specificity = c(1.00, 0.99, 0.95, 0.80, 0.70, 0.95)
)

cohort_mortality_path <- "code/tables/mortality_by_cohort_with_n.tsv"

print_effect_directions <- function(or_vector) {
  direction_label <- function(or_value) {
    if (is.na(or_value)) return("NA")
    if (abs(or_value - 1) < 1e-6) return("NULL")
    if (or_value > 1) return("HARMFUL")
    "PROTECTIVE"
  }

  message("Effect direction sanity check:")
  for (group in names(or_vector)) {
    message(
      sprintf(
        "  Subphenotype %s: %s (OR = %.2f)",
        group,
        direction_label(or_vector[[group]]),
        or_vector[[group]]
      )
    )
  }
}

get_parameters <- function() {
  list(
    subphenotype_table = subphenotype_table,
    sample_sizes = sample_sizes,
    ordinal_sample_sizes = ordinal_sample_sizes,
    accuracy_grid = accuracy_grid,
    alpha_primary = alpha_primary,
    alpha_bonferroni = alpha_bonferroni,
    ordinal_baseline = ordinal_baseline,
    enrichment_tests = enrichment_tests,
    cohort_mortality_path = cohort_mortality_path
  )
}
