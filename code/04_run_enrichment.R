## =====================================================================
##  04_run_enrichment.R -- Supplement 1 enrichment calculations
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table

or_vector <- setNames(subphenotypes$or_arrest_shrunk, subphenotypes$subphenotype)
prevalence <- setNames(subphenotypes$prevalence, subphenotypes$subphenotype)
baseline <- setNames(subphenotypes$baseline_mortality, subphenotypes$subphenotype)

output_root <- file.path("results", "supp")
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)

get_or_p1 <- function(or_value, p0) {
  (or_value * p0) / (1 - p0 + (or_value * p0))
}

calculate_mixed_cohort <- function(target_group, sensitivity, specificity,
                                   or_vector, p0_vector, prevalence) {
  p_vec <- prevalence / sum(prevalence)
  p_target <- p_vec[target_group]
  p_non_target <- 1 - p_target

  non_target_groups <- names(p_vec)[names(p_vec) != target_group]
  weights_mix <- p_vec[non_target_groups] / p_non_target

  p0_target <- p0_vector[target_group]
  p1_target <- get_or_p1(or_vector[target_group], p0_target)

  p0_mix <- sum(p0_vector[non_target_groups] * weights_mix)
  p1_mix <- sum(get_or_p1(or_vector[non_target_groups], p0_vector[non_target_groups]) * weights_mix)

  tp_rate <- sensitivity * p_target
  fp_rate <- (1 - specificity) * p_non_target
  enrol_rate <- tp_rate + fp_rate

  if (enrol_rate == 0) {
    return(list(beta_obs = 0, p0_obs = 0, p1_obs = 0, enrol_rate = 0, p_true = NA_real_))
  }

  p_true <- tp_rate / enrol_rate
  p0_obs <- p_true * p0_target + (1 - p_true) * p0_mix
  p1_obs <- p_true * p1_target + (1 - p_true) * p1_mix

  if (p0_obs %in% c(0, 1) || p1_obs %in% c(0, 1)) {
    beta_obs <- 0
  } else {
    beta_obs <- log((p1_obs / (1 - p1_obs)) / (p0_obs / (1 - p0_obs)))
  }

  list(beta_obs = beta_obs, p0_obs = p0_obs, p1_obs = p1_obs, enrol_rate = enrol_rate, p_true = p_true)
}

calc_required_n <- function(beta_obs, p0_obs, p1_obs, alpha = 0.05, power = 0.80) {
  if (abs(beta_obs) < 1e-6 || p0_obs %in% c(0, 1) || p1_obs %in% c(0, 1)) {
    return(Inf)
  }
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta <- qnorm(power)
  var_term <- (1 / (p0_obs * (1 - p0_obs))) + (1 / (p1_obs * (1 - p1_obs)))
  2 * (z_alpha + z_beta)^2 * var_term / (beta_obs^2)
}

scenario_grid <- expand_grid(
  params$enrichment_tests,
  target_group = subphenotypes$subphenotype[subphenotypes$subphenotype != "A"]
)

enrichment_results <- pmap_dfr(scenario_grid, function(test_type, sensitivity, specificity, target_group) {
  cohort_props <- calculate_mixed_cohort(
    target_group = target_group,
    sensitivity = sensitivity,
    specificity = specificity,
    or_vector = or_vector,
    p0_vector = baseline,
    prevalence = prevalence
  )

  n_total <- calc_required_n(cohort_props$beta_obs, cohort_props$p0_obs, cohort_props$p1_obs)
  nns <- n_total / cohort_props$enrol_rate

  tibble(
    test_type = test_type,
    sensitivity = sensitivity,
    specificity = specificity,
    target_group = target_group,
    n_required = ceiling(n_total),
    nns = ceiling(nns),
    enrol_rate = cohort_props$enrol_rate,
    true_log_or = log(or_vector[target_group]),
    observed_log_or = cohort_props$beta_obs,
    bias_log_or = cohort_props$beta_obs - log(or_vector[target_group])
  )
})

write_tsv(
  enrichment_results,
  file.path(output_root, "tables", "supp1_enrichment_nns_nnr.tsv")
)

effect_plot <- enrichment_results %>%
  mutate(
    observed_or = exp(observed_log_or),
    true_or = exp(true_log_or)
  ) %>%
  ggplot(aes(x = test_type, y = observed_or, color = target_group)) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_hline(aes(yintercept = true_or, color = target_group), linetype = "dashed") +
  scale_y_continuous(trans = "log10") +
  labs(
    title = "Supplement 1: Enrichment bias (expected OR, N=5,000)",
    x = "Test type",
    y = "Observed OR (log scale)",
    color = "Subphenotype"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  file.path(output_root, "plots", "figureS1_enrichment_bias.pdf"),
  effect_plot,
  width = 11,
  height = 6
)
