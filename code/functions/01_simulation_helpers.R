## =====================================================================
##  01_simulation_helpers.R -- Simulation utilities
## =====================================================================

library(dplyr)
library(tibble)
library(purrr)
library(broom)

calculate_treated_risk <- function(or_value, p0) {
  (or_value * p0) / (1 - p0 + (or_value * p0))
}

simulate_binary_trial <- function(or_vector, prevalence, baseline_mortality, n, seed = 123) {
  set.seed(seed)

  group_labels <- names(or_vector)
  prevalence <- prevalence[group_labels]
  baseline_mortality <- baseline_mortality[group_labels]

  tibble(
    id = seq_len(n),
    group = sample(group_labels, size = n, replace = TRUE, prob = prevalence),
    treatment = rbinom(n, 1, 0.5)
  ) %>%
    mutate(
      or_value = or_vector[group],
      p0 = baseline_mortality[group],
      p1 = calculate_treated_risk(or_value, p0),
      outcome = rbinom(n, 1, if_else(treatment == 1, p1, p0))
    ) %>%
    select(id, group, treatment, outcome)
}

derive_survivor_split <- function(ordinal_baseline_probs) {
  survivor_probs <- ordinal_baseline_probs[-1]
  survivor_probs / sum(survivor_probs)
}

build_ordinal_probs_from_death <- function(p_dead, survivor_split) {
  c(p_dead, (1 - p_dead) * survivor_split)
}

get_po_ordinal_probs <- function(p_dead, survivor_split, or_value) {
  control_probs <- build_ordinal_probs_from_death(p_dead, survivor_split)
  cumulative_probs <- cumsum(control_probs)[-length(control_probs)]
  intercepts <- qlogis(cumulative_probs)

  to_probs <- function(eta) {
    cumulative <- plogis(eta)
    probs <- diff(c(0, cumulative, 1))
    probs <- pmax(probs, 0)
    probs / sum(probs)
  }

  list(
    control = to_probs(intercepts),
    treatment = to_probs(intercepts + log(or_value))
  )
}

misclassify_groups <- function(true_group, prevalence, accuracy, seed = 123) {
  set.seed(seed)
  group_labels <- names(prevalence)
  prevalence <- prevalence[group_labels]

  correct <- runif(length(true_group)) < accuracy
  observed <- true_group

  if (any(!correct)) {
    observed[!correct] <- vapply(
      true_group[!correct],
      function(group_label) {
        other_groups <- group_labels[group_labels != group_label]
        other_weights <- prevalence[other_groups]
        sample(other_groups, size = 1, prob = other_weights)
      },
      character(1)
    )
  }

  observed
}

fit_logistic_or <- function(data) {
  if (nrow(data) < 10 || length(unique(data$treatment)) < 2 || length(unique(data$outcome)) < 2) {
    return(tibble(log_or_hat = NA_real_, se = NA_real_))
  }

  model <- tryCatch(
    glm(outcome ~ treatment, data = data, family = binomial()),
    error = function(e) NULL
  )
  if (is.null(model)) return(tibble(log_or_hat = NA_real_, se = NA_real_))

  tidy_model <- broom::tidy(model)
  treatment_row <- tidy_model %>% filter(term == "treatment")
  if (nrow(treatment_row) != 1) {
    return(tibble(log_or_hat = NA_real_, se = NA_real_))
  }

  tibble(log_or_hat = treatment_row$estimate, se = treatment_row$std.error)
}

simulate_ordinal_trial <- function(or_vector, prevalence, baseline_mortality, survivor_split, n, seed = 123) {
  set.seed(seed)

  group_labels <- names(or_vector)
  prevalence <- prevalence[group_labels]
  baseline_mortality <- baseline_mortality[group_labels]

  tibble(
    id = seq_len(n),
    group = sample(group_labels, size = n, replace = TRUE, prob = prevalence),
    treatment = rbinom(n, 1, 0.5)
  ) %>%
    mutate(
      probs = map2(group, treatment, function(group_label, trt) {
        po_probs <- get_po_ordinal_probs(
          p_dead = baseline_mortality[[group_label]],
          survivor_split = survivor_split,
          or_value = or_vector[[group_label]]
        )
        if (trt == 1) {
          po_probs$treatment
        } else {
          po_probs$control
        }
      }),
      outcome = map_int(probs, ~ sample(seq_along(.x), size = 1, prob = .x))
    ) %>%
    select(id, group, treatment, outcome)
}

fit_polr_or <- function(data) {
  if (nrow(data) < 20 || length(unique(data$treatment)) < 2 || length(unique(data$outcome)) < 2) {
    return(tibble(log_or_hat = NA_real_, se = NA_real_))
  }

  data <- data %>%
    mutate(outcome = ordered(outcome, levels = sort(unique(outcome))))

  model <- tryCatch(
    MASS::polr(outcome ~ treatment, data = data, method = "logistic", Hess = TRUE),
    error = function(e) NULL
  )
  if (is.null(model)) return(tibble(log_or_hat = NA_real_, se = NA_real_))

  coef_table <- summary(model)$coefficients
  if (!"treatment" %in% rownames(coef_table)) {
    return(tibble(log_or_hat = NA_real_, se = NA_real_))
  }

  # MASS::polr uses logit(P(Y <= k)) = theta_k - beta * x.
  # Our ordinal simulation uses log_or as the shift in logit(P(Y <= k)) for treatment,
  # so the comparable log_or is -beta.
  tibble(
    log_or_hat = -coef_table["treatment", "Value"],
    se = coef_table["treatment", "Std. Error"]
  )
}
