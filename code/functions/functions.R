
# -----------------------------------------------------------------------------
# Helper function to fit GLM and extract results safely
# Handles cases with insufficient data or model convergence errors
# -----------------------------------------------------------------------------
fit_glm_safe <- function(data) {
  pacman::p_load(logistf) # Ensure logistf is available
  if (nrow(data) < 10 || length(unique(data$treatment)) < 2 || length(unique(data$success)) < 2) {
    return(tibble(beta = NA_real_, se = NA_real_, pval = NA_real_))
  }
  # Try Firth's penalized logistic regression first
  tryCatch({
    data_glm <- data %>% mutate(treatment = factor(treatment))
    model <- logistf::logistf(success ~ treatment, data = data_glm)
    # Find the treatment coefficient (should be 'treatment1')
    idx <- which(names(model$coefficients) == "treatment1")
    if (length(idx) == 1) {
      tibble(
        beta = model$coefficients[[idx]],
        se = sqrt(diag(model$var))[idx],
        pval = model$prob[idx]
      )
    } else {
      tibble(beta = NA_real_, se = NA_real_, pval = NA_real_)
    }
  }, error = function(e) {
    # Fallback to standard GLM if logistf fails
    tryCatch({
      data_glm <- data %>% mutate(treatment = factor(treatment))
      model <- glm(success ~ treatment, data = data_glm, family = binomial())
      tidy_model <- broom::tidy(model)
      treatment_row <- tidy_model %>% filter(str_detect(term, "^treatment"))
      if (nrow(treatment_row) == 1) {
        tibble(
          beta = treatment_row$estimate,
          se = treatment_row$std.error,
          pval = treatment_row$p.value
        )
      } else {
        tibble(beta = NA_real_, se = NA_real_, pval = NA_real_)
      }
    }, error = function(e) {
      tibble(beta = NA_real_, se = NA_real_, pval = NA_real_)
    })
  })
}

# -----------------------------------------------------------------------------
# Function to misclassify group
# -----------------------------------------------------------------------------
misclassify_group <- function(sim_data, accuracy, freq_vector, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  n <- nrow(sim_data)
  group_labels <- names(freq_vector)
  if(is.null(group_labels)) group_labels <- LETTERS[1:length(freq_vector)]
  names(freq_vector) <- group_labels # Ensure names are set
  freq_normalized <- freq_vector / sum(freq_vector)

  # Ensure input data has 'group' column
  if (!"group" %in% names(sim_data)) stop("Input data must have a 'group' column.")

  sim_data_with_assigned <- sim_data |>
    dplyr::mutate(
      is_correct = rbinom(dplyr::n(), 1, accuracy),
      assigned_group = dplyr::case_when(
        is_correct == 1 ~ group,
        # Ensure sampling uses the normalized, named vector
        TRUE ~ sample(names(freq_normalized), dplyr::n(), replace = TRUE, prob = freq_normalized)
      )
    ) %>%
    select(-is_correct) # Remove intermediate column
  return(sim_data_with_assigned)
}


    #' Simulate Trial Data with Subgroups and Heterogeneous Effects
    #'
    #' Generates individual patient data for a simulated two-arm trial, 
    #' incorporating subgroup membership, subgroup-specific baseline risks (p0),
    #' and subgroup-specific treatment effects (odds ratios).
    #'
    #' @param or_vector Numeric vector of odds ratios for the treatment effect in each subgroup.
    #' @param freq_vector Numeric vector of frequencies (proportions) for each subgroup. Must sum to 1 or will be normalized.
    #' @param n Integer, the total number of participants to simulate.
    #' @param p0_vector Numeric vector of baseline event probabilities in the control group for each subgroup. If NULL, defaults to 0.3 for all.
    #' @param seed Integer, seed for random number generation for reproducibility.
    #'
    #' @return A tibble with columns: id, group (true subgroup), treatment (0/1), success (0/1 outcome).
    #' @export
    #'
    #' @examples
    #' or_example <- c(1, 18.8, 0.79, 1.4, 0.3)
    #' freq_example <- c(0.15, 0.13, 0.36, 0.18, 0.18)
    #' p0_example <- c(0.12, 0.001, 0.08, 0.16, 0.015)
    #' sim_data <- simulate_trial_data(or_example, freq_example, 1000, p0_example)
    #' head(sim_data)
    simulate_trial_data <- function(or_vector, freq_vector, n, p0_vector = NULL, seed = 123) {
      
      # --- Input Validation ---
      # Set seed for reproducibility
      set.seed(seed)
      
      # Ensure OR vector and frequency vector are of the same length
      if(length(or_vector) != length(freq_vector)) {
        stop("OR vector and frequency vector must have the same length.")
      }

      # Get group labels from names, ensuring all vectors are named for reliable matching
      group_labels <- names(or_vector)
      if(is.null(group_labels)) {
          group_labels <- LETTERS[1:length(or_vector)]
          names(or_vector) <- group_labels
      }
      if(is.null(names(freq_vector))) names(freq_vector) <- group_labels
      if(!is.null(p0_vector) && is.null(names(p0_vector))) names(p0_vector) <- group_labels
      
      # Normalize frequency vector to sum to 1
      freq_normalized <- freq_vector / sum(freq_vector)
      
      # If p0_vector is not provided, set to default 0.3 for all groups
      if(is.null(p0_vector)) {
        p0_vector <- setNames(rep(0.3, length(or_vector)), group_labels)
      } else {
        if(length(p0_vector) != length(or_vector)) {
          stop("p0_vector must be the same length as or_vector and freq_vector.")
        }
      }
      
      # --- Data Simulation ---
      # Assign true groups to participants based on frequencies
      data <- tibble::tibble(
        id = 1:n,
        group = sample( # Sample subgroup labels based on provided frequencies
          group_labels,
          size = n,
          replace = TRUE,
          prob = freq_normalized
        )
      ) |>
        # Assign treatment randomly (0=control, 1=treatment) with equal probability
        dplyr::mutate(treatment = rbinom(n, 1, 0.5)) |>
        # Map the true OR and baseline p0 to each participant based on their group
        dplyr::mutate(
          or = or_vector[group], # Match by name for robustness
          p0 = p0_vector[group]  # Match by name for robustness
        ) |>
        # Calculate individual probability of the outcome ('success')
        # based on baseline risk (p0), treatment assignment, and subgroup OR
        dplyr::mutate(
          logit_p0 = qlogis(p0),             # Baseline risk on logit scale
          log_or = log(or),                  # Treatment effect on log scale
          logit_p1 = logit_p0 + log_or * treatment, # Risk on logit scale under assigned treatment
          prob_success = plogis(logit_p1)    # Convert back to probability scale
        ) |>
        # Simulate the binary outcome based on the calculated probability
        dplyr::mutate(
          success = rbinom(n, 1, prob_success)
        ) |>
        # Select and return relevant columns
        dplyr::select(id, group, treatment, success)
      
      return(data)
    }
    
    # estimate_effect function removed as estimate_effect_misclassify(accuracy=1) covers it.
    
    # plot_effects function removed as plotting is handled within the QMD file.
    
    #' Replicate Simulations with Misclassification
    #'
    #' Runs the simulation (`simulate_trial_data`) and estimation 
    #' (`estimate_effect_misclassify`) multiple times (`n_reps`) for a given 
    #' set of parameters, collecting results.
    #'
    #' @param or_vector Numeric vector of true odds ratios per subgroup.
    #' @param freq_vector Numeric vector of true frequencies per subgroup.
    #' @param n Integer, total sample size per simulation run.
    #' @param p0_vector Numeric vector of baseline event probabilities (control group) per subgroup.
    #' @param accuracy Numeric, the probability of correct subgroup classification (0 to 1).
    #' @param seed Integer, base seed for reproducibility across repetitions.
    #' @param n_reps Integer, the number of simulation repetitions to run.
    #'
    #' @return A tibble containing the combined results from all repetitions, 
    #'   with an added column 'k' indicating the repetition number.
    #' @export
    replicate_sims <- function(or_vector, freq_vector, n, p0_vector = NULL, accuracy = 1, seed = 123, n_reps = 100) {
      # Set base seed for the overall replication process
      set.seed(seed)
      
      # Use map_dfr to iterate n_reps times, running one simulation + estimation per iteration
      # Seeds are adjusted within the loop for reproducibility of each individual run if needed
      purrr::map_dfr(1:n_reps, function(k) {
        # Simulate data for this repetition (use adjusted seed)
        sim_data <- simulate_trial_data(or_vector, freq_vector, n, p0_vector, seed = seed + k) 
        
        # Misclassify the data first
        misclassified_data <- misclassify_group(sim_data, accuracy = accuracy, freq_vector = freq_vector, seed = seed + k + n_reps)
        
        # Estimate effects with potential misclassification
        estimate_effect_misclassify(misclassified_data, or_vector, accuracy = accuracy) |>
          dplyr::mutate(k = k) # Add repetition number
      })
    }
    
    #' Estimate Treatment Effects with Misclassification
    #'
    #' Simulates subgroup misclassification based on a given accuracy, then 
    #' estimates treatment effects (log odds ratios) via logistic regression 
    #' within each *assigned* (potentially incorrect) subgroup and overall.
    #' Calculates MSE against the true ORs.
    #'
    #' @param sim_data A tibble generated by `simulate_trial_data`, containing true subgroup labels.
    #' @param or_vector The true odds ratio vector used to generate `sim_data`.
    #' @param accuracy Numeric, the probability of correct subgroup classification (0 to 1).
    #' @param seed Integer, seed for random number generation for reproducibility of the misclassification step.
    #'
    #' @return A tibble summarizing results per assigned group: beta, se, pval, n, true_or, true_beta, or, or_lower, or_upper, mse, accuracy.
    #' @export
    estimate_effect_misclassify <- function(sim_data, or_vector, accuracy = 1) {
      
      # --- Input Validation & Setup ---
      # Compute true marginal overall log-OR for pooled data
      true_overall_beta <- glm(success ~ treatment, data = sim_data, family = binomial()) %>%
        broom::tidy() %>%
        dplyr::filter(term == "treatment") %>%
        dplyr::pull(estimate)
  
      # Ensure required packages are loaded
      pacman::p_load(tidyverse, broom)
      
      # Ensure data has required columns
      required_cols <- c("group", "treatment", "success", "assigned_group")
      if(!all(required_cols %in% names(sim_data))) {
        stop("Data must contain columns: group, treatment, success, and assigned_group")
      }
      
      # Get sorted unique group labels
      group_labels <- sort(unique(sim_data$group))
      
      # --- Estimate effects based on ASSIGNED group ---
      results <- sim_data |>
        dplyr::bind_rows(sim_data |> dplyr::mutate(assigned_group = "Overall")) |>
        dplyr::group_by(assigned_group) |>
        tidyr::nest() |>
        dplyr::mutate(
          model_results = purrr::map(data, fit_glm_safe), # Use robust GLM fitting
          n = purrr::map_int(data, nrow),
          true_or = dplyr::if_else(
            assigned_group == "Overall",
            exp(true_overall_beta),
            or_vector[match(assigned_group, group_labels)]
          ),
          true_beta = log(true_or)
        ) |>
        tidyr::unnest(model_results) |>
        dplyr::mutate(
          or = exp(beta),
          or_lower = exp(beta - 1.96 * se),
          or_upper = exp(beta + 1.96 * se),
          mse = (beta - true_beta)^2,
          accuracy = accuracy
        ) |>
        dplyr::rename(group = assigned_group) |>
        dplyr::select(group, beta, se, pval, n, true_or, true_beta, or, or_lower, or_upper, mse, accuracy) |>
        dplyr::arrange(match(group, c(group_labels, "Overall")))
  
      return(results)
}

# =============================================================================
# AIM 3 - SENSITIVITY/SPECIFICITY BASED FUNCTIONS
# =============================================================================

#' Misclassify group based on sensitivity and specificity for a target group
misclassify_group_sens_spec <- function(sim_data, target_group, sensitivity, specificity, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  
  sim_data |>
    dplyr::mutate(
      is_target = (group == target_group),
      test_positive = dplyr::case_when(
        is_target  ~ rbinom(dplyr::n(), 1, sensitivity),
        !is_target ~ rbinom(dplyr::n(), 1, 1 - specificity)
      ),
      assigned_group = ifelse(test_positive == 1, target_group, "Not_Target")
    )
}

#' Run a single enrichment scenario using sens/spec
run_enrichment_scenario_sens_spec <- function(n_screened, target_group, sensitivity, specificity, or_vector, p0_vector, freq_vector, n_reps, base_seed) {
    results <- purrr::map(1:n_reps, ~{
      pool_data <- simulate_trial_data(or_vector, freq_vector, n = n_screened, p0_vector = p0_vector, seed = base_seed + .x)
      tested_pool <- misclassify_group_sens_spec(pool_data, target_group, sensitivity, specificity, seed = base_seed + .x + n_reps)
      enrolled_cohort <- tested_pool %>% dplyr::filter(assigned_group == target_group)
      
      if (nrow(enrolled_cohort) < 20) {
        return(list(pval = 1, n_enrolled = nrow(enrolled_cohort), beta_hat = NA))
      }
      model <- fit_glm_safe(enrolled_cohort)
      return(list(pval = model$pval, n_enrolled = nrow(enrolled_cohort), beta_hat = model$beta))
    })

    p_values <- sapply(results, `[[`, "pval")
    n_values <- sapply(results, `[[`, "n_enrolled")
    beta_hats <- sapply(results, `[[`, "beta_hat")
    power <- mean(p_values < 0.05, na.rm = TRUE)
    
    return(list(power = power, mean_n_enrolled = mean(n_values, na.rm = TRUE), beta_hats = beta_hats, p_values = p_values))
}


#' Find NNS for a scenario using adaptive search (sens/spec version)
find_nns_for_scenario_sens_spec <- function(target_group, sensitivity, specificity, seed,
                                  or_vector, p0_vector, freq_vector, n_reps_per_calc,
                                  target_power = 0.8,
                                  initial_nns = 4000,
                                  max_nns = 1000000,
                                  max_iter = 8) {

  message(sprintf("\n[%s/%s] Searching: Group '%s', Sens: %s, Spec: %s, Scenario: %s",
                  get_global_counter(), get_total_scenarios(),
                  target_group, sensitivity, specificity, get_current_scenario()))
  
  power_calc_wrapper <- function(nns, sd) {
    res <- run_enrichment_scenario_sens_spec(
      n_screened = nns, target_group = target_group, 
      sensitivity = sensitivity, specificity = specificity, 
      or_vector = or_vector, p0_vector = p0_vector, freq_vector = freq_vector,
      n_reps = n_reps_per_calc, base_seed = sd
    )
    return(res$power)
  }

  message("  Phase 1: Finding power bracket...")
  low_nns <- high_nns <- initial_nns
  low_power <- high_power <- power_calc_wrapper(initial_nns, seed)
  message(sprintf("    NNS=%s -> Power=%.0f%%", format(initial_nns, big.mark=","), low_power * 100))
  
  if (low_power >= target_power) {
    while(low_power >= target_power && low_nns > 1) {
      high_nns <- low_nns; high_power <- low_power
      low_nns <- floor(low_nns / 2)
      if (low_nns == 0) low_nns <- 1
      low_power <- power_calc_wrapper(low_nns, seed)
      message(sprintf("    NNS=%s -> Power=%.0f%%", format(low_nns, big.mark=","), low_power * 100))
      if (low_nns == 1 && low_power >= target_power) break
    }
  } else {
    while(high_power < target_power && high_nns <= max_nns) {
      low_nns <- high_nns; low_power <- high_power
      high_nns <- high_nns * 2
      high_power <- power_calc_wrapper(high_nns, seed)
      message(sprintf("    NNS=%s -> Power=%.0f%%", format(high_nns, big.mark=","), high_power * 100))
    }
  }
  
  if (high_power < target_power) {
    message(paste0("  Power goal not reached below ", format(max_nns, big.mark=","), " NNS."))
    return(tibble(target_group = target_group, sensitivity = sensitivity, specificity = specificity, nns_needed = NA_integer_, nnr_corresponding = NA_real_))
  }
  
  message("  Phase 2: Refining NNS with binary search...")
  for(i in 1:max_iter) {
    mid_nns <- floor(mean(c(low_nns, high_nns)))
    if (mid_nns %in% c(low_nns, high_nns)) break
    mid_power <- power_calc_wrapper(mid_nns, seed)
    message(sprintf("    Testing NNS=%s -> Power=%.0f%%", format(mid_nns, big.mark=","), mid_power * 100))
    if (mid_power < target_power) low_nns <- mid_nns else high_nns <- mid_nns
  }
  
  final_nns <- high_nns
  final_scenario_result <- run_enrichment_scenario_sens_spec(
      n_screened = final_nns, target_group = target_group, sensitivity = sensitivity, specificity = specificity,
      or_vector = or_vector, p0_vector = p0_vector, freq_vector = freq_vector,
      n_reps = n_reps_per_calc * 2, base_seed = seed # Use same seed as power search
  )
  
  # Calculate bias and MSE from the final, high-replication run
  beta_target <- log(or_vector[target_group])
  beta_hats <- final_scenario_result$beta_hats
  final_p_values <- final_scenario_result$p_values
  
  bias <- mean(beta_hats - beta_target, na.rm = TRUE)
  mse <- mean((beta_hats - beta_target)^2, na.rm = TRUE)
  
  # Calculate additional metrics
  mean_beta_hat <- mean(beta_hats, na.rm = TRUE)
  sd_beta_hat <- sd(beta_hats, na.rm = TRUE)
  rmse <- sqrt(mse)
  
  # Calculate proportion of trials with wrong direction (sign flipping)
  wrong_direction <- mean(sign(beta_hats) != sign(beta_target), na.rm = TRUE)
  
  # Calculate proportion of trials with significant results in wrong direction
  significant_wrong <- mean(final_p_values < 0.05 & sign(beta_hats) != sign(beta_target), na.rm = TRUE)

  return(tibble(
    target_group = target_group, sensitivity = sensitivity, specificity = specificity,
    nns_needed = final_nns, nnr_corresponding = final_scenario_result$mean_n_enrolled,
    true_beta = beta_target, mean_beta_hat = mean_beta_hat, bias = bias, mse = mse, rmse = rmse,
    sd_beta_hat = sd_beta_hat, wrong_direction = wrong_direction, significant_wrong = significant_wrong
  ))
}

# Global counters and data for Aim 3
.global_counter <- 0
.total_scenarios <- 0
.current_scenario <- "N/A"
reset_global_counter <- function(total) { .global_counter <<- 0; .total_scenarios <<- total }
increment_global_counter <- function() { .global_counter <<- .global_counter + 1 }
set_current_scenario <- function(name) { .current_scenario <<- name }
get_global_counter <- function() .global_counter
get_total_scenarios <- function() .total_scenarios
get_current_scenario <- function() .current_scenario
or_mortality <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
freq_arrest <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)
p0_arrest_raw <- c(A = 7/60, B = 0/52, C = 11/138, D = 11/69, E = 1/69)
n0_B <- 52; p0_B_corrected <- 0.5 / (n0_B + 0.5); p0_arrest_adjusted <- p0_arrest_raw; p0_arrest_adjusted["B"] <- p0_B_corrected
