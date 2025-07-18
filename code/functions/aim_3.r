# --- Aim 3: Enrichment Trial Simulation (Direct Screening - map version) ---

# --- Libraries ---
library(tidyverse)
library(broom)
# library(future) # Uncomment for parallel processing
# library(furrr)  # Uncomment for parallel processing

# --- Simulation Parameters ---
TARGET_SUBGROUP <- "C" # Choose "B", "C", or "E"
TEST_ACCURACY <- 0.9  # Choose accuracy (e.g., 1.0, 0.95, 0.90, 0.80)
N_REPS <- 20          # Number of replicates (INCREASE TO 1000+ FOR REAL RESULTS)
# Define a range of SCREENING sizes
SCREENED_SIZES <- c(1000, 5000, 10000, 20000, 40000, 60000, 80000, 100000, 150000, 200000)
ALPHA_AIM3 <- 0.05
BASE_SEED <- 3000

# --- Population & Effect Parameters (Copied from QMD) ---
or_mortality       <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
# Updated to use overall prevalence (both treatment and control arms combined)
# From Swets et al. CID Supplementary Table 2:
# Group A: 60+55=115, Group B: 52+55=107, Group C: 138+138=276, Group D: 69+52=121, Group E: 69+70=139
# Total: 115+107+276+121+139 = 758
freq_arrest        <- c(A = 115/758, B = 107/758, C = 276/758, D = 121/758, E = 139/758)
# Updated mortality data based on paper (overall mortality across both arms)
# From Swets et al. CID Supplementary Table 2:
# Group A: (13+12)/(60+55)=25/115=21.7%, Group B: (0+8)/(52+55)=8/107=7.5%, 
# Group C: (29+24)/(138+138)=53/276=19.2%, Group D: (11+11)/(69+52)=22/121=18.2%, 
# Group E: (3+1)/(69+70)=4/139=2.9%
p0_arrest_raw      <- c(A = 25/115, B = 8/107, C = 53/276, D = 22/121, E = 4/139)
n0_B <- 107
p0_B_corrected <- 0.5 / (n0_B + 0.5)
p0_arrest_adjusted <- p0_arrest_raw
p0_arrest_adjusted["B"] <- p0_B_corrected

# --- Helper Functions (Based on user provided code) ---

# Function to simulate trial data
simulate_trial_data <- function(or_vector, freq_vector, n, p0_vector = NULL, seed = 123) {
  set.seed(seed)
  if(length(or_vector) != length(freq_vector)) {
    stop("OR vector and frequency vector must have the same length.")
  }
  # Ensure vectors have names for reliable matching
  group_labels <- names(or_vector)
  if(is.null(group_labels)) group_labels <- LETTERS[1:length(or_vector)]
  names(or_vector) <- group_labels
  names(freq_vector) <- group_labels
  if(!is.null(p0_vector) && is.null(names(p0_vector))) names(p0_vector) <- group_labels

  freq_normalized <- freq_vector / sum(freq_vector)

  if(is.null(p0_vector)) {
    p0_vector <- setNames(rep(0.3, length(or_vector)), group_labels)
  } else {
     if(length(p0_vector) != length(or_vector)) {
        stop("p0_vector must be the same length as or_vector and freq_vector.")
     }
  }

  data <- tibble::tibble(
    id = 1:n,
    group = sample(
      group_labels,
      size = n,
      replace = TRUE,
      prob = freq_normalized
    )
  ) |>
    dplyr::mutate(treatment = rbinom(n, 1, 0.5)) |>
    dplyr::mutate(
      or = or_vector[group], # Match by name
      p0 = p0_vector[group]  # Match by name
    ) |>
    dplyr::mutate(
      logit_p0 = qlogis(p0),
      log_or = log(or),
      logit_p1 = logit_p0 + log_or * treatment,
      prob_success = plogis(logit_p1)
    ) |>
    dplyr::mutate(
      success = rbinom(n, 1, prob_success)
    ) |>
    dplyr::select(id, group, treatment, success)
  return(data)
}

# Function to misclassify group (Based on logic within user's estimate_effect_misclassify)
misclassify_group <- function(sim_data, accuracy, freq_vector, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  n <- nrow(sim_data)
  group_labels <- names(freq_vector)
  if(is.null(group_labels)) group_labels <- LETTERS[1:length(freq_vector)]
  names(freq_vector) <- group_labels
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


# Function to fit GLM safely
fit_glm_safe <- function(data) {
  if (nrow(data) < 10 || length(unique(data$treatment)) < 2 || length(unique(data$success)) < 2) {
    return(tibble(beta = NA_real_, se = NA_real_, pval = NA_real_))
  }
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
}

# --- Simulation Loop ---

# Function to run one replicate
run_single_enrichment_rep <- function(rep_id, n_screened, target_group, accuracy,
                                      or_vector, p0_vector, freq_vector,
                                      base_seed) {
  seed <- base_seed + rep_id
  # 1. Simulate screened cohort
  screened_data <- simulate_trial_data(
    or_vector = or_vector,
    freq_vector = freq_vector,
    n = n_screened,
    p0_vector = p0_vector,
    seed = seed
  )
  # 2. Apply test (misclassify)
  tested_data <- misclassify_group(
      screened_data,
      accuracy = accuracy,
      freq_vector = freq_vector, # Pass named freq_vector
      seed = seed + N_REPS # Use different seed offset
      )

  # 3. Select cohort for randomization ("test positives")
  randomized_cohort <- tested_data %>% filter(assigned_group == target_group)
  n_randomized_actual <- nrow(randomized_cohort)

  # 4. Analyze if cohort is large enough
  if(n_randomized_actual >= 10) {
     analysis_results <- fit_glm_safe(randomized_cohort)
  } else {
     analysis_results <- tibble(beta = NA_real_, se = NA_real_, pval = NA_real_)
  }

  # Return results
  tibble(
      n_randomized_actual = n_randomized_actual
      ) %>% bind_cols(analysis_results)
}

# --- Run Simulation using map_dfr ---
message(paste("Running Aim 3 Simulation for Target:", TARGET_SUBGROUP, "Accuracy:", TEST_ACCURACY))

# Use map_dfr to iterate over SCREENED_SIZES
results_aim3_standalone <- map_dfr(SCREENED_SIZES, function(n_scr) {
    message(paste("  N_Screened =", n_scr))
    # Inner map_dfr for replicates
    map_dfr(1:N_REPS, ~run_single_enrichment_rep(
        rep_id = .x,
        n_screened = n_scr,
        target_group = TARGET_SUBGROUP,
        accuracy = TEST_ACCURACY,
        or_vector = or_mortality,
        p0_vector = p0_arrest_adjusted,
        freq_vector = freq_arrest,
        base_seed = BASE_SEED + which(SCREENED_SIZES == n_scr) * N_REPS # Adjust seed per screened size
        ), .id = "rep") %>%
    mutate(n_screened = n_scr) # Add n_screened column after replicates are done
}, .id = NULL) # No need for outer .id


# --- Process Results ---
summary_aim3_standalone <- results_aim3_standalone %>%
  mutate(pval = as.numeric(pval)) %>%
  group_by(n_screened) %>%
  summarise(
    power = mean(pval < ALPHA_AIM3, na.rm = TRUE),
    mean_n_randomized = mean(n_randomized_actual, na.rm = TRUE),
    sd_n_randomized = sd(n_randomized_actual, na.rm = TRUE),
    n_reps_analyzed = sum(!is.na(pval)),
    .groups = "drop"
  )

# Find NNS required for 80% power
nns_enrich_standalone <- summary_aim3_standalone %>%
  arrange(n_screened) %>%
  filter(power >= 0.8) %>%
  slice(1) # First n_screened where power >= 0.8

# Get corresponding NNR
nnr_enrich_standalone <- if (nrow(nns_enrich_standalone) > 0) {
                            round(nns_enrich_standalone$mean_n_randomized)
                         } else {
                            NA # Indicate NNR not determined if power not reached
                         }

# --- Output Results ---
print(paste("Target Subgroup:", TARGET_SUBGROUP))
print(paste("Test Accuracy:", TEST_ACCURACY))
print("Summary Results by N Screened:")
print(summary_aim3_standalone)

if (nrow(nns_enrich_standalone) > 0) {
  nns_val <- scales::comma(nns_enrich_standalone$n_screened)
  nnr_val <- scales::comma(nnr_enrich_standalone)
  print(paste("Estimated NNS to achieve ~80% power:", nns_val))
  print(paste("Estimated Average NNR at that NNS:", nnr_val))
} else {
  print(paste("80% power not achieved within the tested N Screened range up to", scales::comma(max(SCREENED_SIZES))))
}

 Optional: Plot Power vs N Screened
 ggplot(summary_aim3_standalone, aes(x = n_screened, y = power)) +
   geom_line() + geom_point() +
   geom_hline(yintercept = 0.8, linetype = "dashed") +
   scale_x_continuous(labels = scales::comma) +
   scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
   labs(title = paste("Aim 3 Power Curve: Target =", TARGET_SUBGROUP, ", Accuracy =", TEST_ACCURACY),
        x = "Number Needed to Screen (NNS)",
        y = "Power") +
   theme_minimal()

pacman::p_load(furrr)

