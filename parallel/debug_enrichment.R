# Debug script to check enrichment scenario
pacman::p_load(tidyverse, broom, furrr, scales)
source("../code/functions/functions.R")

# Test parameters
or_arrest <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
freq_arrest <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)
p0_arrest_raw <- c(A = 7/60, B = 0/52, C = 11/138, D = 11/69, E = 1/69)
p0_B_corrected <- 0.5 / (52 + 0.5)
p0_arrest_adjusted <- p0_arrest_raw
p0_arrest_adjusted["B"] <- p0_B_corrected

# Test with perfect classification
target_group <- "B"
sensitivity <- 1.0
specificity <- 1.0
n_screened <- 2375  # From the result
seed <- 4001

set.seed(seed)

# Simulate one scenario
pool_data <- simulate_trial_data(or_arrest, freq_arrest, n = n_screened, p0_vector = p0_arrest_adjusted, seed = seed)
tested_pool <- misclassify_group_sens_spec(pool_data, target_group, sensitivity, specificity, seed = seed + 200)
enrolled_cohort <- tested_pool %>% dplyr::filter(assigned_group == target_group)

print("Pool composition:")
print(table(pool_data$group))

print("Test results:")
print(table(tested_pool$assigned_group))

print("Enrolled cohort composition:")
print(table(enrolled_cohort$group))

print("Enrolled cohort size:")
print(nrow(enrolled_cohort))

# Check if there's any contamination
contamination <- enrolled_cohort %>% filter(group != target_group)
print("Contamination (non-B patients in enrolled cohort):")
print(nrow(contamination))
if(nrow(contamination) > 0) {
  print(table(contamination$group))
}

# Fit model to enrolled cohort
if(nrow(enrolled_cohort) >= 20) {
  model <- fit_glm_safe(enrolled_cohort)
  print("Model results:")
  print(model)
  
  # Check true vs estimated
  true_beta <- log(or_arrest[target_group])
  print(paste("True beta:", true_beta))
  print(paste("Estimated beta:", model$beta))
  print(paste("Bias:", model$beta - true_beta))
} 