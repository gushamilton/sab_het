# Debug script to check enrichment scenario
pacman::p_load(tidyverse, broom, furrr, scales)
source("../code/functions/functions.R")

# Test parameters
or_arrest <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
# Updated to use overall prevalence (both treatment and control arms combined)
# From Swets et al. CID Supplementary Table 2:
# Group A: 60+55=115, Group B: 52+55=107, Group C: 138+138=276, Group D: 69+52=121, Group E: 69+70=139
# Total: 115+107+276+121+139 = 758
freq_arrest <- c(A = 115/758, B = 107/758, C = 276/758, D = 121/758, E = 139/758)
# Updated mortality data based on paper (overall mortality across both arms)
# From Swets et al. CID Supplementary Table 2:
# Group A: (13+12)/(60+55)=25/115=21.7%, Group B: (0+8)/(52+55)=8/107=7.5%, 
# Group C: (29+24)/(138+138)=53/276=19.2%, Group D: (11+11)/(69+52)=22/121=18.2%, 
# Group E: (3+1)/(69+70)=4/139=2.9%
p0_arrest_raw <- c(A = 25/115, B = 8/107, C = 53/276, D = 22/121, E = 4/139)
p0_B_corrected <- 0.5 / (107 + 0.5)
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