.# Detailed debug script to find the source of bias
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
n_screened <- 2375
seed <- 4001
n_reps_per_calc <- 20

set.seed(seed)

print("=== TRACING find_nns_for_scenario_sens_spec FUNCTION ===")

# Step 1: Run the final scenario with the same parameters as the function
print("Step 1: Running final scenario with n_reps_per_calc * 2 = 40 replicates...")

final_scenario_result <- run_enrichment_scenario_sens_spec(
    n_screened = n_screened, 
    target_group = target_group, 
    sensitivity = sensitivity, 
    specificity = specificity,
    or_vector = or_arrest, 
    p0_vector = p0_arrest_adjusted, 
    freq_vector = freq_arrest,
    n_reps = n_reps_per_calc * 2, 
    base_seed = seed
)

print("Final scenario results:")
print(paste("Power:", final_scenario_result$power))
print(paste("Mean N enrolled:", final_scenario_result$mean_n_enrolled))
print("Beta hats (first 10):")
print(head(final_scenario_result$beta_hats, 10))

# Step 2: Calculate metrics exactly as the function does
print("\nStep 2: Calculating metrics...")
beta_target <- log(or_arrest[target_group])
beta_hats <- final_scenario_result$beta_hats
final_p_values <- final_scenario_result$p_values

print(paste("True beta (log OR):", beta_target))
print(paste("Number of valid beta_hats:", sum(!is.na(beta_hats))))
print(paste("Number of NA beta_hats:", sum(is.na(beta_hats))))

# Check for extreme values
print("Beta hat summary statistics:")
print(summary(beta_hats))

# Calculate bias and MSE
bias <- mean(beta_hats - beta_target, na.rm = TRUE)
mse <- mean((beta_hats - beta_target)^2, na.rm = TRUE)
mean_beta_hat <- mean(beta_hats, na.rm = TRUE)
sd_beta_hat <- sd(beta_hats, na.rm = TRUE)
rmse <- sqrt(mse)

print(paste("Mean beta hat:", mean_beta_hat))
print(paste("Bias:", bias))
print(paste("MSE:", mse))
print(paste("RMSE:", rmse))
print(paste("SD beta hat:", sd_beta_hat))

# Check for wrong direction
wrong_direction <- mean(sign(beta_hats) != sign(beta_target), na.rm = TRUE)
significant_wrong <- mean(final_p_values < 0.05 & sign(beta_hats) != sign(beta_target), na.rm = TRUE)

print(paste("Wrong direction proportion:", wrong_direction))
print(paste("Significant wrong proportion:", significant_wrong))

# Step 3: Look at individual replicates that might be causing issues
print("\nStep 3: Examining individual replicates...")
for(i in 1:min(10, length(beta_hats))) {
  if(!is.na(beta_hats[i])) {
    print(paste("Replicate", i, "beta_hat:", beta_hats[i], "bias:", beta_hats[i] - beta_target))
  }
}

# Step 4: Check if there are any extreme outliers
print("\nStep 4: Checking for extreme outliers...")
beta_hats_clean <- beta_hats[!is.na(beta_hats)]
if(length(beta_hats_clean) > 0) {
  q99 <- quantile(beta_hats_clean, 0.99)
  q01 <- quantile(beta_hats_clean, 0.01)
  print(paste("99th percentile:", q99))
  print(paste("1st percentile:", q01))
  
  extreme_high <- sum(beta_hats_clean > q99 + 3 * (q99 - q01))
  extreme_low <- sum(beta_hats_clean < q01 - 3 * (q99 - q01))
  print(paste("Extreme high outliers:", extreme_high))
  print(paste("Extreme low outliers:", extreme_low))
}

# Step 5: Examine individual replicates in detail
print("\nStep 5: Examining individual replicates in detail...")

# Look at a few replicates with extreme values
extreme_indices <- which(beta_hats > 10)
normal_indices <- which(beta_hats <= 10 & beta_hats > 0)

print(paste("Number of extreme replicates (beta > 10):", length(extreme_indices)))
print(paste("Number of normal replicates (0 < beta <= 10):", length(normal_indices)))

if(length(extreme_indices) > 0) {
  print("Examining an extreme replicate...")
  rep_idx <- extreme_indices[1]
  
  # Re-run this specific replicate
  set.seed(seed + rep_idx)
  pool_data <- simulate_trial_data(or_arrest, freq_arrest, n = n_screened, p0_vector = p0_arrest_adjusted, seed = seed + rep_idx)
  tested_pool <- misclassify_group_sens_spec(pool_data, target_group, sensitivity, specificity, seed = seed + rep_idx + n_reps_per_calc * 2)
  enrolled_cohort <- tested_pool %>% dplyr::filter(assigned_group == target_group)
  
  print(paste("Extreme replicate", rep_idx, "enrolled cohort size:", nrow(enrolled_cohort)))
  print("Enrolled cohort composition:")
  print(table(enrolled_cohort$group))
  
  if(nrow(enrolled_cohort) >= 20) {
    model <- fit_glm_safe(enrolled_cohort)
    print(paste("Model beta:", model$beta))
    print(paste("Model p-value:", model$pval))
    
    # Check treatment x outcome table
    print("Treatment x Outcome table:")
    print(table(enrolled_cohort$treatment, enrolled_cohort$success))
  }
}

if(length(normal_indices) > 0) {
  print("Examining a normal replicate...")
  rep_idx <- normal_indices[1]
  
  # Re-run this specific replicate
  set.seed(seed + rep_idx)
  pool_data <- simulate_trial_data(or_arrest, freq_arrest, n = n_screened, p0_vector = p0_arrest_adjusted, seed = seed + rep_idx)
  tested_pool <- misclassify_group_sens_spec(pool_data, target_group, sensitivity, specificity, seed = seed + rep_idx + n_reps_per_calc * 2)
  enrolled_cohort <- tested_pool %>% dplyr::filter(assigned_group == target_group)
  
  print(paste("Normal replicate", rep_idx, "enrolled cohort size:", nrow(enrolled_cohort)))
  print("Enrolled cohort composition:")
  print(table(enrolled_cohort$group))
  
  if(nrow(enrolled_cohort) >= 20) {
    model <- fit_glm_safe(enrolled_cohort)
    print(paste("Model beta:", model$beta))
    print(paste("Model p-value:", model$pval))
    
    # Check treatment x outcome table
    print("Treatment x Outcome table:")
    print(table(enrolled_cohort$treatment, enrolled_cohort$success))
  }
} 