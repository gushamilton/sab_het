# Debug script to check vector formats
pacman::p_load(tidyverse, broom, furrr, scales)
source("../code/functions/functions.R")

# Define parameters
or_arrest <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
or_conservative <- c(A = 1.0, B = 2.0, C = 0.7, D = 1.2, E = 0.8)

scenario_definitions <- tribble(
  ~scenario_name, ~or_vector,
  "ARREST",        or_arrest,
  "Conservative",  or_conservative
)

freq_arrest        <- c(A = 60/388, B = 52/388, C = 138/388, D = 69/388, E = 69/388)
p0_arrest_raw      <- c(A = 7/60, B = 0/52, C = 11/138, D = 11/69, E = 1/69)
p0_B_corrected     <- 0.5 / (52 + 0.5)
p0_arrest_adjusted <- p0_arrest_raw
p0_arrest_adjusted["B"] <- p0_B_corrected

sens_spec_scenarios <- tribble(
  ~test_type,                  ~sensitivity, ~specificity,
  "Perfect (100%)",            1.00,         1.00,
  "Near-Perfect",              0.99,         0.99,
  "High Sens/High Spec",       0.95,         0.95,
  "High Sens/Low Spec",        0.95,         0.70,
  "Low Sens/High Spec",        0.70,         0.95,
  "Balanced/Moderate",         0.80,         0.80
)
target_groups_aim3 <- c("B", "C", "D", "E")

all_aim3_scenarios <- expand_grid(
  scenario_definitions,
  sens_spec_scenarios,
  target_group = target_groups_aim3
) %>% mutate(seed = 1:n() + 4000)

# Check scenario 1
scenario_row <- all_aim3_scenarios[1, ]
print("Scenario row:")
print(scenario_row)

print("OR vector type and length:")
print(class(scenario_row$or_vector))
print(length(scenario_row$or_vector))
print(scenario_row$or_vector)

print("Freq vector type and length:")
print(class(freq_arrest))
print(length(freq_arrest))
print(freq_arrest)

print("P0 vector type and length:")
print(class(p0_arrest_adjusted))
print(length(p0_arrest_adjusted))
print(p0_arrest_adjusted)

# Try unlisting
or_vector <- unlist(scenario_row$or_vector)
print("Unlisted OR vector:")
print(class(or_vector))
print(length(or_vector))
print(or_vector)

# Test simulate_trial_data
print("Testing simulate_trial_data...")
test_data <- simulate_trial_data(or_vector, freq_arrest, n = 100, p0_vector = p0_arrest_adjusted, seed = 123)
print("Success! Data generated:")
print(head(test_data)) 