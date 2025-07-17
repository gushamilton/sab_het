# Script to identify scenario parameters
library(tidyverse)

# Define the scenarios (same as in aim3_parallel.R)
or_arrest <- c(A = 1.0, B = 18.8, C = 0.79, D = 1.4, E = 0.3)
or_conservative <- c(A = 1.0, B = 2.0, C = 0.7, D = 1.2, E = 0.8)

scenario_definitions <- tribble(
  ~scenario_name, ~or_vector,
  "ARREST",        or_arrest,
  "Conservative",  or_conservative
)

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

# Show all scenarios
print("All scenarios:")
print(all_aim3_scenarios %>% select(scenario_name, test_type, target_group, sensitivity, specificity))

# Show missing scenarios (20, 22, 37, 45)
missing_scenarios <- c(20, 22, 37, 45)
print("\nMissing scenarios:")
print(all_aim3_scenarios[missing_scenarios, ] %>% select(scenario_name, test_type, target_group, sensitivity, specificity)) 