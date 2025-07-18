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
p0_B_corrected     <- 0.5 / (107 + 0.5)
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