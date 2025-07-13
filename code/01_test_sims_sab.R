
pacman::p_load(tidyverse,  data.table, vroom, ggdist)
theme_set(theme_ggdist())

source("code/functions/functions.R")

# Example usage given the data on death in ARREST

# Define odds ratios and frequencies
or_example <- c(1, 18.8, 0.79, 1.4, 0.3) #odds ratio for death in each group
# or_example <- c(0.35,2.9,0.17,0.64,0.31)
freq_example <- c(60, 52, 138, 69, 69)  # Example frequencies in each group
p0_example <- c(0.3, 0.35, 0.25, 0.4, 0.2) # baseline probabilities of death
n_example <-1e3  # Sample size

# Define baseline probabilities for each group


# Generate simulated data
simulated_data <- simulate_trial_data(or_example, freq_example, n_example, p0_example, seed = 2)
simulated_data
# View the first few rows of the simulated data
head(simulated_data)

# Estimate effect and plot results
results <- estimate_effect_misclassify(simulated_data, or_example, accuracy = 0.2)
plot_effects(results)




# Generate simulated data
simulated_data <- simulate_trial_data(or_example, freq_example, n_example, p0_example, seed = 2)
simulated_data
# View the first few rows of the simulated data
head(simulated_data)

# Estimate effect and plot results
results <- estimate_effect(simulated_data, or_example, accuracy = 0.2)
plot_effects(results)

# do it at scale to calculate power

results_reps <- replicate_sims(or_example, freq_example, n_example, p0_example, n_reps = 100, seed = 2423, accuracy = 1)
results_reps  |> 
group_by(group) |> 
summarise(mean(mse), power = mean(pval < 0.05/5))


# now using a misclassified variable
results_reps <- replicate_sims(or_example, freq_example, n_example, p0_example, n_reps = 100, seed = 2423, accuracy = 0.2)
results_reps  |> 
group_by(group) |> 
summarise(mean(mse), power = mean(pval < 0.05/5))


over_misclass <- map_dfr(seq(0.5, 1, by = 0.1), ~replicate_sims(or_example, freq_example, n = 2000, p0_example, n_reps = 100, seed = 2423, accuracy = .x))



over_misclass |> 
group_by(group, accuracy) |> 
summarise(mean(mse), power = mean(pval < 0.05/5))

over_misclass |> 
group_by(group, accuracy) |> 
summarise(mean(mse), power = mean(pval < 0.05/5)) |> 
ggplot(aes(x = accuracy, y = power, color = group)) +
geom_point() +
geom_line() 


over_misclass <- map_dfr(seq(0.5, 1, by = 0.1), ~replicate_sims(or_example, freq_example, 5e3, p0_example, n_reps = 200, seed = 2423, accuracy = .x))
over_misclass |> 
group_by(group, accuracy) |> 
summarise(mean(mse), power = mean(pval < 0.05/5))

over_misclass |> 
group_by(group, accuracy) |> 
summarise(mean(mse), power = mean(pval < 0.05/5)) |> 
ggplot(aes(x = accuracy, y = power, color = group)) +
geom_point() +
geom_line() 


over_misclass |> 
group_by(group, accuracy) |> 
summarise(
  mean(mse), 
  power = mean(pval < 0.05/5),
  wrong_dir = mean(sign(beta) != sign(true_beta)),
  wrong_dir_sig = mean(sign(beta) != sign(true_beta) & pval < 0.05)
) |> 
ggplot(aes(x = accuracy, y = wrong_dir, color = group)) +
geom_point() +
geom_line() 

# Now I want to do it for 4 groups with more realistic values 
# Define OR distribution and proportions for subgroups

# Define proportions and baseline mortality for four subgroups
proportions <- c(0.4, 0.3, 0.2, 0.1)  # Example proportions
p0_distribution <- c(0.1, 0.2, 0.3, 0.4)  # Varying baseline probabilities


# Choose ORs for the first three subgroups
or1 <- 1.6
or2 <- 1.25
or3 <- 0.65

# Calculate the required OR4 to balance the overall OR to 1
log_or1 <- log(or1)
log_or2 <- log(or2)
log_or3 <- log(or3)

# Calculate log OR4
log_or4 <- -(log_or1 * proportions[1] + log_or2 * proportions[2] + log_or3 * proportions[3]) / proportions[4]
or4 <- exp(log_or4)

or4 <- 0.35
# Define the complete OR distribution and calculate Z scores
or_distribution <- c(or1, or2, or3, or4)


# Calculate weighted sum of log ORs
weighted_log_or <- sum(log(or_distribution) * proportions)

# Calculate the overall OR
overall_or <- exp(weighted_log_or)

# Print the results
print(data.frame(
  Subgroup = paste0("Group ", 1:4),
  OR = or_distribution,
  Proportion = proportions,
  p0 = p0_distribution
))
print(overall_or)  # Should output 1


# Generate simulated data with adjusted ORs and p0
simulated_data <- simulate_trial_data(
  or_vector = or_distribution,
  freq_vector = proportions,
  n = 1e6,
  p0_vector = p0_distribution,
  seed = 42
)

# Estimate effect and plot results
results <- estimate_effect(simulated_data, or_distribution)
plot_effects(results)

# Further simulations as needed...


tidy(glm(treatment ~ success, data = simulated_data, family = binomial))
# View the first few rows of the simulated data
head(simulated_data)

# Estimate effect and plot results


# Run simulations across different accuracy levels
over_misclass_new <- map_dfr(
  seq(0.5, 1, by = 0.05), 
  ~replicate_sims(
    or_vector = or_distribution,
    freq_vector = proportions,
    n = 5000,
    p0_vector = p0_distribution,
    n_reps = 300,
    seed = 2423,
    accuracy = .x
  )
)

# Analyze results
over_misclass_new |> 
  group_by(group, accuracy) |> 
  summarise(
    mse = mean(mse),
    power = mean(pval < 0.05/5),
    wrong_dir = mean(sign(beta) != sign(true_beta)),
    wrong_dir_sig = mean(sign(beta) != sign(true_beta) & pval < 0.05),
    .groups = "drop"
  )

# Plot power across accuracy levels
over_misclass_new |> 
  group_by(group, accuracy) |> 
  summarise(
    power = mean(pval < 0.05/5),
    .groups = "drop"
  ) |> 
  ggplot(aes(x = accuracy, y = power, color = group)) +
  geom_point() +
  geom_line() +
  labs(
    x = "Classification Accuracy",
    y = "Power",
    title = "Statistical Power vs Classification Accuracy"
  )


over_misclass_new |> 
  group_by(group, accuracy) |> 
  summarise(
        wrong_dir = mean(sign(beta) != sign(true_beta)),
    .groups = "drop"
  ) |> 
  ggplot(aes(x = accuracy, y = wrong_dir, color = group)) +
  geom_point() +
  geom_line() +
  labs(
    x = "Classification Accuracy",
    y = "Wrong direction",
    title = "Wrong direction vs Classification Accuracy"
  )
