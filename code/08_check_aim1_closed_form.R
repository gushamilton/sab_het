# Compare Aim 1 simulation vs closed-form approximation

pacman::p_load(tidyverse)

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_sim <- file.path(results_root, "tables", "aim1_power_summary.tsv")
output_cf <- file.path(results_root, "tables", "aim1_closed_form_power.tsv")
output_cmp <- file.path(results_root, "tables", "aim1_closed_form_vs_sim.tsv")

if (!file.exists(input_sim)) {
  stop("Missing input: ", input_sim)
}

source("code/common_parameters.R")
params <- get_common_params()

sim <- readr::read_tsv(input_sim, show_col_types = FALSE)

# pick OR vector based on scenario name
get_or_vector <- function(name) {
  if (name == "ARREST_raw") return(params$or_arrest)
  if (name == "ARREST_shrunk") return(params$or_arrest_shrunk_k05)
  if (name == "Conservative") return(params$or_conservative)
  stop("Unknown scenario: ", name)
}

alpha <- 0.05 / 5
z_alpha <- qnorm(1 - alpha / 2)

# closed-form power per subgroup
calc_power <- function(n_total, group, or_vec, p0_vec, freq_vec) {
  f <- freq_vec[group]
  n_sub <- n_total * f
  n0 <- n_sub / 2
  n1 <- n_sub / 2
  or <- or_vec[group]
  p0 <- p0_vec[group]
  p1 <- (or * p0) / (1 - p0 + or * p0)
  log_or <- log(or)
  var_term <- (1 / (n0 * p0 * (1 - p0))) + (1 / (n1 * p1 * (1 - p1)))
  z <- abs(log_or) / sqrt(var_term)
  # two-sided power
  power <- 1 - (pnorm(z_alpha - z) - pnorm(-z_alpha - z))
  power
}

cf <- sim %>%
  mutate(
    or_vector = map(scenario_name, get_or_vector),
    p0_vector = list(params$p0_overall),
    freq_vector = list(params$freq_arrest)
  ) %>%
  rowwise() %>%
  mutate(
    power_cf = if (group == "Overall") NA_real_ else calc_power(n_total, group, or_vector, p0_vector, freq_vector)
  ) %>%
  ungroup() %>%
  select(scenario_name, n_total, group, power_cf, scenario_set)

write_tsv(cf, output_cf)

cmp <- sim %>%
  left_join(cf, by = c("scenario_name", "n_total", "group", "scenario_set")) %>%
  mutate(diff = power - power_cf)

write_tsv(cmp, output_cmp)

summary <- cmp %>%
  filter(group != "Overall") %>%
  summarise(
    mean_abs_diff = mean(abs(diff), na.rm = TRUE),
    max_abs_diff = max(abs(diff), na.rm = TRUE)
  )

print(summary)
