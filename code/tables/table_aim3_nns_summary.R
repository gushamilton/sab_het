# Create Aim 3 NNS summary table (simulation vs closed form)

pacman::p_load(tidyverse)

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "supp")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_sim <- file.path(results_root, "tables", "aim3_sens_spec_summary.tsv")
input_closed <- file.path(results_root, "tables", "aim3_closed_form_summary.tsv")
output_file <- file.path(results_root, "tables", "aim3_nns_table.tsv")

if (!file.exists(input_sim)) {
  stop("Missing input: ", input_sim)
}
if (!file.exists(input_closed)) {
  stop("Missing input: ", input_closed)
}

sim <- readr::read_tsv(input_sim, show_col_types = FALSE) %>%
  transmute(
    scenario_name,
    test_type,
    target_group,
    method = "Simulation",
    NNS = nns_needed,
    NNR = nnr_corresponding,
    Bias = bias
  )

closed <- readr::read_tsv(input_closed, show_col_types = FALSE) %>%
  transmute(
    scenario_name,
    test_type,
    target_group,
    method = "Closed Form",
    NNS = nns_closed_form,
    NNR = nnr_closed_form,
    Bias = bias_closed_form
  )

combined <- bind_rows(sim, closed) %>%
  arrange(scenario_name, target_group, test_type, method)

readr::write_tsv(combined, output_file)
message("Saved: ", output_file)
