# Create Aim 3 supplementary table S2 (simulation summary by test type)

pacman::p_load(tidyverse, gt, scales)

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_file <- file.path(results_root, "tables", "aim3_sens_spec_summary.tsv")
output_tsv <- file.path(results_root, "tables", "aim3_table_S2.tsv")
output_html <- file.path(results_root, "tables", "aim3_table_S2.html")

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

aim3 <- readr::read_tsv(input_file, show_col_types = FALSE) %>%
  filter(scenario_name == "ARREST_raw", target_group %in% c("B", "C", "D", "E")) %>%
  mutate(
    NNS = ifelse(is.na(nns_needed), ">1,000,000", scales::comma(round(nns_needed, 0))),
    NNR = ifelse(is.na(nnr_corresponding), "-", scales::comma(round(nnr_corresponding, 0))),
    Bias = ifelse(is.na(bias), "-", sprintf("%.3f", bias)),
    MSE = ifelse(is.na(mse), "-", sprintf("%.3f", mse)),
    Wrong_dir = ifelse(is.na(wrong_direction), "-", sprintf("%.1f", wrong_direction))
  ) %>%
  select(target_group, test_type, NNS, NNR, Bias, MSE, Wrong_dir) %>%
  arrange(target_group, test_type)

readr::write_tsv(aim3, output_tsv)

gt_tbl <- aim3 %>%
  gt(groupname_col = "target_group") %>%
  cols_label(
    test_type = "Test Type",
    NNS = "NNS",
    NNR = "NNR",
    Bias = "Bias",
    MSE = "MSE",
    Wrong_dir = "Wrong dir (%)"
  ) %>%
  tab_header(title = "Table S2: Aim 3 simulation summary by test type") %>%
  tab_options(table.width = "100%")

gtsave(gt_tbl, output_html)
message("Saved: ", output_tsv)
message("Saved: ", output_html)
