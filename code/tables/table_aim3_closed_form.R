# Create Aim 3 closed-form table (NNS/NNR only)

pacman::p_load(tidyverse, gt)

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_file <- file.path(results_root, "tables", "aim3_closed_form_summary.tsv")
output_tsv <- file.path(results_root, "tables", "aim3_closed_form_table.tsv")
output_html <- file.path(results_root, "tables", "aim3_closed_form_table.html")

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

closed <- readr::read_tsv(input_file, show_col_types = FALSE) %>%
  filter(target_group %in% c("B", "C", "D", "E")) %>%
  mutate(
    NNS = scales::comma(nns_closed_form),
    NNR = scales::comma(nnr_closed_form)
  ) %>%
  select(test_type, target_group, NNS, NNR) %>%
  arrange(target_group, test_type)

readr::write_tsv(closed, output_tsv)

gt_tbl <- closed %>%
  gt(groupname_col = "target_group") %>%
  cols_label(
    test_type = "Test Type",
    NNS = "NNS (Closed Form)",
    NNR = "NNR (Closed Form)"
  ) %>%
  tab_header(title = "Aim 3: Closed-Form NNS and NNR") %>%
  tab_options(table.width = "100%")

gtsave(gt_tbl, output_html)
message("Saved: ", output_tsv)
message("Saved: ", output_html)
