# Create Aim 1 power table (GT)

pacman::p_load(tidyverse, scales, gt)

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_file <- file.path(results_root, "tables", "aim1_power_summary.tsv")
output_tsv <- file.path(results_root, "tables", "aim1_power_table.tsv")
output_html <- file.path(results_root, "tables", "aim1_power_table.html")

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

aim1 <- readr::read_tsv(input_file, show_col_types = FALSE)

target_n <- c(100, 500, 1000, 2000, 3000, 5000, 7500, 10000, 15000, 20000)

wide <- aim1 %>%
  filter(scenario_name == "ARREST_raw", n_total %in% target_n, group %in% c("A", "B", "C", "D", "E")) %>%
  mutate(power = scales::percent(power, accuracy = 1)) %>%
  select(scenario_name, n_total, group, power) %>%
  pivot_wider(names_from = group, values_from = power) %>%
  arrange(n_total)

readr::write_tsv(wide, output_tsv)

gt_tbl <- wide %>%
  select(-scenario_name) %>%
  gt(rowname_col = "n_total") %>%
  cols_label(
    n_total = "Total n",
    A = "A", B = "B", C = "C", D = "D", E = "E"
  ) %>%
  tab_header(
    title = "Statistical power by sample size",
    subtitle = "Power to detect subphenotype-specific effects (Bonferroni alpha = 0.01)"
  ) %>%
  cols_align(align = "center", columns = c(A, B, C, D, E)) %>%
  cols_align(align = "right", columns = n_total)

gtsave(gt_tbl, output_html)
message("Saved: ", output_tsv)
message("Saved: ", output_html)
