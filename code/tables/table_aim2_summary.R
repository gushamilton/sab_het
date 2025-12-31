# Create Aim 2 summary table (GT with spanners)

pacman::p_load(tidyverse, gt)

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
input_file <- file.path(results_root, "tables", "aim2_accuracy_summary.tsv")
output_tsv <- file.path(results_root, "tables", "aim2_accuracy_table.tsv")
output_html <- file.path(results_root, "tables", "aim2_accuracy_table.html")

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

aim2 <- readr::read_tsv(input_file, show_col_types = FALSE)

target_acc <- c(0.70, 0.80, 0.90, 0.95, 0.99, 1.00)

formatted <- aim2 %>%
  mutate(accuracy = round(accuracy, 2)) %>%
  filter(group %in% c("B", "C", "D", "E"), accuracy %in% target_acc) %>%
  mutate(
    power = sprintf("%.1f", 100 * power),
    bias = sprintf("%.3f", bias),
    mse = sprintf("%.3f", mse),
    wrong_dir = sprintf("%.1f", wrong_dir)
  ) %>%
  select(accuracy, group, power, bias, mse, wrong_dir) %>%
  pivot_longer(cols = c(power, bias, mse, wrong_dir), names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = c(group, metric), values_from = value) %>%
  arrange(accuracy)

formatted <- formatted %>%
  rename_with(function(x) {
    x <- gsub("power", "Power", x, fixed = TRUE)
    x <- gsub("bias", "Bias", x, fixed = TRUE)
    x <- gsub("mse", "MSE", x, fixed = TRUE)
    x <- gsub("wrong_dir", "Wrong_dir", x, fixed = TRUE)
    x
  })

readr::write_tsv(formatted, output_tsv)

gt_tbl <- formatted %>%
  gt(rowname_col = "accuracy") %>%
  cols_label(
    accuracy = "Accuracy",
    B_Power = "Power",
    B_Bias = "Bias",
    B_MSE = "MSE",
    B_Wrong_dir = "Wrong dir",
    C_Power = "Power",
    C_Bias = "Bias",
    C_MSE = "MSE",
    C_Wrong_dir = "Wrong dir",
    D_Power = "Power",
    D_Bias = "Bias",
    D_MSE = "MSE",
    D_Wrong_dir = "Wrong dir",
    E_Power = "Power",
    E_Bias = "Bias",
    E_MSE = "MSE",
    E_Wrong_dir = "Wrong dir"
  ) %>%
  tab_spanner(label = "B", columns = starts_with("B_")) %>%
  tab_spanner(label = "C", columns = starts_with("C_")) %>%
  tab_spanner(label = "D", columns = starts_with("D_")) %>%
  tab_spanner(label = "E", columns = starts_with("E_")) %>%
  tab_header(
    title = "Aim 2: Accuracy, power, bias, and error rates by subphenotype"
  ) %>%
  cols_align(align = "center", columns = -accuracy)

gtsave(gt_tbl, output_html)
message("Saved: ", output_tsv)
message("Saved: ", output_html)
