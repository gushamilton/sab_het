## =====================================================================
##  07_build_gt_tables.R -- Render key TSV outputs as gt HTML tables
## =====================================================================

library(tidyverse)

if (!requireNamespace("gt", quietly = TRUE)) {
  message("Package 'gt' is not installed; skipping HTML table rendering.")
  quit(save = "no", status = 0)
}

library(gt)

dir.create("results/main/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/supp/tables", showWarnings = FALSE, recursive = TRUE)

render_gt <- function(input_path, output_path, title, subtitle = NULL, digits = 3) {
  if (!file.exists(input_path)) {
    message(sprintf("Missing input table: %s", input_path))
    return(invisible(NULL))
  }

  tab <- readr::read_tsv(input_path, show_col_types = FALSE) %>%
    gt() %>%
    tab_header(title = title, subtitle = subtitle) %>%
    fmt_number(columns = where(is.numeric), decimals = digits)

  gtsave(tab, output_path)
  message(sprintf("Wrote %s", output_path))
}

render_gt(
  input_path = "results/main/tables/main_binary_summary_N5000_N10000.tsv",
  output_path = "results/main/tables/tableR1_main_binary_summary.html",
  title = md("**R1 Main Analysis: Binary Outcome Summary**"),
  subtitle = "Type I, power, Type S and Type M by sample size and classification accuracy"
)

render_gt(
  input_path = "results/supp/tables/supp1_enrichment_nns_nnr.tsv",
  output_path = "results/supp/tables/tableS1_enrichment_summary.html",
  title = md("**Supplement S1: Enrichment Trial Feasibility**"),
  subtitle = "Expected dilution, N required, and N needed to screen by test performance"
)

render_gt(
  input_path = "results/supp/tables/supp2_cohort_prevalence_mortality.tsv",
  output_path = "results/supp/tables/tableS2_cohort_inputs.html",
  title = md("**Supplement S2: Cohort Inputs**"),
  subtitle = "Subphenotype prevalence and mortality assumptions by cohort"
)

render_gt(
  input_path = "results/supp/tables/supp2_power_typeS_fixedN.tsv",
  output_path = "results/supp/tables/tableS2_power_typeS_fixedN.html",
  title = md("**Supplement S2: Fixed-N Cohort Results**"),
  subtitle = "Power and directional error under cohort-specific prevalence/event rates"
)

render_gt(
  input_path = "results/supp/tables/supp3_power_binary_vs_ordinal_N5000.tsv",
  output_path = "results/supp/tables/tableS3_binary_vs_ordinal_N5000.html",
  title = md("**Supplement S3: Binary vs Ordinal at N=5000**"),
  subtitle = "Power, Type S and Type M comparison"
)

render_gt(
  input_path = "results/supp/tables/supp3_ordinal_shift_delta.tsv",
  output_path = "results/supp/tables/tableS3_ordinal_shift_delta.html",
  title = md("**Supplement S3: Ordinal Distribution Shift**"),
  subtitle = "Treatment minus control change in ordinal category probabilities"
)

render_gt(
  input_path = "results/supp/tables/supp3_type1_binary_vs_ordinal.tsv",
  output_path = "results/supp/tables/tableS3_type1_binary_vs_ordinal.html",
  title = md("**Supplement S3: Type I Error (Binary vs Ordinal)**"),
  subtitle = "Null subgroup calibration across sample size"
)

render_gt(
  input_path = "results/supp/tables/supp3_bias_binary_vs_ordinal.tsv",
  output_path = "results/supp/tables/tableS3_bias_binary_vs_ordinal.html",
  title = md("**Supplement S3: Bias Summary (Binary vs Ordinal)**"),
  subtitle = "Mean signed bias, absolute bias, and RMSE of subgroup log-OR estimates"
)

render_gt(
  input_path = "results/supp/tables/supp3_power_vs_bias_tradeoff.tsv",
  output_path = "results/supp/tables/tableS3_power_vs_bias_tradeoff.html",
  title = md("**Supplement S3: Power Gain vs Bias Change**"),
  subtitle = "Checks whether ordinal power gains are accompanied by bias shifts"
)
