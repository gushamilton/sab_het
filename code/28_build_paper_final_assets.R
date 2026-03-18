## =====================================================================
##  28_build_paper_final_assets.R -- Manuscript-ordered assets
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")
params <- get_parameters()
subphenotypes <- params$subphenotype_table

acc <- as.numeric(Sys.getenv("PAPER_ACCURACY", unset = "1.00"))
acc_label <- sprintf("acc%03d", as.integer(round(acc * 100)))
run_tag <- Sys.getenv("RUN_TAG", unset = "")
tag_suffix <- if (nzchar(run_tag)) paste0("_", run_tag) else ""
nonpo_tag <- Sys.getenv("NONPO_RUN_TAG", unset = "paper_final")

out_root <- file.path("results", "paper", "final")
fig_out <- file.path(out_root, "figures")
tab_out <- file.path(out_root, "tables")
dir.create(fig_out, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_out, showWarnings = FALSE, recursive = TRUE)

safe_copy <- function(src, dst) {
  if (!file.exists(src)) stop("Missing source file: ", src, call. = FALSE)
  ok <- file.copy(src, dst, overwrite = TRUE)
  if (!ok) stop("Failed to copy: ", src, " -> ", dst, call. = FALSE)
}

# ---------------------------------------------------------------------
# Main figure order (Figure 2/3 swapped intentionally to match paper)
# ---------------------------------------------------------------------
figure_map <- tibble::tribble(
  ~figure_id, ~source, ~target, ~purpose,
  "Figure 1", "results/core_figures/plots/figure1_accuracy_dilution_patchwork.pdf", "Figure1_misclassification_accuracy.pdf", "Impact of classification accuracy on power and bias",
  "Figure 2", "results/core_figures/plots/figure3.pdf", "Figure2_cohort_conditionality.pdf", "Cohort prevalence/mortality conditionality and required N",
  "Figure 3", "results/core_figures/plots/figure2_enrichment_master.pdf", "Figure3_enrichment_feasibility.pdf", "Enrichment feasibility (NNS/NNR/bias)",
  "Figure 4", "paper/ordinal_outcome_two_panel.pdf", "Figure4_ordinal_category_shift.pdf", "Ordinal outcome category shifts under PO and nonPO",
  "Figure 5", "results/supp/ordinal_nonPO_comparison/plots/ordinal_nonPO_patchwork_paper_final_independent_20260313.pdf", "Figure5_ordinal_binary_summary.pdf", "Ordinal versus binary outcome performance (PO vs nonPO patchwork)"
)

purrr::pwalk(figure_map, function(figure_id, source, target, purpose) {
  safe_copy(source, file.path(fig_out, target))
})
readr::write_tsv(figure_map, file.path(fig_out, "figure_order.tsv"))

# ---------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------
main_metrics <- readr::read_tsv("results/main/data/main_binary_metrics.tsv", show_col_types = FALSE)
supp1 <- readr::read_tsv("results/supp/tables/supp1_enrichment_nns_nnr.tsv", show_col_types = FALSE)
supp2_prev <- readr::read_tsv("results/supp/tables/supp2_cohort_prevalence_mortality.tsv", show_col_types = FALSE)
supp2_req <- readr::read_tsv("results/supp/tables/supp2_required_n_by_cohort.tsv", show_col_types = FALSE)

supp3_metrics <- readr::read_tsv(
  file.path("results", "supp", "tables", paste0("supp3_binary_vs_ordinal_metrics_", acc_label, tag_suffix, ".tsv")),
  show_col_types = FALSE
)
supp3_type1 <- readr::read_tsv(
  file.path("results", "supp", "tables", paste0("supp3_type1_binary_vs_ordinal_", acc_label, tag_suffix, ".tsv")),
  show_col_types = FALSE
)
supp3_bias <- readr::read_tsv(
  file.path("results", "supp", "tables", paste0("supp3_bias_binary_vs_ordinal_", acc_label, tag_suffix, ".tsv")),
  show_col_types = FALSE
)
supp3_tradeoff <- readr::read_tsv(
  file.path("results", "supp", "tables", paste0("supp3_power_vs_bias_tradeoff_", acc_label, tag_suffix, ".tsv")),
  show_col_types = FALSE
)

nonpo_metrics <- readr::read_tsv(
  file.path("results", "supp", "ordinal_nonPO_comparison", "tables", paste0("ordinal_nonPO_comparison_metrics_", nonpo_tag, ".tsv")),
  show_col_types = FALSE
)
nonpo_bias <- readr::read_tsv(
  file.path("results", "supp", "ordinal_nonPO_comparison", "tables", paste0("ordinal_nonPO_comparison_bias_", nonpo_tag, ".tsv")),
  show_col_types = FALSE
)
nonpo_cal_path <- file.path(
  "results", "supp", "ordinal_nonPO_comparison", "tables",
  paste0("ordinal_nonPO_comparison_calibration_", nonpo_tag, ".tsv")
)
if (file.exists(nonpo_cal_path)) {
  nonpo_cal <- readr::read_tsv(nonpo_cal_path, show_col_types = FALSE)
} else {
  warning("Calibration table not found for run tag ", nonpo_tag, "; writing empty Table S7.")
  nonpo_cal <- tibble::tibble(
    dgm = character(),
    group = character(),
    model_type = character(),
    log_or_true = double(),
    log_or_hat = double(),
    raw_bias = double(),
    prop_bias = double()
  )
}

available_nonpo_n <- sort(unique(nonpo_metrics$sample_size))
if (3000 %in% available_nonpo_n) {
  nonpo_key_n <- 3000
} else if (2500 %in% available_nonpo_n) {
  nonpo_key_n <- 2500
} else if (5000 %in% available_nonpo_n) {
  nonpo_key_n <- 5000
} else {
  nonpo_key_n <- min(available_nonpo_n, na.rm = TRUE)
}

# ---------------------------------------------------------------------
# Main tables
# ---------------------------------------------------------------------
table1 <- subphenotypes %>%
  transmute(
    subphenotype,
    label,
    prevalence,
    baseline_mortality,
    or_arrest_raw,
    or_arrest_shrunk
  ) %>%
  arrange(subphenotype)
readr::write_tsv(table1, file.path(tab_out, "Table1_main_parameters.tsv"))

# Main Table 2: ordinal-focused PO vs nonPO comparison (power + bias)
table2 <- nonpo_metrics %>%
  filter(group != "A", sample_size == nonpo_key_n, accuracy %in% c(0.7, 1.0)) %>%
  left_join(
    nonpo_bias %>%
      select(sample_size, accuracy, dgm, group, model_type, mean_prop_bias, rmse_log),
    by = c("sample_size", "accuracy", "dgm", "group", "model_type")
  ) %>%
  arrange(dgm, accuracy, group, model_type) %>%
  select(
    sample_size, accuracy, dgm, group, model_type,
    power, type_s, type_m, mean_prop_bias, rmse_log
  )
readr::write_tsv(table2, file.path(tab_out, "Table2_ordinal_PO_vs_nonPO_power_bias.tsv"))

# ---------------------------------------------------------------------
# Supplementary tables (ordered to follow manuscript flow)
# ---------------------------------------------------------------------
tableS1 <- main_metrics %>%
  filter(alpha == params$alpha_primary, sample_size %in% c(5000, 10000)) %>%
  mutate(group_label = subphenotypes$label[match(group, subphenotypes$subphenotype)]) %>%
  select(sample_size, accuracy, group, group_label, type1, power, type_s, type_m) %>%
  arrange(sample_size, accuracy, group)
readr::write_tsv(tableS1, file.path(tab_out, "TableS1_main_binary_operating_characteristics.tsv"))

tableS2 <- supp2_req %>%
  left_join(
    supp2_prev %>% select(cohort, subphenotype, prevalence, mortality),
    by = c("cohort", "subphenotype")
  ) %>%
  mutate(n_total_required = n_required / prevalence) %>%
  select(cohort, subphenotype, prevalence, mortality, n_required, n_total_required) %>%
  arrange(cohort, subphenotype)
readr::write_tsv(tableS2, file.path(tab_out, "TableS2_cohort_required_n.tsv"))

tableS3 <- supp1 %>%
  mutate(
    observed_or = exp(observed_log_or),
    true_or = exp(true_log_or)
  ) %>%
  select(
    test_type, sensitivity, specificity, target_group,
    n_required, nns, enrol_rate, true_or, observed_or, bias_log_or
  ) %>%
  arrange(target_group, desc(sensitivity), desc(specificity))
readr::write_tsv(tableS3, file.path(tab_out, "TableS3_enrichment_summary.tsv"))

tableS4 <- supp3_metrics %>%
  filter(sample_size %in% c(2500, 5000, 10000)) %>%
  left_join(
    supp3_bias %>%
      select(sample_size, outcome_type, group, mean_abs_bias_log_or, mean_bias_log_or),
    by = c("sample_size", "outcome_type", "group")
  ) %>%
  arrange(sample_size, outcome_type, group) %>%
  select(sample_size, outcome_type, group, type1, power, type_s, type_m, mean_abs_bias_log_or, mean_bias_log_or)
readr::write_tsv(tableS4, file.path(tab_out, "TableS4_ordinal_binary_metrics_bias_key_n.tsv"))

tableS5 <- supp3_type1 %>%
  arrange(sample_size, outcome_type)
readr::write_tsv(tableS5, file.path(tab_out, "TableS5_ordinal_type1.tsv"))

tableS6 <- supp3_tradeoff %>%
  filter(sample_size %in% c(2500, 5000, 10000)) %>%
  arrange(sample_size, group)
readr::write_tsv(tableS6, file.path(tab_out, "TableS6_ordinal_power_bias_tradeoff_key_n.tsv"))

tableS7 <- nonpo_cal %>%
  arrange(dgm, group, model_type)
readr::write_tsv(tableS7, file.path(tab_out, "TableS7_ordinal_PO_nonPO_calibration.tsv"))

table_map <- tibble::tribble(
  ~table_id, ~path, ~purpose,
  "Table 1", "results/paper/final/tables/Table1_main_parameters.tsv", "Main simulation assumptions",
  "Table 2", "results/paper/final/tables/Table2_ordinal_PO_vs_nonPO_power_bias.tsv", "Main ordinal result: PO vs nonPO, power and bias",
  "Table S1", "results/paper/final/tables/TableS1_main_binary_operating_characteristics.tsv", "Detailed main binary operating characteristics",
  "Table S2", "results/paper/final/tables/TableS2_cohort_required_n.tsv", "Cohort-conditional required N",
  "Table S3", "results/paper/final/tables/TableS3_enrichment_summary.tsv", "Enrichment feasibility detail",
  "Table S4", "results/paper/final/tables/TableS4_ordinal_binary_metrics_bias_key_n.tsv", "Ordinal vs binary key-N metrics + bias",
  "Table S5", "results/paper/final/tables/TableS5_ordinal_type1.tsv", "Ordinal/binary Type I calibration",
  "Table S6", "results/paper/final/tables/TableS6_ordinal_power_bias_tradeoff_key_n.tsv", "Ordinal power-bias tradeoff",
  "Table S7", "results/paper/final/tables/TableS7_ordinal_PO_nonPO_calibration.tsv", "PO vs nonPO model calibration"
)
readr::write_tsv(table_map, file.path(tab_out, "table_order.tsv"))

# Keep provenance pointers to canonical raw sources for secondary analyses.
raw_provenance <- tibble::tribble(
  ~dataset, ~path,
  "Main binary raw", "results/main/data/main_binary_raw.tsv",
  "Main binary metrics", "results/main/data/main_binary_metrics.tsv",
  "Supplement 3 raw", file.path("results", "supp", "data", paste0("supp3_binary_vs_ordinal_raw_", acc_label, tag_suffix, ".tsv.gz")),
  "Supplement 3 metrics", file.path("results", "supp", "tables", paste0("supp3_binary_vs_ordinal_metrics_", acc_label, tag_suffix, ".tsv")),
  "PO/nonPO raw", file.path("results", "supp", "ordinal_nonPO_comparison", "data", paste0("ordinal_nonPO_comparison_raw_", nonpo_tag, ".tsv.gz"))
)
readr::write_tsv(raw_provenance, file.path(tab_out, "raw_data_provenance.tsv"))
