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
nonpo_tag <- Sys.getenv("NONPO_RUN_TAG", unset = "ordinal6_figure5_20260325")

out_root <- file.path("results", "paper", "final")
fig_out <- file.path(out_root, "figures")
tab_out <- file.path(out_root, "tables")
dir.create(fig_out, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_out, showWarnings = FALSE, recursive = TRUE)

invisible(file.remove(
  file.path(tab_out, "Table2_ordinal_PO_vs_nonPO_power_bias.tsv"),
  file.path(fig_out, "Figure4_ordinal_category_shift_6point_subgroupB.pdf")
))

safe_copy <- function(src, dst) {
  if (!file.exists(src)) stop("Missing source file: ", src, call. = FALSE)
  ok <- file.copy(src, dst, overwrite = TRUE)
  if (!ok) stop("Failed to copy: ", src, " -> ", dst, call. = FALSE)
}

fmt_int <- function(x) {
  ifelse(is.na(x), "", scales::comma(as.integer(round(x))))
}

fmt_prob3 <- function(x) {
  ifelse(is.na(x), "", sprintf("%.3f", x))
}

fmt_pct1 <- function(x) {
  ifelse(is.na(x), "", sprintf("%.1f", 100 * x))
}

fmt_num2 <- function(x) {
  ifelse(is.na(x), "", sprintf("%.2f", x))
}

fmt_num3 <- function(x) {
  ifelse(is.na(x), "", sprintf("%.3f", x))
}

fmt_num2_or_inf <- function(x) {
  ifelse(is.infinite(x), "Inf", ifelse(is.na(x), "", sprintf("%.2f", x)))
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
  "Figure 5", file.path("results", "supp", "ordinal_nonPO_comparison", "plots", paste0("ordinal_nonPO_patchwork_", nonpo_tag, ".pdf")), "Figure5_ordinal_binary_summary.pdf", "Ordinal versus binary outcome performance (PO vs nonPO patchwork)",
  "Figure S1", "results/paper/final/figures/FigureS1_ordinal_power_5pt_vs_6pt.pdf", "FigureS1_ordinal_power_5pt_vs_6pt.pdf", "Sensitivity of ordinal power to 5-point versus 6-point outcome scales"
)

purrr::pwalk(figure_map %>% filter(!str_detect(source, "^results/paper/final/figures/")), function(figure_id, source, target, purpose) {
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
    `Subphenotype` = subphenotype,
    `Clinical label` = label,
    `Frequency in ARREST placebo arm (%)` = fmt_pct1(prevalence),
    `Baseline 84-day mortality across both arms (%)` = fmt_pct1(baseline_mortality),
    `Original OR for 84-day mortality` = fmt_num2(or_arrest_raw),
    `Conservative OR for 84-day mortality` = fmt_num2(or_arrest_shrunk)
  ) %>%
  arrange(`Subphenotype`)
readr::write_tsv(table1, file.path(tab_out, "Table1_main_parameters.tsv"))

table2_main <- main_metrics %>%
  filter(sample_size %in% c(500, 1000, 2000, 3000, 5000, 10000, 20000)) %>%
  filter(alpha == params$alpha_primary, accuracy == 1.0, group != "A") %>%
  select(sample_size, group, power) %>%
  mutate(power = fmt_prob3(power)) %>%
  pivot_wider(names_from = group, values_from = power) %>%
  mutate(
    A = "-",
    `Total trial size` = fmt_int(sample_size)
  ) %>%
  select(`Total trial size`, A, B, C, D, E) %>%
  arrange(as.integer(gsub(",", "", `Total trial size`)))
readr::write_tsv(table2_main, file.path(tab_out, "Table2_main_power_grid.tsv"))

# ---------------------------------------------------------------------
# Supplementary tables
# ---------------------------------------------------------------------
tableS1 <- main_metrics %>%
  filter(alpha == params$alpha_primary, sample_size %in% c(5000, 10000)) %>%
  mutate(group_label = subphenotypes$label[match(group, subphenotypes$subphenotype)]) %>%
  transmute(
    `Total trial size` = fmt_int(sample_size),
    `Classification accuracy` = if_else(accuracy == 1.0, "100%", paste0(round(100 * accuracy), "%")),
    `Subphenotype` = group,
    `Clinical label` = group_label,
    `Type I error` = fmt_prob3(type1),
    `Power` = fmt_prob3(power),
    `Type S error` = fmt_prob3(type_s),
    `Type M error` = fmt_num2(type_m)
  ) %>%
  arrange(as.integer(gsub(",", "", `Total trial size`)), `Classification accuracy`, `Subphenotype`)
readr::write_tsv(tableS1, file.path(tab_out, "TableS1_main_binary_operating_characteristics.tsv"))

tableS2 <- supp2_req %>%
  left_join(
    supp2_prev %>% select(cohort, subphenotype, prevalence, mortality),
    by = c("cohort", "subphenotype")
  ) %>%
  mutate(n_total_required = n_required / prevalence) %>%
  transmute(
    `Cohort` = cohort,
    `Subphenotype` = subphenotype,
    `Prevalence (%)` = fmt_pct1(prevalence),
    `Baseline mortality (%)` = fmt_pct1(mortality),
    `Required N within subgroup` = fmt_num2_or_inf(n_required),
    `Required total trial N` = fmt_num2_or_inf(n_total_required)
  ) %>%
  arrange(`Cohort`, `Subphenotype`)
readr::write_tsv(tableS2, file.path(tab_out, "TableS2_cohort_required_n.tsv"))

tableS3 <- supp1 %>%
  mutate(
    observed_or = exp(observed_log_or),
    true_or = exp(true_log_or),
    sensitivity = fmt_pct1(sensitivity),
    specificity = fmt_pct1(specificity),
    n_required = fmt_int(n_required),
    nns = fmt_int(nns),
    enrol_rate = fmt_prob3(enrol_rate),
    true_or = fmt_num2(true_or),
    observed_or = fmt_num2(observed_or),
    bias_log_or = fmt_num3(bias_log_or)
  ) %>%
  transmute(
    `Test-performance scenario` = test_type,
    `Sensitivity (%)` = sensitivity,
    `Specificity (%)` = specificity,
    `Target subphenotype` = target_group,
    `Required N randomised` = n_required,
    `Number needed to screen` = nns,
    `Enrolment rate` = enrol_rate,
    `True OR` = true_or,
    `Observed OR after enrichment` = observed_or,
    `Bias in log-OR` = bias_log_or
  ) %>%
  arrange(`Target subphenotype`, desc(`Sensitivity (%)`), desc(`Specificity (%)`))
readr::write_tsv(tableS3, file.path(tab_out, "TableS3_enrichment_summary.tsv"))

tableS4 <- supp3_metrics %>%
  filter(sample_size %in% c(2500, 5000, 10000)) %>%
  left_join(
    supp3_bias %>%
      select(sample_size, outcome_type, group, mean_abs_bias_log_or, mean_bias_log_or),
    by = c("sample_size", "outcome_type", "group")
  ) %>%
  mutate(
    outcome_type = recode(outcome_type, Binary = "Binary (death)", Ordinal = "Ordinal")
  ) %>%
  arrange(sample_size, outcome_type, group) %>%
  transmute(
    `Total trial size` = fmt_int(sample_size),
    `Outcome type` = outcome_type,
    `Subphenotype` = group,
    `Type I error` = fmt_prob3(type1),
    `Power` = fmt_prob3(power),
    `Type S error` = fmt_prob3(type_s),
    `Type M error` = fmt_num2(type_m),
    `Mean absolute bias in log-OR` = fmt_num3(mean_abs_bias_log_or),
    `Mean signed bias in log-OR` = fmt_num3(mean_bias_log_or)
  )
readr::write_tsv(tableS4, file.path(tab_out, "TableS4_ordinal_binary_metrics_bias_key_n.tsv"))

tableS5 <- supp3_type1 %>%
  mutate(
    outcome_type = recode(outcome_type, Binary = "Binary (death)", Ordinal = "Ordinal")
  ) %>%
  arrange(sample_size, outcome_type) %>%
  transmute(
    `Total trial size` = fmt_int(sample_size),
    `Outcome type` = outcome_type,
    `Type I error` = fmt_prob3(type1),
    `Classification accuracy` = if_else(accuracy == 1.0, "100%", paste0(round(100 * accuracy), "%"))
  )
readr::write_tsv(tableS5, file.path(tab_out, "TableS5_ordinal_type1.tsv"))

tableS6 <- supp3_tradeoff %>%
  filter(sample_size %in% c(2500, 5000, 10000)) %>%
  arrange(sample_size, group) %>%
  transmute(
    `Total trial size` = fmt_int(sample_size),
    `Subphenotype` = group,
    `Binary power` = fmt_prob3(Binary),
    `Ordinal power` = fmt_prob3(Ordinal),
    `Absolute power gain` = fmt_prob3(power_gain_abs),
    `Relative power gain` = fmt_num2(power_gain_rel),
    `Binary mean absolute bias in log-OR` = fmt_num3(mean_abs_bias_log_or_Binary),
    `Ordinal mean absolute bias in log-OR` = fmt_num3(mean_abs_bias_log_or_Ordinal),
    `Binary mean signed bias in log-OR` = fmt_num3(mean_bias_log_or_Binary),
    `Ordinal mean signed bias in log-OR` = fmt_num3(mean_bias_log_or_Ordinal),
    `Change in absolute bias` = fmt_num3(delta_abs_bias),
    `Change in signed bias` = fmt_num3(delta_signed_bias)
  )
readr::write_tsv(tableS6, file.path(tab_out, "TableS6_ordinal_power_bias_tradeoff_key_n.tsv"))

tableS7 <- nonpo_cal %>%
  mutate(
    dgm = recode(dgm, PO = "Proportional-odds", nonPO = "Death-only non-proportional"),
    prop_bias = if_else(abs(log_or_true) < 1e-12, NA_real_, prop_bias)
  ) %>%
  arrange(dgm, group, model_type) %>%
  transmute(
    `Data-generating mechanism` = dgm,
    `Subphenotype` = group,
    `Model` = model_type,
    `True log-OR` = fmt_num3(log_or_true),
    `Estimated log-OR` = fmt_num3(log_or_hat),
    `Raw bias` = fmt_num3(raw_bias),
    `Proportional bias` = ifelse(is.na(prop_bias), "", fmt_num3(prop_bias))
  )
readr::write_tsv(tableS7, file.path(tab_out, "TableS7_ordinal_PO_nonPO_calibration.tsv"))

tableS8 <- nonpo_metrics %>%
  filter(group != "A", sample_size == nonpo_key_n, accuracy %in% c(0.7, 1.0)) %>%
  left_join(
    nonpo_bias %>%
      select(sample_size, accuracy, dgm, group, model_type, mean_prop_bias, rmse_log),
    by = c("sample_size", "accuracy", "dgm", "group", "model_type")
  ) %>%
  mutate(
    accuracy = if_else(accuracy == 1.0, "100%", "70%"),
    dgm = recode(dgm, PO = "Proportional-odds", nonPO = "Death-only non-proportional")
  ) %>%
  arrange(dgm, accuracy, group, model_type) %>%
  transmute(
    `Total trial size` = fmt_int(sample_size),
    `Classification accuracy` = accuracy,
    `Data-generating mechanism` = dgm,
    `Subphenotype` = group,
    `Model` = model_type,
    `Power` = fmt_prob3(power),
    `Type S error` = fmt_prob3(type_s),
    `Type M error` = fmt_num2(type_m),
    `Proportional bias` = fmt_num3(mean_prop_bias),
    `RMSE on log-OR scale` = fmt_num3(rmse_log)
  )
readr::write_tsv(tableS8, file.path(tab_out, "TableS8_ordinal_PO_vs_nonPO_power_bias.tsv"))

table_map <- tibble::tribble(
  ~table_id, ~path, ~purpose,
  "Table 1", "results/paper/final/tables/Table1_main_parameters.tsv", "Main simulation assumptions",
  "Table 2", "results/paper/final/tables/Table2_main_power_grid.tsv", "Main binary power grid used in manuscript text",
  "Table S1", "results/paper/final/tables/TableS1_main_binary_operating_characteristics.tsv", "Detailed main binary operating characteristics",
  "Table S2", "results/paper/final/tables/TableS2_cohort_required_n.tsv", "Cohort-conditional required N",
  "Table S3", "results/paper/final/tables/TableS3_enrichment_summary.tsv", "Enrichment feasibility detail",
  "Table S4", "results/paper/final/tables/TableS4_ordinal_binary_metrics_bias_key_n.tsv", "Ordinal vs binary key-N metrics + bias",
  "Table S5", "results/paper/final/tables/TableS5_ordinal_type1.tsv", "Ordinal/binary Type I calibration",
  "Table S6", "results/paper/final/tables/TableS6_ordinal_power_bias_tradeoff_key_n.tsv", "Ordinal power-bias tradeoff",
  "Table S7", "results/paper/final/tables/TableS7_ordinal_PO_nonPO_calibration.tsv", "PO vs nonPO model calibration",
  "Table S8", "results/paper/final/tables/TableS8_ordinal_PO_vs_nonPO_power_bias.tsv", "Ordinal PO vs nonPO power and bias comparison"
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
