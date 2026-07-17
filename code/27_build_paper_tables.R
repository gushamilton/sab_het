## =====================================================================
##  27_build_paper_tables.R -- Consolidated paper table pack
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table

accuracy <- as.numeric(Sys.getenv("ACCURACY", unset = "1.00"))
acc_label <- sprintf("acc%03d", as.integer(round(accuracy * 100)))
run_tag <- Sys.getenv("RUN_TAG", unset = "")
tag_suffix <- if (nzchar(run_tag)) paste0("_", run_tag) else ""

paper_tables_dir <- file.path("results", "paper", "tables")
dir.create(paper_tables_dir, showWarnings = FALSE, recursive = TRUE)

main_metrics_path <- file.path("results", "main", "data", "main_binary_metrics.tsv")
supp1_path <- file.path("results", "supp", "tables", "supp1_enrichment_nns_nnr.tsv")
supp2_prev_path <- file.path("results", "supp", "tables", "supp2_cohort_prevalence_mortality.tsv")
supp2_req_path <- file.path("results", "supp", "tables", "supp2_required_n_by_cohort.tsv")
supp3_metrics_path <- file.path("results", "supp", "tables", paste0("supp3_binary_vs_ordinal_metrics_", acc_label, tag_suffix, ".tsv"))
supp3_type1_path <- file.path("results", "supp", "tables", paste0("supp3_type1_binary_vs_ordinal_", acc_label, tag_suffix, ".tsv"))
supp3_bias_path <- file.path("results", "supp", "tables", paste0("supp3_bias_binary_vs_ordinal_", acc_label, tag_suffix, ".tsv"))
supp3_gain_path <- file.path("results", "supp", "tables", paste0("supp3_power_gain_", acc_label, tag_suffix, ".tsv"))
supp3_tradeoff_path <- file.path("results", "supp", "tables", paste0("supp3_power_vs_bias_tradeoff_", acc_label, tag_suffix, ".tsv"))

required <- c(
  main_metrics_path, supp1_path, supp2_prev_path, supp2_req_path,
  supp3_metrics_path, supp3_type1_path, supp3_bias_path, supp3_gain_path, supp3_tradeoff_path
)

missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop(
    "Cannot build paper table pack. Missing inputs:\n",
    paste0(" - ", missing, collapse = "\n"),
    call. = FALSE
  )
}

main_metrics <- readr::read_tsv(main_metrics_path, show_col_types = FALSE)
supp1 <- readr::read_tsv(supp1_path, show_col_types = FALSE)
supp2_prev <- readr::read_tsv(supp2_prev_path, show_col_types = FALSE)
supp2_req <- readr::read_tsv(supp2_req_path, show_col_types = FALSE)
supp3_metrics <- readr::read_tsv(supp3_metrics_path, show_col_types = FALSE)
supp3_type1 <- readr::read_tsv(supp3_type1_path, show_col_types = FALSE)
supp3_bias <- readr::read_tsv(supp3_bias_path, show_col_types = FALSE)
supp3_gain <- readr::read_tsv(supp3_gain_path, show_col_types = FALSE)
supp3_tradeoff <- readr::read_tsv(supp3_tradeoff_path, show_col_types = FALSE)

table1 <- subphenotypes %>%
  transmute(
    subphenotype,
    label,
    prevalence,
    baseline_mortality,
    or_arrest_shrunk
  ) %>%
  arrange(subphenotype)
readr::write_tsv(table1, file.path(paper_tables_dir, "table1_main_parameters.tsv"))

table2 <- main_metrics %>%
  filter(
    alpha == params$alpha_primary,
    sample_size %in% c(5000, 10000)
  ) %>%
  transmute(
    sample_size,
    accuracy,
    subphenotype = group,
    group_label,
    type1,
    power,
    type_s,
    type_m
  ) %>%
  arrange(sample_size, accuracy, subphenotype)
readr::write_tsv(table2, file.path(paper_tables_dir, "table2_main_binary_operating_characteristics.tsv"))

table3 <- main_metrics %>%
  filter(alpha == params$alpha_primary, sample_size == 10000) %>%
  transmute(
    accuracy,
    subphenotype = group,
    group_label,
    type1,
    power,
    type_s,
    type_m
  ) %>%
  arrange(subphenotype, accuracy)
readr::write_tsv(table3, file.path(paper_tables_dir, "table3_accuracy_impact_n10000.tsv"))

tableS1 <- supp1 %>%
  mutate(
    observed_or = exp(observed_log_or),
    true_or = exp(true_log_or)
  ) %>%
  select(
    test_type, sensitivity, specificity, target_group,
    n_required, nns, enrol_rate, true_or, observed_or, bias_log_or
  ) %>%
  arrange(target_group, desc(sensitivity), desc(specificity))
readr::write_tsv(tableS1, file.path(paper_tables_dir, "tableS1_enrichment_summary.tsv"))

tableS2 <- supp2_req %>%
  left_join(
    supp2_prev %>% select(cohort, subphenotype, prevalence, mortality),
    by = c("cohort", "subphenotype")
  ) %>%
  mutate(
    n_total_required = n_required / prevalence
  ) %>%
  select(cohort, subphenotype, prevalence, mortality, n_required, n_total_required) %>%
  arrange(cohort, subphenotype)
readr::write_tsv(tableS2, file.path(paper_tables_dir, "tableS2_cohort_required_n.tsv"))

tableS3 <- supp3_metrics %>%
  filter(sample_size %in% c(2500, 5000, 10000)) %>%
  select(sample_size, outcome_type, group, type1, power, type_s, type_m, accuracy) %>%
  arrange(sample_size, outcome_type, group)
readr::write_tsv(tableS3, file.path(paper_tables_dir, "tableS3_ordinal_binary_metrics_key_n.tsv"))

tableS4 <- supp3_type1 %>%
  left_join(
    supp3_gain %>% filter(sample_size %in% c(2500, 5000, 10000)),
    by = "sample_size"
  ) %>%
  left_join(
    supp3_bias %>%
      filter(sample_size %in% c(2500, 5000, 10000)) %>%
      select(sample_size, outcome_type, group, mean_abs_bias_log_or),
    by = c("sample_size", "outcome_type", "group")
  ) %>%
  arrange(sample_size, outcome_type, group)
readr::write_tsv(tableS4, file.path(paper_tables_dir, "tableS4_ordinal_summary_key_n.tsv"))

tableS5 <- supp3_tradeoff %>%
  filter(sample_size %in% c(2500, 5000, 10000)) %>%
  arrange(sample_size, group)
readr::write_tsv(tableS5, file.path(paper_tables_dir, "tableS5_power_bias_tradeoff_key_n.tsv"))

table_index <- tibble::tribble(
  ~table_id, ~path, ~purpose,
  "Table 1", "results/paper/tables/table1_main_parameters.tsv", "Primary simulation inputs and shrunk effect assumptions",
  "Table 2", "results/paper/tables/table2_main_binary_operating_characteristics.tsv", "Main operating characteristics at N=5,000 and N=10,000",
  "Table 3", "results/paper/tables/table3_accuracy_impact_n10000.tsv", "Accuracy sensitivity at fixed N=10,000",
  "Table S1", "results/paper/tables/tableS1_enrichment_summary.tsv", "Enrichment feasibility, NNR/NNS and dilution",
  "Table S2", "results/paper/tables/tableS2_cohort_required_n.tsv", "Cohort prevalence/mortality and required subgroup/total N",
  "Table S3", "results/paper/tables/tableS3_ordinal_binary_metrics_key_n.tsv", "Binary vs ordinal metrics at key sample sizes",
  "Table S4", "results/paper/tables/tableS4_ordinal_summary_key_n.tsv", "Combined ordinal summary: Type I, power gain, and bias",
  "Table S5", "results/paper/tables/tableS5_power_bias_tradeoff_key_n.tsv", "Power-gain vs bias-change trade-off"
)
readr::write_tsv(table_index, file.path(paper_tables_dir, "table_index.tsv"))
