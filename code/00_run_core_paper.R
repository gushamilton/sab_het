## =====================================================================
##  00_run_core_paper.R -- End-to-end core paper pipeline
## =====================================================================

library(tidyverse)

run_with_env <- function(script_path, env = list()) {
  old <- Sys.getenv(names(env), unset = NA_character_)
  on.exit({
    for (nm in names(env)) {
      prev <- old[[nm]]
      if (is.na(prev)) {
        Sys.unsetenv(nm)
      } else {
        Sys.setenv(structure(prev, names = nm))
      }
    }
  }, add = TRUE)

  if (length(env) > 0) {
    do.call(Sys.setenv, env)
  }

  message(sprintf("[core] Running %s", script_path))
  source(script_path, local = new.env(parent = globalenv()))
}

check_outputs <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "Missing expected outputs:\n",
      paste0(" - ", missing, collapse = "\n"),
      call. = FALSE
    )
  }
}

quick_mode <- tolower(Sys.getenv("CORE_QUICK", unset = "0")) %in% c("1", "true", "yes")
run_acc <- as.numeric(Sys.getenv("CORE_ACCURACY", unset = "1.00"))

dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("results", "paper"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("results", "paper", "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("results", "paper", "figures"), showWarnings = FALSE, recursive = TRUE)

if (quick_mode) {
  message("[core] Quick mode enabled.")
}

main_env <- list(
  N_REPS_MAIN = if (quick_mode) "50" else Sys.getenv("N_REPS_MAIN", unset = "2000"),
  SAMPLE_SIZES = if (quick_mode) "500,1000" else Sys.getenv("SAMPLE_SIZES", unset = "")
)

cohort_env <- list(
  N_REPS_COHORT = if (quick_mode) "50" else Sys.getenv("N_REPS_COHORT", unset = "1000"),
  COHORT_FIXED_N = if (quick_mode) "1000" else Sys.getenv("COHORT_FIXED_N", unset = "5000")
)

fig1_env <- list(
  FIG1_N_REPS = if (quick_mode) "50" else Sys.getenv("FIG1_N_REPS", unset = "1000"),
  FIG1_N = if (quick_mode) "1000" else Sys.getenv("FIG1_N", unset = "20000")
)

ordinal_env <- list(
  ACCURACY = sprintf("%.2f", run_acc),
  N_REPS_ORDINAL = if (quick_mode) "50" else Sys.getenv("N_REPS_ORDINAL", unset = "1000"),
  ORDINAL_SAMPLE_SIZES = if (quick_mode) "500,1000,2000" else Sys.getenv("ORDINAL_SAMPLE_SIZES", unset = ""),
  RUN_ORDER = Sys.getenv("RUN_ORDER", unset = "binary,ordinal"),
  RUN_TAG = Sys.getenv("RUN_TAG", unset = "")
)

run_with_env("code/03_run_main_binary.R", env = main_env)
run_with_env("code/04_run_enrichment.R")
run_with_env("code/05_run_cohort_variation.R", env = cohort_env)
run_with_env("code/11_run_figure1_accuracy_scan.R", env = fig1_env)
run_with_env("code/12_plot_figure2_master.R")
run_with_env("code/13_plot_figure3_cohort_summary.R")
run_with_env("code/06_run_ordinal_comparison.R", env = ordinal_env)
run_with_env("code/26_build_core_figure4_ordinal.R", env = list(ACCURACY = sprintf("%.2f", run_acc)))
run_with_env("code/27_build_paper_tables.R", env = list(ACCURACY = sprintf("%.2f", run_acc)))

acc_label <- sprintf("acc%03d", as.integer(round(run_acc * 100)))
tag_suffix <- if (nzchar(Sys.getenv("RUN_TAG", unset = ""))) paste0("_", Sys.getenv("RUN_TAG")) else ""

expected_outputs <- c(
  "results/core_figures/plots/figure1_accuracy_dilution_patchwork.pdf",
  "results/core_figures/plots/figure2_enrichment_master.pdf",
  "results/core_figures/plots/figure3.pdf",
  "results/core_figures/plots/figure4.pdf",
  "results/main/tables/main_binary_summary_N5000_N10000.tsv",
  "results/supp/tables/supp1_enrichment_nns_nnr.tsv",
  "results/supp/tables/supp2_required_n_by_cohort.tsv",
  file.path("results", "supp", "tables", paste0("supp3_binary_vs_ordinal_metrics_", acc_label, tag_suffix, ".tsv")),
  file.path("results", "paper", "tables", "table1_main_parameters.tsv"),
  file.path("results", "paper", "tables", "table2_main_binary_operating_characteristics.tsv"),
  file.path("results", "paper", "tables", "tableS4_ordinal_summary_key_n.tsv")
)

check_outputs(expected_outputs)

manifest <- tibble::tibble(
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  quick_mode = quick_mode,
  accuracy = run_acc,
  path = expected_outputs
)

manifest_path <- file.path("results", "paper", "tables", "core_paper_manifest.tsv")
readr::write_tsv(manifest, manifest_path)
message(sprintf("[core] Complete. Manifest: %s", manifest_path))
