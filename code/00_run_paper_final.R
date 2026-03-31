## =====================================================================
##  00_run_paper_final.R -- Final manuscript production pipeline
## =====================================================================

library(tidyverse)

run_with_env <- function(script_path, env = list()) {
  old <- setNames(
    lapply(names(env), function(nm) Sys.getenv(nm, unset = NA_character_)),
    names(env)
  )
  on.exit({
    for (nm in names(env)) {
      if (is.na(old[[nm]][[1]])) {
        Sys.unsetenv(nm)
      } else {
        do.call(Sys.setenv, setNames(list(old[[nm]][[1]]), nm))
      }
    }
  }, add = TRUE)

  if (length(env) > 0) do.call(Sys.setenv, env)
  message(sprintf("[paper-final] Running %s", script_path))
  source(script_path, local = new.env(parent = globalenv()))
}

quick_mode <- tolower(Sys.getenv("PAPER_QUICK", unset = "0")) %in% c("1", "true", "yes")
acc <- as.numeric(Sys.getenv("PAPER_ACCURACY", unset = "1.00"))
nonpo_tag <- Sys.getenv("NONPO_RUN_TAG", unset = "ordinal6_figure5_20260325")
skip_nonpo <- tolower(Sys.getenv("PAPER_SKIP_NONPO", unset = "0")) %in% c("1", "true", "yes")
skip_assets <- tolower(Sys.getenv("PAPER_SKIP_ASSETS", unset = "0")) %in% c("1", "true", "yes")

if (quick_mode) message("[paper-final] Quick mode enabled.")

# Preserve all raw/intermediate outputs for downstream secondary work.
dir.create(file.path("results", "paper", "final"), showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

main_env <- list(
  N_REPS_MAIN = if (quick_mode) "50" else Sys.getenv("N_REPS_MAIN", unset = "2000"),
  SAMPLE_SIZES = if (quick_mode) "500,1000" else Sys.getenv("SAMPLE_SIZES", unset = "")
)

cohort_env <- list(
  N_REPS_COHORT = if (quick_mode) "50" else Sys.getenv("N_REPS_COHORT", unset = "1000"),
  COHORT_FIXED_N = if (quick_mode) "1000" else Sys.getenv("COHORT_FIXED_N", unset = "5000"),
  SAVE_COHORT_RAW = "1"
)

fig1_env <- list(
  FIG1_N_REPS = if (quick_mode) "50" else Sys.getenv("FIG1_N_REPS", unset = "1000"),
  FIG1_N = if (quick_mode) "1000" else Sys.getenv("FIG1_N", unset = "20000")
)

ordinal_env <- list(
  ACCURACY = sprintf("%.2f", acc),
  ORDINAL_POINTS = Sys.getenv("ORDINAL_POINTS", unset = "6"),
  N_REPS_ORDINAL = if (quick_mode) "50" else Sys.getenv("N_REPS_ORDINAL", unset = "1000"),
  ORDINAL_SAMPLE_SIZES = if (quick_mode) "500,1000,2000" else Sys.getenv("ORDINAL_SAMPLE_SIZES", unset = ""),
  RUN_ORDER = Sys.getenv("RUN_ORDER", unset = "binary,ordinal"),
  RUN_TAG = Sys.getenv("RUN_TAG", unset = "")
)

nonpo_env <- list(
  RUN_TAG = nonpo_tag,
  N_REPS_DOOR = if (quick_mode) "100" else Sys.getenv("N_REPS_NONPO", unset = "1000"),
  CALIBRATION_N = if (quick_mode) "50000" else Sys.getenv("CALIBRATION_N", unset = "200000"),
  SAMPLE_SIZES = if (quick_mode) "2500" else Sys.getenv("NONPO_SAMPLE_SIZES", unset = "")
)

run_with_env("code/03_run_main_binary.R", env = main_env)
run_with_env("code/04_run_enrichment.R")
run_with_env("code/05_run_cohort_variation.R", env = cohort_env)

run_with_env("code/11_run_figure1_accuracy_scan.R", env = fig1_env)
run_with_env("code/12_plot_figure2_master.R")
run_with_env("code/13_plot_figure3_cohort_summary.R")

run_with_env("code/06_run_ordinal_comparison.R", env = ordinal_env)
run_with_env("code/26_build_core_figure4_ordinal.R", env = list(ACCURACY = sprintf("%.2f", acc)))
run_with_env("code/30_build_paper_ordinal_shift_figure.R")

if (!skip_nonpo) {
  # PO vs non-PO analysis retained for manuscript Table 2 and supplementary detail.
  run_with_env("code/15_run_ordinal_nonPO_comparison.R", env = nonpo_env)
  run_with_env("code/22_summarise_relative_rmse_nonpo.R", env = list(RUN_TAG = nonpo_tag))

  # Optional manuscript-facing nonPO diagnostic plots.
  build_nonpo_plots <- tolower(Sys.getenv("PAPER_BUILD_NONPO_PLOTS", unset = "1")) %in% c("1", "true", "yes")
  if (build_nonpo_plots) {
    run_with_env("code/20_plot_ordinal_po_nonpo_patchwork.R", env = list(RUN_TAG = nonpo_tag))
    run_with_env("code/21_plot_ordinal_po_nonpo_centered_bias_patchwork.R", env = list(RUN_TAG = nonpo_tag))
    run_with_env("code/23_plot_split_bias_by_model_nonpo.R", env = list(RUN_TAG = nonpo_tag))
  }
} else {
  message("[paper-final] Skipping nonPO simulation and nonPO plot builds.")
}

if (!skip_assets) {
  run_with_env(
    "code/31_plot_ordinal_power_5pt_vs_6pt.R",
    env = list(
      RUN_TAG_5PT = Sys.getenv("RUN_TAG_5PT", unset = "paper_final_5pt_20260324"),
      RUN_TAG_6PT = nonpo_tag
    )
  )
  run_with_env(
    "code/28_build_paper_final_assets.R",
    env = list(
      PAPER_ACCURACY = sprintf("%.2f", acc),
      NONPO_RUN_TAG = nonpo_tag
    )
  )
} else {
  message("[paper-final] Skipping final manuscript asset build.")
}

message("[paper-final] Complete.")
