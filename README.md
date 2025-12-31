# SAB HTE simulation study

This repo runs simulation and closed-form analyses for heterogenous treatment effects
in Staphylococcus aureus bacteremia (SAB) subphenotypes.

## Key choices (current)
- Baseline risk: overall mortality across both arms (p0_overall).
- Subphenotype prevalences: overall counts across both arms.
- Main scenario: ARREST raw ORs.
- Supplementary scenarios: ARREST ORs with log-OR shrinkage k=0.5.
- Aim 2: fixed N = 10,000, accuracy grid 0.70 to 1.00.
- Aim 3: sens/spec grid in `code/common_parameters.R`, adaptive NNS search capped at 1,000,000.

## Structure
- `code/` main scripts and shared parameters
- `code/functions/` core simulation utilities
- `scripts/` runnable helpers (smoke test, run-all)
- `results/main/` primary outputs (ARREST_shrunk)
- `results/supp/` supplementary outputs (ARREST_raw + Conservative)
- `paper/` manuscript QMD and rendered HTML

## How to run
Smoke test (fast, main only):
```
./scripts/00_smoke_test.sh
```

Main analyses (Aim 1/2/3 + closed-form + combine):
```
./scripts/01_run_all_main.sh
```

Supplementary analyses:
```
./scripts/02_run_all_supp.sh
```

You can also run individual aims:
```
SCENARIO_SET=main Rscript code/03_run_aim1.R
SCENARIO_SET=main Rscript code/04_run_aim2.R
SCENARIO_SET=main Rscript code/05_run_aim3.R
SCENARIO_SET=main Rscript code/02_closed_form_aim3.R
SCENARIO_SET=main Rscript code/06_create_aim3_data.R
```

### Useful env vars
- `SCENARIO_SET`: `main` or `supp` (default: `main`)
- `N_REPS_AIM1_2`: Aim 1/2 reps (default: 1000)
- `AIM1_MODE`: `closed_form` (default) or `sim`
- `AIM1_VALIDATE`: set to `1` to run a small sim check for Aim 1
- `AIM1_VALIDATE_SIZES`: comma list for validation N (default: `1000,5000,10000,20000`)
- `N_REPS_AIM1_VALIDATE`: validation reps (default: 200)
- `N_REPS_AIM2`: Aim 2 reps (default: 1000)
- `N_REPS_AIM3`: Aim 3 reps (default: 200)
- `TRUE_OR_N`: sample size for Aim 2 overall OR calibration (default: 10,000,000)
- `AIM3_TEST_MODE=1`: limit Aim 3 to a single scenario for smoke tests
- `PRINT_GT=0`: suppress GT output in closed-form script

## Outputs (key files)
- Aim 1: `results/*/tables/aim1_power_summary.tsv`, `results/*/plots/aim1_power_vs_samplesize.pdf`
- Aim 2: `results/*/tables/aim2_accuracy_summary.tsv`, `results/*/plots/aim2_power_vs_accuracy.pdf`
- Aim 3 sim: `results/*/tables/aim3_sens_spec_summary.tsv`, `results/*/plots/aim3_nns_summary.pdf`
- Aim 3 closed form: `results/*/tables/aim3_closed_form_summary.tsv`
- Aim 3 combined: `results/*/tables/aim3_combined_summary.tsv`

## Notes
- Aim 3 simulations can take the longest; check logs in `results/*/logs/`.
- The manuscript is rendered from `paper/sab_hte_paper_draft_1.qmd`.
