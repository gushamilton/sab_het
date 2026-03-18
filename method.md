# Methods (R1)

## High-Level Methods (Main Manuscript)

We used a simulation framework to evaluate how imperfect subphenotype classification affects detection of heterogeneous treatment effects (HTE) in SAB.

The analysis was intentionally simplified into one main question and two supplementary design questions.

1. Main analysis:
   We simulated randomized trial data with subgroup-specific treatment effects, then imposed imperfect subphenotype classification, and estimated subgroup effects using observed (potentially misclassified) subphenotypes.
   The primary outputs were Type I error (null subgroup), power (correct-direction detection), and Type M error (magnitude exaggeration among statistically significant findings).

2. Supplementary analysis 1 (enrichment feasibility):
   We evaluated enrichment designs under varying sensitivity/specificity to quantify effect dilution and screening burden (number needed to randomize and number needed to screen).

3. Supplementary analysis 2 (ordinal outcomes):
   We compared binary versus ordinal endpoints under the same subgroup-effect structure and misclassification settings to assess gains in precision and their relationship to bias.

Overall, the design tests two core manuscript points: accurate subphenotyping is essential for credible subgroup inference, and outcome design (especially ordinal outcomes) strongly affects inferential efficiency.

---

## Supplementary Methods (Technical Detail)

### 1. Simulation Inputs and Parameterization

Subphenotypes were represented as groups A-E with fixed:
- prevalence,
- baseline event risk,
- subgroup-specific treatment odds ratios.

The main scenario used shrunk ARREST-based subgroup ORs (`or_arrest_shrunk`), with group A as the null subgroup.

Core parameter definitions are in:
- `code/01_parameters.R`

### 2. Data-Generating Mechanisms

#### 2.1 Binary outcome generator

For each replicate:
- participants were assigned to subphenotypes by multinomial sampling,
- treatment was randomized 1:1,
- subgroup-specific treated risk was derived from baseline risk and subgroup OR,
- outcomes were sampled from Bernoulli distributions.

Implemented in:
- `simulate_binary_trial()` in `code/functions/01_simulation_helpers.R`

#### 2.2 Misclassification mechanism

Observed subphenotype labels were generated from true labels using an accuracy parameter.
When misclassification occurred, participants were reassigned to another subgroup with probability proportional to subgroup prevalence.

Implemented in:
- `misclassify_groups()` in `code/functions/01_simulation_helpers.R`

#### 2.3 Ordinal outcome generator

Ordinal outcomes were simulated under a proportional-odds structure with fixed baseline category probabilities and subgroup-specific treatment log-OR shifts.
This preserves a common ordinal scale while allowing subgroup-specific treatment effects.

Implemented in:
- `simulate_ordinal_trial()` in `code/functions/01_simulation_helpers.R`

### 3. Estimands and Performance Metrics

Subgroup treatment effects were estimated on the log-OR scale:
- binary models: logistic regression,
- ordinal models: proportional-odds logistic regression.

Metrics were computed by subgroup across replicates:
- Type I error: probability of statistical significance under the null subgroup,
- Power: probability of significance in the correct direction,
- Type S error: probability of significance in the wrong direction,
- Type M error: average absolute exaggeration ratio among significant estimates.

Implemented in:
- `fit_logistic_or()`, `fit_polr_or()` in `code/functions/01_simulation_helpers.R`
- `summarise_metrics()` in `code/functions/02_metrics.R`

### 4. R1 Analytic Components

#### 4.1 Main binary analysis (R1 primary)

Grid-based simulation over:
- sample size,
- classification accuracy.

Outputs include:
- replicate-level raw estimates,
- metric summaries,
- core power and Type S plots.

Script:
- `code/03_run_main_binary.R`

#### 4.2 Enrichment feasibility supplement

Closed-form calculations were used to:
- form expected enriched-cohort mixtures under test sensitivity/specificity,
- estimate diluted observed treatment effects,
- derive required sample size and number needed to screen.

Script:
- `code/04_run_enrichment.R`

#### 4.3 Cohort prevalence/event-rate variation supplement

Cohort-specific prevalence and mortality inputs were used to:
- recalculate required N,
- simulate fixed-N subgroup performance,
- optionally save replicate-level raw data for post hoc analyses.

Script:
- `code/05_run_cohort_variation.R`

#### 4.4 Ordinal versus binary supplement

Binary and ordinal scenarios were run on a dense sample-size grid.
Additional outputs include:
- explicit Type I comparison (binary vs ordinal),
- mean signed and absolute bias summaries,
- power-gain versus bias-change diagnostics,
- expected ordinal category-shift plots.

Replicate-level raw data are written for downstream custom plotting.

Script:
- `code/06_run_ordinal_comparison.R`

### 5. Tabulation and Presentation Outputs

TSV outputs are rendered into HTML summary tables using `gt`.

Script:
- `code/07_build_gt_tables.R`

### 6. Execution and Reproducibility

Local/smoke execution wrappers:
- `scripts/00_smoke_test.sh`
- `scripts/01_run_main.sh`
- `scripts/02_run_supp.sh`

Cluster execution wrappers:
- `slurm/01_main_binary.sh`
- `slurm/02_supplements.sh`
- `slurm/03_ordinal_acc100.sh`
- `slurm/04_r1_full_highres.sh` (high-rep, dense-grid, raw-output run)

To avoid provenance confusion between runs, pre-run outputs can be archived by cutoff timestamp:
- `scripts/03_archive_pre_run_outputs.sh`

### 7. Statistical Framing

The current implementation is frequentist and reports operating characteristics (Type I, power, Type S, Type M, bias summaries).
Bayesian design framing is not required for the current R1 implementation and can be layered later if needed.
