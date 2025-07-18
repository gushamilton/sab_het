# Derivation of the Closed-Form (Analytic) Solution for Aim 3

This document outlines the mathematical derivation for the "closed-form" or analytic solution used to calculate the Number Needed to Screen (NNS) and Number Needed to Randomize (NNR) for an enrichment trial, as presented in Aim 3 of the study.

The core idea is to determine the properties of the final *enrolled cohort* (which is a mix of true-positive and false-positive patients) and then use a standard sample size formula to find the NNR required to achieve a desired statistical power for the *observed (diluted)* treatment effect in that cohort. The NNS is then derived from the NNR and the overall enrollment rate.

## 1. Defining Inputs

We start with the following known parameters:

-   **Subgroup Parameters:**
    -   `p_j`: The population prevalence of each subgroup `j`.
    -   `p0_j`: The baseline event rate (risk in the control arm) for each subgroup `j`.
    -   `OR_j`: The true odds ratio for the treatment effect in each subgroup `j`.
-   **Target and Test Parameters:**
    -   `target_group`: The specific subgroup we are enriching for (e.g., "B").
    -   `sensitivity (sens)`: The probability that the test correctly identifies a patient from the `target_group`.
    -   `specificity (spec)`: The probability that the test correctly identifies a patient who is *not* from the `target_group`.
-   **Statistical Parameters:**
    -   `alpha`: The desired Type I error rate (e.g., 0.05).
    -   `power`: The desired statistical power (e.g., 0.80).

## 2. Characterizing the Enrolled Cohort

The cohort enrolled in the trial is a mixture of two types of patients:
1.  **True Positives (TP):** Patients who are truly in the `target_group` and test positive.
2.  **False Positives (FP):** Patients who are *not* in the `target_group` but test positive.

Our first step is to calculate the properties of this mixed cohort.

### Step 2a: Properties of the "Contaminant Mix"

First, we define a conceptual "contaminant mix" which represents the average properties of all non-target subgroups.

-   **Prevalence of the target group:** `P(Target) = p_target`
-   **Prevalence of the non-target (contaminant) pool:** `P(Non-Target) = 1 - p_target`

The baseline event rate (`p0_mix`) and treated event rate (`p1_mix`) of this contaminant pool are the weighted averages of the rates of the individual non-target subgroups, where the weights are their prevalences *within the non-target pool*.

-   `p0_mix = sum(p_j * p0_j) / P(Non-Target)` for all `j != target_group`
-   `p1_mix = sum(p_j * p1_j) / P(Non-Target)` for all `j != target_group`, where `p1_j` is calculated from `p0_j` and `OR_j`.

### Step 2b: Composition of the Enrolled Cohort

Next, we determine the proportion of the total screened population that will be enrolled (the enrollment rate) and the composition of that enrolled group.

-   **True Positive Rate (as a fraction of total population):** `P(Test+|Target) * P(Target) = sens * p_target`
-   **False Positive Rate (as a fraction of total population):** `P(Test+|Non-Target) * P(Non-Target) = (1 - spec) * (1 - p_target)`

The **overall enrollment rate** is the sum of these two probabilities:
-   `enrol_rate = (sens * p_target) + ((1 - spec) * (1 - p_target))`

Within the enrolled cohort, the proportion of patients who are true targets (`p_true_in_enrolled`) is:
-   `p_true_in_enrolled = (sens * p_target) / enrol_rate`

The proportion who are contaminants is simply `1 - p_true_in_enrolled`.

### Step 2c: Observed (Diluted) Event Rates and Effect Size

The event rates we expect to *observe* in our trial (`p0_obs` and `p1_obs`) are the weighted averages of the rates from the true-positive and false-positive components of the enrolled cohort.

-   `p0_obs = (p_true_in_enrolled * p0_target) + ((1 - p_true_in_enrolled) * p0_mix)`
-   `p1_obs = (p_true_in_enrolled * p1_target) + ((1 - p_true_in_enrolled) * p1_mix)`

From these observed rates, we can calculate the **observed (diluted) log-odds ratio (`beta_obs`)**:
-   `beta_obs = log( (p1_obs / (1 - p1_obs)) / (p0_obs / (1 - p0_obs)) )`

This `beta_obs` is the effect size that our enrichment trial is actually powered to detect. It is "diluted" from the true `beta_target` due to the inclusion of contaminant patients.

## 3. Sample Size Calculation

With the parameters of the enrolled cohort now defined, we can use a standard sample size formula for a two-group comparison of proportions (based on the log-odds ratio) with 1:1 randomization.

### Step 3a: Number Needed to Randomize (NNR)

The formula for the total sample size (`N_total`, which is our NNR) is:

`NNR = 2 * (Z_alpha/2 + Z_power)^2 * Var(beta_obs) / (beta_obs)^2`

Where:
-   `Z_alpha/2` is the Z-score for the two-sided alpha level (e.g., `qnorm(1 - 0.05/2) = 1.96`).
-   `Z_power` is the Z-score for the desired power (e.g., `qnorm(0.80) = 0.84`).
-   `Var(beta_obs)` is the variance of the observed log-odds ratio, approximated by:
    -   `Var(beta_obs) = 1 / (p0_obs * (1 - p0_obs)) + 1 / (p1_obs * (1 - p1_obs))`

The initial factor of `2` in the NNR formula accounts for the two arms of the trial (treatment and control).

### Step 3b: Number Needed to Screen (NNS)

The NNS is the total number of patients we must screen to get the required number of patients for randomization (NNR). This is simply the NNR divided by the overall enrollment rate.

-   `NNS = ceil(NNR / enrol_rate)`

We take the ceiling to ensure we have a whole number of patients.

This completes the derivation of the analytic solution used to provide a rapid, simulation-free estimate of the required sample sizes for the enrichment trial scenarios. 