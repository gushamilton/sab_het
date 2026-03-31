# Supplementary Methods

## ADEMP Framework

### Aims

We used a simulation study to evaluate how imperfect subphenotype classification influences detection of heterogeneous treatment effects (HTE) in *Staphylococcus aureus* bacteraemia (SAB).

The primary aim was to quantify effects of misclassification on subgroup-specific inference for binary mortality outcomes. Prespecified supplementary aims were to evaluate: (i) sensitivity of sample-size requirements to cohort-specific prevalence and baseline risk, (ii) feasibility of enrichment strategies under imperfect classifiers, and (iii) performance of binary versus ordinal analyses under proportional-odds and non-proportional data-generating mechanisms.

### Data-Generating Mechanisms

#### Population structure and treatment effects

Participants were assigned to one of five subphenotypes (A-E) according to specified prevalences. Subphenotype-specific baseline mortality risks and treatment effects were parameterized from published SAB subphenotype analyses and ARREST-derived estimates. To reduce winner's-curse inflation, subgroup treatment effects were conservatively shrunk by 50% on the log-odds scale in primary analyses.

Within each replicate, treatment was randomized 1:1.

#### Misclassification mechanism

Observed subphenotype labels were generated from true labels using a scalar classification-accuracy parameter. With probability equal to accuracy, labels were correct; otherwise individuals were reassigned to one of the remaining subphenotypes with probability proportional to subphenotype prevalence.

#### Binary-outcome mechanism

For binary analyses, death was generated from subgroup-specific baseline risk and treatment effect on the log-odds scale.

#### Ordinal-outcome mechanisms

The publication ordinal analyses used a six-level outcome, ordered from worst to best:

1. Death
2. ICU/ventilated
3. Still hospitalised
4. Discharged to rehab
5. Discharged with complications
6. Discharged well

Two ordinal data-generating mechanisms were considered:

1. Proportional-odds (PO): treatment applies a common shift across all cumulative logits.
2. Non-proportional (death-only): treatment changes death odds only; the conditional distribution among non-death categories is held fixed and rescaled to sum to \(1 - p_{death,treatment}\).

### Estimands

The estimands were defined on the log-odds scale.

- Binary model estimand: subgroup-specific treatment effect for death.
- Ordinal model estimand under PO: subgroup-specific common proportional-odds effect.
- Ordinal model estimand under non-PO: pseudo-true proportional-odds coefficient under model misspecification.

### Methods

For each replicate and subgroup, treatment effects were estimated using:

- logistic regression for binary death outcomes,
- proportional-odds ordinal logistic regression for the six-level outcome.

For the core misclassification-accuracy analysis (Figure 1), total trial size was fixed at 20,000 participants while classification accuracy varied from 70% to 100%.

For the ordinal PO versus non-PO comparison, the main manuscript figures use a six-level ordinal outcome, with the main power comparison shown at n = 3,000 and the illustrative category-shift figure generated for a representative subgroup using the same six-level scale. The detailed numerical comparison is provided in Supplementary Table S8. A five-level version is retained as a sensitivity analysis and summarised in Supplementary Figure S1.

Power was defined using two-sided hypothesis testing with effect-direction consistency. Type I error was assessed in null-effect subgroups.

### Performance Measures

Primary performance measures were:

- Power: proportion of replicates with statistically significant effect in the true direction.
- Type I error: proportion of significant effects under a true null subgroup effect.
- Bias: difference between estimated and true subgroup log-effect.

Additional supplementary summaries included mean signed bias, mean absolute bias, proportional bias, and root-mean-square error on the log-odds scale, together with Type S and Type M error summaries where relevant.

## Supplementary Notes

All analyses were fully prespecified and script-generated, with figures and tables produced directly from simulation outputs to minimize transcription error. Publication-facing supplementary tables were formatted from the canonical simulation outputs with rounded values and explicit column labels for reporting clarity.
