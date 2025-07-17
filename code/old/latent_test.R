## ---- quick_latent_debug.R -----------------------------------------------
## Assumes draw_patient_batch(), compute_tau(), convert_or(), p0_raw,
##      or_vec are already in your workspace (they are in your script).

library(flexmix)

set.seed(1)

# ----- pick one scenario ---------------------------------------------------
target_group   <- "B"
snr            <- 3
n_to_randomize <- 10000           # small-ish so it runs fast
sens_target    <- 0.95

# ----- generate one pilot pool and threshold ------------------------------
pilot_pool <- draw_patient_batch(1e5, target_group, snr)
tau        <- compute_tau(pilot_pool, "sens", sens_target, NA_real_)

# ----- probabilistic recruitment (Design C/D/E) ---------------------------
new_batch <- draw_patient_batch(5e4, target_group, snr)
p_enroll  <- new_batch$p
cohort_prob <- new_batch[rbinom(nrow(new_batch), 1, p_enroll) == 1, ]
cohort_prob <- head(cohort_prob, n_to_randomize)

# outcome + treatment for mixture model
cohort_prob$treat <- rbinom(nrow(cohort_prob), 1, 0.5)
p0 <- p0_raw[cohort_prob$group]
p1 <- convert_or(or_vec[cohort_prob$group], p0)
cohort_prob$y <- rbinom(nrow(cohort_prob), 1,
                        ifelse(cohort_prob$treat == 1, p1, p0))


d <- as.data.frame(cohort_prob)

## binomial response must be 2‑column matrix
d$y_mat <- cbind(d$y, 1 - d$y)

# ----- fit 2-component mixture --------------------------------------------
mix_fit <- flexmix(
  y_mat ~ treat | p,          # <‑‑ use the 2‑col matrix here
  data  = d,
  k     = 2,
  model = FLXMRglm(family = "binomial")
)
if (inherits(mix_fit, "try-error")) {
  cat("\n*** flexmix failed — root cause below ***\n")
  print(mix_fit)
} else {
  print(mix_fit)
  cat("\nComponent-specific treatment log-ORs:\n")
  print(parameters(mix_fit)["treat", ])
}
## 1. turn the parameters slot into a 2‑column matrix
coef_mat <- sapply(mix_fit@components, function(cmp) cmp[[1]]@parameters$coef)

##    rownames: "(Intercept)" "treat"
##    colnames: "Comp.1" "Comp.2"

## 2. pick the component whose patients have higher mean biomarker score
comp_means  <- tapply(d$p, clusters(mix_fit), mean)
target_comp <- names(which.max(comp_means))

## 3. extract the treatment coefficient
beta_target <- coef_mat["treat", which.max(comp_means)]

beta_target    # this is the latent‑class es

