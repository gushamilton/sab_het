## =====================================================================
##  15_run_ordinal_nonPO_comparison.R
##  Compare ordinal (polr) and binary (death cut-point) under
##  proportional-odds (PO) and non-proportional-odds (nonPO) DGMs.
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")
source("code/functions/01_simulation_helpers.R")
source("code/functions/02_metrics.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table

or_vector <- setNames(subphenotypes$or_arrest_shrunk, subphenotypes$subphenotype)
prevalence <- setNames(subphenotypes$prevalence, subphenotypes$subphenotype)
baseline <- setNames(subphenotypes$baseline_mortality, subphenotypes$subphenotype)
ordinal_points <- as.integer(Sys.getenv("ORDINAL_POINTS", unset = "6"))
if (ordinal_points == 5) {
  ordinal_baseline <- params$ordinal_baseline$prevalence
  survivor_split <- derive_survivor_split(ordinal_baseline)
} else if (ordinal_points == 6) {
  # Legacy paper-facing 6-point setup: death plus five non-death levels.
  survivor_split <- c(0.08, 0.12, 0.20, 0.25, 0.35)
} else {
  stop("ORDINAL_POINTS must be 5 or 6.", call. = FALSE)
}
ordinal_levels <- seq_len(ordinal_points)

sample_sizes <- params$sample_sizes
sample_override <- Sys.getenv("SAMPLE_SIZES", unset = "")
if (nzchar(sample_override)) {
  sample_sizes <- as.integer(strsplit(sample_override, ",")[[1]])
}

accuracy_grid <- c(0.70, 1.00)
n_reps <- as.integer(Sys.getenv("N_REPS_DOOR", unset = "1000"))
seed_base <- as.integer(Sys.getenv("SEED_BASE", unset = "15001"))
rep_offset <- as.integer(Sys.getenv("REP_OFFSET", unset = "0"))
alpha <- params$alpha_primary
calibration_n <- as.integer(Sys.getenv("CALIBRATION_N", unset = "200000"))
skip_cal <- tolower(Sys.getenv("SKIP_CALIBRATION", unset = "0")) %in% c("1", "true", "yes")
calibration_only <- tolower(Sys.getenv("CALIBRATION_ONLY", unset = "0")) %in% c("1", "true", "yes")
run_tag <- Sys.getenv("RUN_TAG", unset = "")
tag_suffix <- if (nzchar(run_tag)) paste0("_", run_tag) else ""

stopifnot(abs(sum(survivor_split) - 1) < 1e-8)

output_root <- "results/supp/ordinal_nonPO_comparison"
for (d in c("data", "tables", "plots")) {
  dir.create(file.path(output_root, d), showWarnings = FALSE, recursive = TRUE)
}

build_cum_control <- function(p_dead, survivor_split) {
  nondeath <- cumsum((1 - p_dead) * survivor_split)
  c(p_dead, p_dead + head(nondeath, -1))
}

simulate_po_ordinal_trial <- function(or_vector, prevalence, baseline, n, seed) {
  set.seed(seed)
  groups <- names(or_vector)

  dat <- tibble(
    id = seq_len(n),
    group = sample(groups, size = n, replace = TRUE, prob = prevalence),
    treatment = rbinom(n, 1, 0.5)
  )

  out <- map2_int(dat$group, dat$treatment, function(g, t) {
    p0 <- baseline[[g]]
    beta <- log(or_vector[[g]])
    cum_ctrl <- build_cum_control(p0, survivor_split)
    theta <- qlogis(cum_ctrl)
    cum_t <- plogis(theta + beta * t)
    probs <- diff(c(0, cum_t, 1))
    probs <- pmax(probs, 0)
    probs <- probs / sum(probs)
    sample(ordinal_levels, size = 1, prob = probs)
  })

  dat %>% mutate(outcome_ord = out)
}

rescaled_nonpo_probs <- function(p_dead, survivor_split, or_death) {
  cum_ctrl <- build_cum_control(p_dead, survivor_split)
  p_ctrl <- diff(c(0, cum_ctrl, 1))

  odds0 <- p_dead / (1 - p_dead)
  odds1 <- odds0 * or_death
  p_dead_t <- odds1 / (1 + odds1)

  surv_w <- p_ctrl[-1] / sum(p_ctrl[-1])
  p_trt <- c(p_dead_t, (1 - p_dead_t) * surv_w)

  list(control = p_ctrl, treatment = p_trt)
}

simulate_nonPO_trial <- function(or_vector, prevalence, baseline, n, seed) {
  # Non-PO (rescaled survivors): treatment changes death odds only; conditional
  # distribution across non-death levels is held fixed.
  set.seed(seed)
  groups <- names(or_vector)

  dat <- tibble(
    id = seq_len(n),
    group = sample(groups, size = n, replace = TRUE, prob = prevalence),
    treatment = rbinom(n, 1, 0.5)
  )

  out <- map2_int(dat$group, dat$treatment, function(g, t) {
    p0 <- baseline[[g]]
    d <- rescaled_nonpo_probs(p0, survivor_split, or_vector[[g]])
    probs <- if (t == 1) d$treatment else d$control
    sample(ordinal_levels, size = 1, prob = probs)
  })

  dat %>% mutate(outcome_ord = out)
}

probs_from_po_params <- function(theta, beta, trt) {
  cum <- plogis(theta + beta * trt)
  p <- diff(c(0, cum, 1))
  pmax(p, 1e-12)
}

compute_nonPO_polr_truth <- function(p_dead, survivor_split, or_death) {
  d <- rescaled_nonpo_probs(p_dead, survivor_split, or_death)
  p0 <- d$control
  p1 <- d$treatment

  nll <- function(par) {
    t1 <- par[1]
    inc <- exp(par[2:(length(p0) - 1)])
    theta <- c(t1, t1 + cumsum(inc))
    beta <- par[length(p0)]

    q0 <- probs_from_po_params(theta, beta, trt = 0)
    q1 <- probs_from_po_params(theta, beta, trt = 1)

    -(sum(p0 * log(q0)) + sum(p1 * log(q1)))
  }

  init_theta <- qlogis(cumsum(p0)[1:(length(p0) - 1)])
  init_par <- c(
    init_theta[1],
    log(pmax(diff(init_theta), 1e-3)),
    log(or_death)
  )

  fit <- optim(init_par, nll, method = "BFGS", control = list(maxit = 2000))
  as.numeric(fit$par[length(p0)])
}

true_log_or_polr_nonPO <- setNames(
  map_dbl(names(or_vector), function(g) {
    compute_nonPO_polr_truth(baseline[[g]], survivor_split, or_vector[[g]])
  }),
  names(or_vector)
)

# ── DOOR analysis: retained for reference but excluded from main simulation ──
#
# Under proportional odds (PO), DOOR (win odds) and polr (proportional odds
# logistic regression) are asymptotically equivalent: both extract the same
# rank-based information from the data, and empirical simulations confirmed
# near-identical proportional bias and power under both perfect and imperfect
# (70%) classification. DOOR is therefore redundant given polr under our DGM.
#
# DOOR would be preferred over polr only under non-PO with opposing treatment
# effects across the outcome scale (e.g. mortality benefit but functional harm),
# as polr would then be misspecified while DOOR remains consistent for its
# estimand. We do not simulate this scenario as it requires assumed effect
# magnitudes at non-death outcome levels for which no empirical data exist.
#
# Reference: Hamilton et al. simulation (internal), confirmed in plots:
#   door_comparison_dist_n20000_r50.pdf
#   door_comparison_bias_n20000_r50.pdf
#   door_comparison_estimands_n20000_r50.pdf
#
# compute_true_log_win_odds <- function(p_dead, survivor_split, or_value) {
#   cum_ctrl <- build_cum_control(p_dead, survivor_split)
#   theta <- qlogis(cum_ctrl)
#   beta <- log(or_value)
#   cum_trt <- plogis(theta + beta)
#   p_ctrl <- diff(c(0, cum_ctrl, 1))
#   p_trt <- diff(c(0, cum_trt, 1))
#   k <- length(p_ctrl)
#   idx <- seq_len(k)
#   n_win <- sum(outer(idx, idx, function(t, c) (t > c) * p_trt[t] * p_ctrl[c]))
#   n_loss <- sum(outer(idx, idx, function(t, c) (t < c) * p_trt[t] * p_ctrl[c]))
#   log(n_loss / n_win)
# }
#
# fit_door_win_odds <- function(data) {
#   trt <- data %>% filter(treatment == 1) %>% pull(outcome)
#   ctrl <- data %>% filter(treatment == 0) %>% pull(outcome)
#   n1 <- length(trt); n0 <- length(ctrl)
#   if (n1 < 5 || n0 < 5) {
#     return(tibble(log_or_hat = NA_real_, se = NA_real_))
#   }
#   wins <- sum(outer(trt, ctrl, `>`))
#   losses <- sum(outer(trt, ctrl, `<`))
#   if (wins == 0 || losses == 0) {
#     return(tibble(log_or_hat = NA_real_, se = NA_real_))
#   }
#   log_wo <- log(losses / wins)
#   se <- sqrt(1 / wins + 1 / losses)
#   tibble(log_or_hat = log_wo, se = se)
# }

fit_all_models_by_group <- function(sim_data) {
  groups <- sort(unique(sim_data$observed_group))

  map_dfr(groups, function(g) {
    d <- sim_data %>% filter(observed_group == g)

    polr_fit <- fit_polr_or(
      d %>% transmute(treatment, outcome = outcome_ord)
    ) %>% mutate(model_type = "Ordinal (polr)")

    bin_fit <- fit_logistic_or(
      d %>% transmute(treatment, outcome = as.integer(outcome_ord == 1))
    ) %>% mutate(model_type = "Binary (death)")

    bind_rows(polr_fit, bin_fit) %>% mutate(group = g)
  })
}

run_scenario <- function(n, accuracy, dgm, n_reps, seed_base, rep_offset = 0L) {
  sim_fn <- if (dgm == "PO") {
    function(...) simulate_po_ordinal_trial(...)
  } else {
    function(...) simulate_nonPO_trial(...)
  }

  map_dfr(seq_len(n_reps), function(k) {
    rep_id <- rep_offset + k

    sim <- sim_fn(
      or_vector = or_vector,
      prevalence = prevalence,
      baseline = baseline,
      n = n,
      seed = seed_base + rep_id
    )

    sim$observed_group <- misclassify_groups(
      true_group = sim$group,
      prevalence = prevalence,
      accuracy = accuracy,
      seed = seed_base + 10000L + rep_id
    )

    fit_all_models_by_group(sim) %>%
      mutate(rep_id = rep_id, sample_size = n, accuracy = accuracy, dgm = dgm)
  })
}

calibrate_once <- function(n, dgm, seed = 999) {
  sim_fn <- if (dgm == "PO") {
    function(...) simulate_po_ordinal_trial(...)
  } else {
    function(...) simulate_nonPO_trial(...)
  }

  sim <- sim_fn(or_vector = or_vector, prevalence = prevalence, baseline = baseline, n = n, seed = seed)
  sim$observed_group <- sim$group

  fits <- fit_all_models_by_group(sim)

  fits %>%
    mutate(
      dgm = dgm,
      log_or_true = case_when(
        dgm == "PO" & model_type == "Ordinal (polr)" ~ log(or_vector[group]),
        dgm == "PO" & model_type == "Binary (death)" ~ log(or_vector[group]),
        dgm == "nonPO" & model_type == "Binary (death)" ~ log(or_vector[group]),
        dgm == "nonPO" & model_type == "Ordinal (polr)" ~ true_log_or_polr_nonPO[group]
      ),
      raw_bias = log_or_hat - log_or_true,
      prop_bias = raw_bias / abs(log_or_true)
    ) %>%
    select(dgm, group, model_type, log_or_true, log_or_hat, raw_bias, prop_bias)
}

if (!skip_cal) {
  calibration <- bind_rows(
    calibrate_once(calibration_n, "PO"),
    calibrate_once(calibration_n, "nonPO")
  )
  write_tsv(
    calibration,
    file.path(output_root, "tables", paste0("ordinal_nonPO_comparison_calibration", tag_suffix, ".tsv"))
  )
}

if (calibration_only) {
  message("Calibration-only mode complete.")
  quit(save = "no", status = 0)
}

dgm_grid <- expand_grid(
  sample_size = sample_sizes,
  accuracy = accuracy_grid,
  dgm = c("PO", "nonPO")
)

message("Running simulation grid...")

raw_results <- pmap_dfr(dgm_grid, function(sample_size, accuracy, dgm) {
  message(sprintf("  DGM = %s, N = %d, accuracy = %.2f", dgm, sample_size, accuracy))
  run_scenario(sample_size, accuracy, dgm, n_reps, seed_base, rep_offset = rep_offset)
}) %>%
  mutate(
    log_or_true = case_when(
      dgm == "PO" & model_type == "Ordinal (polr)" ~ log(or_vector[group]),
      dgm == "PO" & model_type == "Binary (death)" ~ log(or_vector[group]),
      dgm == "nonPO" & model_type == "Binary (death)" ~ log(or_vector[group]),
      dgm == "nonPO" & model_type == "Ordinal (polr)" ~ true_log_or_polr_nonPO[group]
    )
  )

write_tsv(
  raw_results,
  file.path(output_root, "data", paste0("ordinal_nonPO_comparison_raw", tag_suffix, ".tsv.gz"))
)

metric_results <- raw_results %>%
  group_by(sample_size, accuracy, dgm, group, model_type) %>%
  group_modify(~ {
    truth <- unique(.x$log_or_true)
    is_null <- abs(truth) < 1e-6
    summarise_metrics(.x, log_or_true = truth, alpha = alpha, is_null = is_null)
  }) %>%
  ungroup()

write_tsv(
  metric_results,
  file.path(output_root, "tables", paste0("ordinal_nonPO_comparison_metrics", tag_suffix, ".tsv"))
)

bias_results <- raw_results %>%
  filter(group != "A") %>%
  mutate(
    raw_bias = log_or_hat - log_or_true,
    prop_bias = raw_bias / abs(log_or_true)
  ) %>%
  group_by(sample_size, accuracy, dgm, group, model_type) %>%
  summarise(
    median_prop_bias = median(prop_bias, na.rm = TRUE),
    mean_prop_bias = mean(prop_bias, na.rm = TRUE),
    rmse_log = sqrt(mean(raw_bias^2, na.rm = TRUE)),
    .groups = "drop"
  )

write_tsv(
  bias_results,
  file.path(output_root, "tables", paste0("ordinal_nonPO_comparison_bias", tag_suffix, ".tsv"))
)

model_colours <- c(
  "Ordinal (polr)" = "#2166ac",
  "Binary (death)" = "#d6604d"
)
dgm_linetype <- c("PO" = "solid", "nonPO" = "dashed")

power_data <- metric_results %>%
  filter(group != "A") %>%
  mutate(
    panel = factor(
      paste0(if_else(accuracy == 0.7, "Accuracy = 70%", "Accuracy = 100%"), "\nDGM = ", dgm),
      levels = c(
        "Accuracy = 70%\nDGM = PO",
        "Accuracy = 70%\nDGM = nonPO",
        "Accuracy = 100%\nDGM = PO",
        "Accuracy = 100%\nDGM = nonPO"
      )
    )
  )

p_power <- power_data %>%
  ggplot(aes(x = sample_size, y = power, colour = model_type, linetype = dgm)) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(group ~ panel, scales = "free_y") +
  scale_colour_manual(values = model_colours) +
  scale_linetype_manual(values = dgm_linetype) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title = "Power under proportional and non-proportional odds DGMs",
    x = "Total trial size",
    y = "Power",
    colour = "Model",
    linetype = "DGM"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_root, "plots", paste0("ordinal_nonPO_comparison_power", tag_suffix, ".pdf")),
  p_power, width = 14, height = 9
)

dist_data <- raw_results %>%
  filter(sample_size == max(sample_sizes), group != "A") %>%
  mutate(
    prop_bias = (log_or_hat - log_or_true) / abs(log_or_true),
    panel = factor(
      paste0(if_else(accuracy == 0.7, "Accuracy = 70%", "Accuracy = 100%"), "\nDGM = ", dgm),
      levels = c(
        "Accuracy = 70%\nDGM = PO",
        "Accuracy = 70%\nDGM = nonPO",
        "Accuracy = 100%\nDGM = PO",
        "Accuracy = 100%\nDGM = nonPO"
      )
    )
  )

p_dist <- dist_data %>%
  ggplot(aes(x = model_type, y = prop_bias, fill = model_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_grid(group ~ panel, scales = "free_y") +
  scale_fill_manual(values = model_colours) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = paste0("Proportional bias distributions at N = ", max(sample_sizes)),
    x = NULL,
    y = "Proportional bias",
    fill = "Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_root, "plots", paste0("ordinal_nonPO_comparison_dist", tag_suffix, ".pdf")),
  p_dist, width = 14, height = 9
)

estimand_tbl <- expand_grid(
  group = names(or_vector),
  dgm = c("PO", "nonPO"),
  model_type = c("Binary (death)", "Ordinal (polr)")
) %>%
  mutate(
    true_log_effect = case_when(
      dgm == "PO" & model_type == "Binary (death)" ~ log(or_vector[group]),
      dgm == "PO" & model_type == "Ordinal (polr)" ~ log(or_vector[group]),
      dgm == "nonPO" & model_type == "Binary (death)" ~ log(or_vector[group]),
      dgm == "nonPO" & model_type == "Ordinal (polr)" ~ true_log_or_polr_nonPO[group]
    ),
    estimand = paste(model_type, dgm, sep = " / ")
  )

p_estimand <- estimand_tbl %>%
  ggplot(aes(x = group, y = true_log_effect, fill = interaction(model_type, dgm, sep = " / "))) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, colour = "grey30") +
  scale_fill_manual(values = c(
    "Binary (death) / PO" = "#d6604d",
    "Ordinal (polr) / PO" = "#2166ac",
    "Binary (death) / nonPO" = scales::muted("#d6604d"),
    "Ordinal (polr) / nonPO" = scales::muted("#2166ac")
  )) +
  labs(
    title = "True estimands under PO vs nonPO",
    x = "Subphenotype",
    y = "True log-effect",
    fill = "Method / DGM"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_root, "plots", paste0("ordinal_nonPO_comparison_estimands", tag_suffix, ".pdf")),
  p_estimand, width = 10, height = 5
)

message("Done. Results in: ", output_root)
