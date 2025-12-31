# Closed-form power utilities for Aim 1

get_or_p1 <- function(or, p0) {
  (or * p0) / (1 - p0 + (or * p0))
}

calc_power_closed_form <- function(n_total, p0, or, prevalence, alpha = 0.05 / 5) {
  n_sub <- n_total * prevalence
  n0 <- n_sub / 2
  n1 <- n_sub / 2
  p1 <- get_or_p1(or, p0)
  log_or <- log(or)
  var_term <- (1 / (n0 * p0 * (1 - p0))) + (1 / (n1 * p1 * (1 - p1)))
  z_alpha <- qnorm(1 - alpha / 2)
  z <- abs(log_or) / sqrt(var_term)
  1 - (pnorm(z_alpha - z) - pnorm(-z_alpha - z))
}

solve_n_required <- function(target_power, p0, or, prevalence, alpha = 0.05 / 5,
                             n_min = 500, n_max = 100000, max_iter = 40) {
  low <- n_min
  high <- n_max
  pow_low <- calc_power_closed_form(low, p0, or, prevalence, alpha)
  pow_high <- calc_power_closed_form(high, p0, or, prevalence, alpha)
  if (pow_low >= target_power) return(low)
  if (pow_high < target_power) return(NA_real_)
  for (i in 1:max_iter) {
    mid <- floor((low + high) / 2)
    if (mid == low || mid == high) break
    pow_mid <- calc_power_closed_form(mid, p0, or, prevalence, alpha)
    if (pow_mid < target_power) {
      low <- mid
    } else {
      high <- mid
    }
  }
  high
}
