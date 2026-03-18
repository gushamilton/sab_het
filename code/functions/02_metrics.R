## =====================================================================
##  02_metrics.R -- Metric calculations for power and errors
## =====================================================================

library(dplyr)

add_significance <- function(results, alpha) {
  z <- qnorm(1 - alpha / 2)
  results %>%
    mutate(
      ci_lower = log_or_hat - z * se,
      ci_upper = log_or_hat + z * se,
      significant = !is.na(log_or_hat) & !is.na(se) & (ci_lower > 0 | ci_upper < 0)
    )
}

summarise_metrics <- function(results, log_or_true, alpha, is_null = FALSE) {
  enriched <- add_significance(results, alpha)

  if (is_null) {
    type1 <- mean(enriched$significant, na.rm = TRUE)
    return(tibble(
      type1 = type1,
      power = NA_real_,
      type_s = NA_real_,
      type_m = NA_real_
    ))
  }

  sign_true <- sign(log_or_true)
  sign_hat <- sign(enriched$log_or_hat)

  power <- mean(enriched$significant & (sign_hat == sign_true), na.rm = TRUE)
  type_s <- mean(enriched$significant & (sign_hat != sign_true), na.rm = TRUE)
  type_m <- enriched %>%
    filter(significant) %>%
    summarise(ratio = mean(abs(log_or_hat) / abs(log_or_true), na.rm = TRUE)) %>%
    pull(ratio)

  tibble(
    type1 = NA_real_,
    power = power,
    type_s = type_s,
    type_m = type_m
  )
}
