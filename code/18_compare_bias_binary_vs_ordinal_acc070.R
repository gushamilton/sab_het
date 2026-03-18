## =====================================================================
##  18_compare_bias_binary_vs_ordinal_acc070.R
##  Formal paired comparison of Binary vs Ordinal bias at accuracy=0.70
## =====================================================================

library(tidyverse)

input_path <- "results/supp/ordinal6_sameeffect/data/ordinal6_sameeffect_raw_n2500_n20000.tsv.gz"
output_root <- "results/supp/ordinal6_sameeffect/tables"
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)

raw <- read_tsv(input_path, show_col_types = FALSE) %>%
  filter(accuracy == 0.7, group != "A") %>%
  mutate(
    signed_err = log_or_hat - log_or_true,
    abs_err = abs(signed_err)
  )

wide <- raw %>%
  select(sample_size, group, rep_id, model_type, signed_err, abs_err) %>%
  pivot_wider(
    names_from = model_type,
    values_from = c(signed_err, abs_err),
    names_sep = "_"
  ) %>%
  filter(!is.na(signed_err_Binary), !is.na(signed_err_Ordinal))

compare_one <- function(df) {
  d_signed <- df$signed_err_Ordinal - df$signed_err_Binary
  d_abs <- df$abs_err_Ordinal - df$abs_err_Binary

  wt_signed <- suppressWarnings(wilcox.test(d_signed, mu = 0, alternative = "two.sided", exact = FALSE))
  wt_abs <- suppressWarnings(wilcox.test(d_abs, mu = 0, alternative = "two.sided", exact = FALSE))

  tibble(
    n_pairs = nrow(df),
    median_signed_binary = median(df$signed_err_Binary, na.rm = TRUE),
    median_signed_ordinal = median(df$signed_err_Ordinal, na.rm = TRUE),
    median_delta_signed = median(d_signed, na.rm = TRUE),
    mean_delta_signed = mean(d_signed, na.rm = TRUE),
    sd_delta_signed = sd(d_signed, na.rm = TRUE),
    p_wilcox_signed = wt_signed$p.value,
    median_abs_binary = median(df$abs_err_Binary, na.rm = TRUE),
    median_abs_ordinal = median(df$abs_err_Ordinal, na.rm = TRUE),
    median_delta_abs = median(d_abs, na.rm = TRUE),
    mean_delta_abs = mean(d_abs, na.rm = TRUE),
    sd_delta_abs = sd(d_abs, na.rm = TRUE),
    p_wilcox_abs = wt_abs$p.value
  )
}

by_group <- wide %>%
  group_by(sample_size, group) %>%
  group_modify(~ compare_one(.x)) %>%
  ungroup() %>%
  mutate(
    p_wilcox_signed_fdr = p.adjust(p_wilcox_signed, method = "BH"),
    p_wilcox_abs_fdr = p.adjust(p_wilcox_abs, method = "BH")
  ) %>%
  arrange(sample_size, group)

overall <- compare_one(wide) %>%
  mutate(sample_size = NA_real_, group = "Overall") %>%
  select(sample_size, group, everything())

write_tsv(by_group, file.path(output_root, "ordinal6_bias_binary_vs_ordinal_acc070_paired.tsv"))
write_tsv(overall, file.path(output_root, "ordinal6_bias_binary_vs_ordinal_acc070_overall.tsv"))

message("Wrote paired comparison tables for accuracy=0.70.")
