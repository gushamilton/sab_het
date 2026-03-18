## =====================================================================
##  29_combine_ordinal_nonPO_chunks.R
##  Combine chunked outputs from 15_run_ordinal_nonPO_comparison.R
## =====================================================================

library(tidyverse)

source("code/01_parameters.R")
source("code/functions/02_metrics.R")

params <- get_parameters()
alpha <- params$alpha_primary

chunk_tags <- Sys.getenv("CHUNK_TAGS", unset = "")
final_tag <- Sys.getenv("FINAL_RUN_TAG", unset = "")

if (!nzchar(chunk_tags)) {
  stop("Set CHUNK_TAGS as comma-separated run tags.")
}
if (!nzchar(final_tag)) {
  stop("Set FINAL_RUN_TAG for merged outputs.")
}

tags <- strsplit(chunk_tags, ",", fixed = TRUE)[[1]] %>%
  trimws() %>%
  discard(~ !nzchar(.x))

if (length(tags) == 0) {
  stop("No valid chunk tags after parsing CHUNK_TAGS.")
}

output_root <- file.path("results", "supp", "ordinal_nonPO_comparison")
dir.create(file.path(output_root, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "plots"), showWarnings = FALSE, recursive = TRUE)

read_chunk_raw <- function(tag) {
  path <- file.path(output_root, "data", paste0("ordinal_nonPO_comparison_raw_", tag, ".tsv.gz"))
  if (!file.exists(path)) stop("Missing chunk raw file: ", path)
  read_tsv(path, show_col_types = FALSE) %>% mutate(chunk_tag = tag)
}

raw_all <- map_dfr(tags, read_chunk_raw) %>%
  arrange(sample_size, accuracy, dgm, group, model_type, rep_id)

write_tsv(
  raw_all %>% select(-chunk_tag),
  file.path(output_root, "data", paste0("ordinal_nonPO_comparison_raw_", final_tag, ".tsv.gz"))
)

metric_results <- raw_all %>%
  group_by(sample_size, accuracy, dgm, group, model_type) %>%
  group_modify(~ {
    truth <- unique(.x$log_or_true)
    is_null <- abs(truth) < 1e-6
    summarise_metrics(.x, log_or_true = truth, alpha = alpha, is_null = is_null)
  }) %>%
  ungroup()

write_tsv(
  metric_results,
  file.path(output_root, "tables", paste0("ordinal_nonPO_comparison_metrics_", final_tag, ".tsv"))
)

bias_results <- raw_all %>%
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
  file.path(output_root, "tables", paste0("ordinal_nonPO_comparison_bias_", final_tag, ".tsv"))
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
  file.path(output_root, "plots", paste0("ordinal_nonPO_comparison_power_", final_tag, ".pdf")),
  p_power, width = 14, height = 9
)

sample_size_focus <- max(raw_all$sample_size, na.rm = TRUE)
dist_data <- raw_all %>%
  filter(sample_size == sample_size_focus, group != "A") %>%
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
    title = paste0("Proportional bias distributions at N = ", sample_size_focus),
    x = NULL,
    y = "Proportional bias",
    fill = "Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_root, "plots", paste0("ordinal_nonPO_comparison_dist_", final_tag, ".pdf")),
  p_dist, width = 14, height = 9
)

or_vector <- setNames(params$subphenotype_table$or_arrest_shrunk, params$subphenotype_table$subphenotype)
estimand_tbl <- raw_all %>%
  distinct(group, dgm, model_type, log_or_true) %>%
  mutate(estimand = paste(model_type, dgm, sep = " / "))

# Fill any missing entries defensively.
template <- expand_grid(
  group = names(or_vector),
  dgm = c("PO", "nonPO"),
  model_type = c("Binary (death)", "Ordinal (polr)")
)

estimand_tbl <- template %>%
  left_join(estimand_tbl, by = c("group", "dgm", "model_type")) %>%
  mutate(
    log_or_true = case_when(
      !is.na(log_or_true) ~ log_or_true,
      model_type == "Binary (death)" ~ log(or_vector[group]),
      TRUE ~ NA_real_
    ),
    estimand = paste(model_type, dgm, sep = " / ")
  )

p_estimand <- estimand_tbl %>%
  ggplot(aes(x = group, y = log_or_true, fill = interaction(model_type, dgm, sep = " / "))) +
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
  file.path(output_root, "plots", paste0("ordinal_nonPO_comparison_estimands_", final_tag, ".pdf")),
  p_estimand, width = 10, height = 5
)

message("Combined nonPO chunks into run tag: ", final_tag)
