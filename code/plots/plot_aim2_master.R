# Patchwork: Aim 2 power, bias, and replicate betas

pacman::p_load(tidyverse, scales, patchwork, ggdist)

theme_set(
  ggdist::theme_ggdist() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "white", color = NA)
    )
)

scenario_set <- Sys.getenv("SCENARIO_SET", unset = "main")
if (!scenario_set %in% c("main", "supp")) {
  stop("SCENARIO_SET must be 'main' or 'supp'")
}

results_root <- file.path("results", scenario_set)
summary_file <- file.path(results_root, "tables", "aim2_accuracy_summary.tsv")
raw_file <- file.path(results_root, "tables", "aim2_raw_results.tsv.gz")
output_file <- file.path(results_root, "plots", "aim2_master.pdf")

if (!file.exists(summary_file)) {
  stop("Missing input: ", summary_file)
}
if (!file.exists(raw_file)) {
  stop("Missing input: ", raw_file)
}

aim2 <- readr::read_tsv(summary_file, show_col_types = FALSE) %>%
  filter(group != "Overall")

acc_breaks <- sort(unique(aim2$accuracy))
acc_limits <- range(acc_breaks)

power_plot <- aim2 %>%
  ggplot(aes(x = accuracy, y = power, color = group)) +
  geom_line() +
  scale_x_continuous(breaks = acc_breaks, limits = acc_limits) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Accuracy", y = "Power", color = "Group")

bias_plot <- aim2 %>%
  ggplot(aes(x = accuracy, y = bias, color = group)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = acc_breaks, limits = acc_limits) +
  labs(x = "Accuracy", y = "Bias (log OR)", color = "Group")

replicates <- readr::read_tsv(raw_file, show_col_types = FALSE) %>%
  mutate(accuracy_round = round(accuracy, 2)) %>%
  filter(group != "Overall", accuracy_round %in% c(1.0, 0.9, 0.8, 0.7)) %>%
  mutate(
    group = factor(group, levels = c("A", "B", "C", "D", "E")),
    accuracy = factor(
      accuracy_round,
      levels = c(1.0, 0.9, 0.8, 0.7),
      labels = c("Accuracy = 1.0", "Accuracy = 0.9", "Accuracy = 0.8", "Accuracy = 0.7")
    )
  )

true_beta_df <- replicates %>%
  distinct(group, accuracy, true_beta)

replicate_plot <- ggplot(replicates, aes(x = group, y = beta, color = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, coef = 1.5) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.15, size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_segment(
    data = true_beta_df,
    aes(
      x = as.numeric(group) - 0.3,
      xend = as.numeric(group) + 0.3,
      y = true_beta,
      yend = true_beta
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.7
  ) +
  facet_wrap(~ accuracy, ncol = 2, strip.position = "top") +
  labs(x = "Subphenotype", y = "Beta (log OR)", color = "Group") +
  theme(strip.text = element_text(size = 11, face = "bold"))

final_plot <- ((power_plot | bias_plot) / replicate_plot) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", heights = c(1, 1.3)) &
  theme(legend.position = "bottom")

ggsave(output_file, final_plot, width = 9, height = 12)
message("Saved: ", output_file)
