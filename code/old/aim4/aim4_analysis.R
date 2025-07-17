# -----------------------------------------------------------------------------
# Aim 4 Analysis Script (Temporary)
#
# This script loads the Aim 4 results and creates plots showing the relationship
# between multipliers and NNS on a log scale.
# -----------------------------------------------------------------------------

# --- 1. SETUP ----------------------------------------------------------------

pacman::p_load(tidyverse, scales, patchwork)

theme_set(theme_minimal() + 
          theme(legend.position = "bottom",
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5)))

# --- 2. LOAD DATA -----------------------------------------------------------

# Load Aim 4 results
results_aim4 <- read_tsv("results/tables/aim4_comparison_summary.tsv")

# Check data structure
cat("Data loaded with", nrow(results_aim4), "rows\n")
cat("Unique test types:", unique(results_aim4$test_type), "\n")
cat("Multiplier range:", range(results_aim4$multiplier), "\n")
cat("NNS range:", range(results_aim4$nns, na.rm = TRUE), "\n")

# --- 3. CREATE PLOTS --------------------------------------------------------

# Create color palette for test types (6 colors)
test_colors <- c("Perfect (100%)" = "#1f77b4", 
                 "Near-Perfect" = "#ff7f0e", 
                 "High Sens/High Spec" = "#2ca02c",
                 "Balanced/Moderate" = "#d62728",
                 "Low Sens/High Spec" = "#9467bd",
                 "High Sens/Low Spec" = "#8c564b")

# Plot 1: Event Rate Multipliers
p1 <- results_aim4 %>%
  filter(scenario == "Event Rate Multiplier") %>%
  ggplot(aes(x = multiplier, y = nns, color = subgroup, shape = test_type)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_log10(labels = comma) +
  scale_color_manual(values = test_colors) +
  scale_shape_manual(values = c(16, 17, 15, 18, 19, 20)) +
  labs(title = "Event Rate Multipliers vs NNS",
       subtitle = "Effect of scaling all baseline event rates",
       x = "Event Rate Multiplier",
       y = "Number Needed to Screen (log scale)",
       color = "Test Type",
       shape = "Test Type") +
  theme(legend.position = "bottom") + facet_wrap(~test_type)
p1
# Plot 2: Proportion Multipliers
p2 <- results_aim4 %>%
  filter(scenario == "Proportion Multiplier") %>%
  ggplot(aes(x = multiplier, y = nns, color = test_type, shape = test_type)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(alpha = 0.5) +
  scale_y_log10(labels = comma) +
  scale_color_manual(values = test_colors) +
  scale_shape_manual(values = c(16, 17, 15, 18, 19, 20)) +
  labs(title = "Proportion Multipliers vs NNS",
       subtitle = "Effect of increasing proportion of events attributed to Group C",
       x = "Proportion Multiplier",
       y = "Number Needed to Screen (log scale)",
       color = "Test Type",
       shape = "Test Type") +
  theme(legend.position = "bottom")
p2
# Combined plot
combined_plot <- p1 + p2 + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# --- 4. SAVE PLOTS ----------------------------------------------------------

# Ensure output directory exists
dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

# Save individual plots
ggsave("results/plots/aim4_event_rate_multipliers.pdf", p1, 
       width = 8, height = 6, dpi = 300)
ggsave("results/plots/aim4_proportion_multipliers.pdf", p2, 
       width = 8, height = 6, dpi = 300)
ggsave("results/plots/aim4_combined_comparison.pdf", combined_plot, 
       width = 12, height = 6, dpi = 300)

cat("Plots saved to results/plots/\n")

# --- 5. SUMMARY STATISTICS --------------------------------------------------

# Summary by scenario and test type
summary_stats <- results_aim4 %>%
  group_by(scenario, test_type) %>%
  summarise(
    n_scenarios = n(),
    nns_min = min(nns, na.rm = TRUE),
    nns_max = max(nns, na.rm = TRUE),
    nns_median = median(nns, na.rm = TRUE),
    nns_mean = mean(nns, na.rm = TRUE),
    bias_range = max(bias_closed_form, na.rm = TRUE) - min(bias_closed_form, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(scenario, test_type)

cat("\n=== SUMMARY STATISTICS ===\n")
print(summary_stats)

# Compare baseline (multiplier = 1) vs maximum multiplier
baseline_comparison <- results_aim4 %>%
  filter(multiplier == 1.0) %>%
  select(scenario, test_type, nns_baseline = nns) %>%
  left_join(
    results_aim4 %>%
      group_by(scenario, test_type) %>%
      filter(multiplier == max(multiplier)) %>%
      select(scenario, test_type, nns_max_mult = nns, max_multiplier = multiplier),
    by = c("scenario", "test_type")
  ) %>%
  mutate(
    nns_ratio = nns_max_mult / nns_baseline,
    log_nns_ratio = log(nns_ratio)
  )

cat("\n=== BASELINE COMPARISON (multiplier = 1.0 vs max) ===\n")
print(baseline_comparison)

# --- 6. DISPLAY PLOTS -------------------------------------------------------

# Show plots in R session
print(p1)
print(p2)
print(combined_plot)

cat("\nAnalysis complete!\n") 