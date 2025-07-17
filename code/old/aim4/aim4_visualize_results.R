# -----------------------------------------------------------------------------
# aim4_visualize_results.R
#
# This script loads the results from the Aim 4 closed-form analyses
# and creates heatmap plots to visualize the NNS under different
# multiplier scenarios for different subgroups.
# -----------------------------------------------------------------------------

# --- 1. SETUP ---
message("--- Loading Packages ---")
pacman::p_load(tidyverse, scales, viridis)

theme_set(theme_minimal() + theme(legend.position = "bottom"))

# Create results directories
dir.create("results/plots", showWarnings = FALSE)

# --- 2. LOAD DATA ---
message("--- Loading Data ---")

results_c <- read_tsv("results/tables/aim4_closed_form_v2_summary.tsv", show_col_types = FALSE) %>%
  mutate(subgroup = "C")
results_b <- read_tsv("results/tables/aim4_closed_form_subgroupB_summary.tsv", show_col_types = FALSE) %>%
  mutate(subgroup = "B")
results_e <- read_tsv("results/tables/aim4_closed_form_subgroupE_summary.tsv", show_col_types = FALSE) %>%
  mutate(subgroup = "E")

# Combine results
all_results <- bind_rows(results_c, results_b, results_e) %>%
  mutate(log10_nns = log10(nns_closed_form))

# --- 3. CREATE PLOT FUNCTION ---
message("--- Creating Plot Function ---")

plot_nns_heatmap <- function(data, target_subgroup, test_type_filter) {
  
  plot_data <- data %>%
    filter(subgroup == target_subgroup, test_type == test_type_filter)
  
  p <- ggplot(plot_data, aes(x = event_rate_multiplier, y = proportion_multiplier, fill = log10_nns)) +
    geom_tile() +
    scale_fill_viridis(name = "log10(NNS)", option = "C", direction = -1, na.value = "grey50") +
    scale_x_continuous(breaks = seq(1, 2, 0.2)) +
    scale_y_continuous(breaks = seq(1, 2, 0.2)) +
    labs(
      title = paste("Log10(NNS) for Subgroup", target_subgroup),
      subtitle = paste("Test Type:", test_type_filter),
      x = "Overall Event Rate Multiplier",
      y = "Proportion Multiplier for Target Subgroup"
    ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  return(p)
}

# --- 4. GENERATE AND SAVE PLOTS ---
message("--- Generating and Saving Plots ---")

# Plot for Subgroup C
plot_c_perfect <- plot_nns_heatmap(all_results, "C", "Perfect (100%)")
plot_c_moderate <- plot_nns_heatmap(all_results, "C", "Balanced/Moderate")

ggsave("results/plots/aim4_heatmap_C_perfect.pdf", plot_c_perfect, width = 8, height = 6)
ggsave("results/plots/aim4_heatmap_C_moderate.pdf", plot_c_moderate, width = 8, height = 6)

# Plot for Subgroup E
plot_e_perfect <- plot_nns_heatmap(all_results, "E", "Perfect (100%)")
plot_e_moderate <- plot_nns_heatmap(all_results, "E", "Balanced/Moderate")

ggsave("results/plots/aim4_heatmap_E_perfect.pdf", plot_e_perfect, width = 8, height = 6)
ggsave("results/plots/aim4_heatmap_E_moderate.pdf", plot_e_moderate, width = 8, height = 6)


message("--- Plots saved to results/plots/ ---") 

# --- 5. CREATE PROPORTION MULTIPLIER PLOT ---
message("--- Creating Proportion Multiplier Line Plots ---")

# Filter for event_rate_multiplier == 1
proportion_only_data <- all_results %>%
  filter(event_rate_multiplier == 1)

# Create the plot
proportion_plot <- ggplot(proportion_only_data,
                          aes(x = proportion_multiplier,
                              y = nns_closed_form,
                              color = test_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10(labels = scales::label_log(),
                breaks = 10^seq(2, 12, by = 2)) +
  scale_color_viridis_d(name = "Test Type") +
  facet_wrap(~ subgroup, scales = "free_y", labeller = labeller(subgroup = function(x) paste("Subgroup", x))) +
  labs(
    title = "Impact of Proportion Multiplier on NNS (Event Rate Multiplier = 1)",
    subtitle = "Comparing different test types across subgroups",
    x = "Proportion of Events Attributed to Subgroup Multiplier",
    y = "NNS (log scale)"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

# Save the plot
ggsave("results/plots/aim4_proportion_multiplier_vs_nns.pdf",
       plot = proportion_plot,
       width = 14,
       height = 8)

message("--- Visualization script finished ---") 