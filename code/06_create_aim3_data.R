# Create Aim 3 combined data for the paper
# This script combines parallel simulation results and closed-form calculations
# to create comprehensive Aim 3 analysis data

library(tidyverse)
library(gt)
library(scales)
library(patchwork)

# Load the data
parallel_results <- read_tsv("parallel/results/tables/aim3_sens_spec_summary.tsv")
closed_form_results <- read_tsv("results/tables/aim3_closed_form_summary.tsv")

# Clean and prepare the data
parallel_clean <- parallel_results %>%
  filter(!is.na(nns_needed)) %>%  # Remove cases where power goal wasn't met
  mutate(
    data_source = "Simulation",
    test_type_clean = case_when(
      test_type == "Perfect (100%)" ~ "Perfect",
      test_type == "Near-Perfect" ~ "Near-Perfect", 
      test_type == "High Sens/High Spec" ~ "High Sens/High Spec",
      test_type == "High Sens/Low Spec" ~ "High Sens/Low Spec",
      test_type == "Low Sens/High Spec" ~ "Low Sens/High Spec",
      test_type == "Balanced/Moderate" ~ "Balanced/Moderate",
      TRUE ~ test_type
    )
  )

closed_form_clean <- closed_form_results %>%
  mutate(
    data_source = "Closed Form",
    test_type_clean = case_when(
      test_type == "Perfect (100%)" ~ "Perfect",
      test_type == "Near-Perfect" ~ "Near-Perfect",
      test_type == "High Sens/High Spec" ~ "High Sens/High Spec", 
      test_type == "High Sens/Low Spec" ~ "High Sens/Low Spec",
      test_type == "Low Sens/High Spec" ~ "Low Sens/High Spec",
      test_type == "Balanced/Moderate" ~ "Balanced/Moderate",
      TRUE ~ test_type
    )
  )

# Combine the datasets for comparison
combined_results <- bind_rows(
  parallel_clean %>% 
    select(scenario_name, test_type_clean, target_group, 
           nns_needed, nnr_corresponding, bias, mse, rmse, 
           wrong_direction, data_source),
  closed_form_clean %>%
    select(scenario_name, test_type_clean, target_group,
           nns_closed_form, nnr_closed_form, bias_closed_form,
           true_beta, predicted_beta, data_source) %>%
    rename(nns_needed = nns_closed_form, 
           nnr_corresponding = nnr_closed_form,
           bias = bias_closed_form) %>%
    mutate(mse = NA, rmse = NA, wrong_direction = NA)
)

# Create summary table for the paper
summary_table <- combined_results %>%
  filter(target_group %in% c("B", "C", "D", "E")) %>%  # Include group D
  group_by(scenario_name, test_type_clean, target_group, data_source) %>%
  summarise(
    NNS = ifelse(all(is.na(nns_needed)), ">100,000", 
                 scales::comma(min(nns_needed, na.rm = TRUE))),
    NNR = ifelse(all(is.na(nnr_corresponding)), "-", 
                 scales::comma(round(mean(nnr_corresponding, na.rm = TRUE), 1))),
    Bias = ifelse(all(is.na(bias)), "-", 
                  sprintf("%.3f", mean(bias, na.rm = TRUE))),
    MSE = ifelse(all(is.na(mse)), "-", 
                 sprintf("%.3f", mean(mse, na.rm = TRUE))),
    Wrong_Dir = ifelse(all(is.na(wrong_direction)), "-", 
                       scales::percent(mean(wrong_direction, na.rm = TRUE), accuracy = 0.1)),
    # Add beta columns
    True_Beta = ifelse(all(is.na(true_beta)), "-", 
                       sprintf("%.3f", mean(true_beta, na.rm = TRUE))),
    Predicted_Beta = ifelse(all(is.na(predicted_beta)), "-", 
                            sprintf("%.3f", mean(predicted_beta, na.rm = TRUE))),
    Empirical_Beta = ifelse(all(is.na(bias)), "-", 
                            sprintf("%.3f", mean(true_beta, na.rm = TRUE) + mean(bias, na.rm = TRUE))),
    .groups = "drop"
  ) %>%
  arrange(scenario_name, target_group, test_type_clean, data_source)

# Create GT table for the paper
gt_table <- summary_table %>%
  gt(groupname_col = c("scenario_name", "target_group")) %>%
  cols_label(
    test_type_clean = "Test Type",
    data_source = "Method",
    NNS = "NNS",
    NNR = "NNR", 
    Bias = "Bias",
    MSE = "MSE",
    Wrong_Dir = "Wrong Dir %",
    True_Beta = "True β",
    Predicted_Beta = "Predicted β",
    Empirical_Beta = "Empirical β"
  ) %>%
  tab_header(title = "Aim 3: NNS, NNR, and Performance Metrics by Method") %>%
  tab_options(table.width = "100%")

# Create comparison plots
# NNS comparison
nns_comparison <- combined_results %>%
  filter(target_group %in% c("B", "C", "D", "E"), !is.na(nns_needed)) %>%
  ggplot(aes(x = test_type_clean, y = nns_needed, color = data_source, group = data_source)) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(scenario_name ~ target_group, scales = "free_y") +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "NNS Comparison: Simulation vs Closed Form",
    x = "Test Type",
    y = "Number Needed to Screen (log scale)",
    color = "Method"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

# NNR comparison  
nnr_comparison <- combined_results %>%
  filter(target_group %in% c("B", "C", "D", "E"), !is.na(nnr_corresponding)) %>%
  ggplot(aes(x = test_type_clean, y = nnr_corresponding, color = data_source, group = data_source)) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(scenario_name ~ target_group, scales = "free_y") +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "NNR Comparison: Simulation vs Closed Form", 
    x = "Test Type",
    y = "Number Needed to Randomize (log scale)",
    color = "Method"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

# Bias comparison (simulation only since closed form has theoretical bias)
bias_comparison <- parallel_clean %>%
  filter(target_group %in% c("B", "C", "D", "E"), !is.na(bias)) %>%
  ggplot(aes(x = test_type_clean, y = bias, color = target_group, group = target_group)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ scenario_name) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Bias by Test Type (Simulation Results)",
    x = "Test Type", 
    y = "Bias (log OR)",
    color = "Target Group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

# Save the results for use in QMD
saveRDS(summary_table, "results/objects/aim3_combined_summary.rds")
saveRDS(gt_table, "results/objects/aim3_combined_gt_table.rds")
saveRDS(nns_comparison, "results/objects/aim3_nns_comparison_plot.rds") 
saveRDS(nnr_comparison, "results/objects/aim3_nnr_comparison_plot.rds")
saveRDS(bias_comparison, "results/objects/aim3_bias_comparison_plot.rds")

# Also save the combined data for potential further analysis
saveRDS(combined_results, "results/objects/aim3_combined_data.rds")

# Save summary table as TSV for easy inspection
write_tsv(summary_table, "results/tables/aim3_combined_summary.tsv")

cat("Aim 3 results created and saved to results/objects/ and results/tables/\n")
cat("Files created:\n")
cat("- results/objects/aim3_combined_summary.rds\n")
cat("- results/objects/aim3_combined_gt_table.rds\n")
cat("- results/objects/aim3_nns_comparison_plot.rds\n")
cat("- results/objects/aim3_nnr_comparison_plot.rds\n")
cat("- results/objects/aim3_bias_comparison_plot.rds\n")
cat("- results/objects/aim3_combined_data.rds\n")
cat("- results/tables/aim3_combined_summary.tsv\n") 