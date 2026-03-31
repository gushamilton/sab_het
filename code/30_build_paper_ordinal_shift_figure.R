## =====================================================================
##  30_build_paper_ordinal_shift_figure.R -- Manuscript Figure 4 artwork
## =====================================================================

library(tidyverse)
library(patchwork)
library(ggdist)

source("code/01_parameters.R")
source("code/functions/01_simulation_helpers.R")

params <- get_parameters()
subphenotypes <- params$subphenotype_table

or_vector <- setNames(subphenotypes$or_arrest_shrunk, subphenotypes$subphenotype)
baseline <- setNames(subphenotypes$baseline_mortality, subphenotypes$subphenotype)
ordinal_points <- as.integer(Sys.getenv("ORDINAL_POINTS", unset = "6"))
if (ordinal_points == 6) {
  ordinal_labels <- c(
    "Dead",
    "ICU/ventilated",
    "Still hospitalised",
    "Discharged to rehab",
    "Discharged with complications",
    "Discharged well"
  )
  survivor_split <- c(0.08, 0.12, 0.20, 0.25, 0.35)
} else if (ordinal_points == 5) {
  ordinal_labels <- params$ordinal_baseline$label
  survivor_split <- derive_survivor_split(params$ordinal_baseline$prevalence)
} else {
  stop("ORDINAL_POINTS must be 5 or 6.", call. = FALSE)
}

target_group <- Sys.getenv("TARGET_GROUP", unset = "B")
if (!(target_group %in% names(or_vector))) {
  stop("TARGET_GROUP must be one of: ", paste(names(or_vector), collapse = ", "), call. = FALSE)
}

rescaled_nonpo_probs <- function(p_dead, survivor_split, or_death) {
  control <- build_ordinal_probs_from_death(p_dead, survivor_split)

  odds0 <- p_dead / (1 - p_dead)
  odds1 <- odds0 * or_death
  p_dead_t <- odds1 / (1 + odds1)

  survivor_weights <- control[-1] / sum(control[-1])
  treatment <- c(p_dead_t, (1 - p_dead_t) * survivor_weights)

  list(control = control, treatment = treatment)
}

po_probs <- get_po_ordinal_probs(
  p_dead = baseline[[target_group]],
  survivor_split = survivor_split,
  or_value = or_vector[[target_group]]
)
nonpo_probs <- rescaled_nonpo_probs(
  p_dead = baseline[[target_group]],
  survivor_split = survivor_split,
  or_death = or_vector[[target_group]]
)

plot_df <- bind_rows(
  tibble(
    panel = "A. Proportional-odds model",
    arm = c("Placebo", "Rifampicin"),
    probs = list(po_probs$control, po_probs$treatment)
  ),
  tibble(
    panel = "B. Death-only non-proportional model",
    arm = c("Placebo", "Rifampicin"),
    probs = list(nonpo_probs$control, nonpo_probs$treatment)
  )
) %>%
  unnest_longer(probs, indices_to = "level", values_to = "probability") %>%
  mutate(
    label = factor(ordinal_labels[level], levels = ordinal_labels),
    arm = factor(arm, levels = c("Rifampicin", "Placebo")),
    panel = factor(panel, levels = c("A. Proportional-odds model", "B. Death-only non-proportional model"))
  )

fill_values <- c(
  "Dead" = "#16324f",
  "ICU/ventilated" = "#3b6c8e",
  "Still hospitalised" = "#6c9fb8",
  "Discharged to rehab" = "#88b4c8",
  "Discharged with complications" = "#a8c6d7",
  "Discharged well" = "#dce9ef"
)

plot_df <- plot_df %>%
  group_by(panel, arm) %>%
  arrange(level, .by_group = TRUE) %>%
  mutate(
    xmin = lag(cumsum(probability), default = 0),
    xmax = cumsum(probability),
    xmid = xmin + probability / 2,
    percent_label = if_else(probability >= 0.075, paste0(round(probability * 100), "%"), "")
  ) %>%
  ungroup()

p <- ggplot(plot_df) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = as.numeric(arm) - 0.28, ymax = as.numeric(arm) + 0.28, fill = label),
    color = "white",
    linewidth = 0.6
  ) +
  geom_text(
    aes(x = xmid, y = as.numeric(arm), label = percent_label, color = label),
    size = 3.2,
    lineheight = 0.9
  ) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = fill_values, name = NULL) +
  scale_color_manual(
    values = c(
      "Dead" = "white",
      "ICU/ventilated" = "white",
      "Still hospitalised" = "white",
      "Discharged to rehab" = "#1f2933",
      "Discharged with complications" = "#1f2933",
      "Discharged well" = "#1f2933"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = c(1, 2),
    labels = c("Rifampicin", "Placebo"),
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  labs(
    x = "Participants",
    y = NULL
  ) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  ggdist::theme_ggdist(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.4, "lines"),
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "plain"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 8)),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave("paper/ordinal_outcome_two_panel.pdf", p, width = 11, height = 4.8)
ggsave("paper/ordinal_outcome_two_panel.png", p, width = 11, height = 4.8, dpi = 300)
