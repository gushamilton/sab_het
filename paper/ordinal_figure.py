import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from pathlib import Path
import os

# ----- Data: legacy 5-level sensitivity distributions -----
# This helper is retained only for the 5-point sensitivity comparison. The
# main manuscript-facing ordinal figures are now built in R from the 6-point
# setup in code/30_build_paper_ordinal_shift_figure.R.
groups = ["A", "B", "C", "D", "E"]
baseline_dead = np.array([0.217, 0.075, 0.192, 0.182, 0.029], dtype=float)
or_shrunk = np.array([1.00, 4.34, 0.89, 1.19, 0.56], dtype=float)
survivor_split = np.array([0.05, 0.10, 0.20, 0.50], dtype=float)
survivor_split = survivor_split / survivor_split.sum()

# Order is worst -> best
categories = [
    "Death", "ICU/ventilated", "Still hospitalised", "Discharged with complications", "Discharged well"
]

def build_cum_control(p_dead):
    nondeath = np.cumsum((1 - p_dead) * survivor_split)
    return np.concatenate(([p_dead], p_dead + nondeath[:-1]))

def probs_from_cum(cum):
    return np.diff(np.concatenate(([0], cum, [1])))

def subgroup_po_probs(p_dead, or_value):
    cum_ctrl = build_cum_control(p_dead)
    theta = np.log(cum_ctrl / (1 - cum_ctrl))
    cum_trt = 1 / (1 + np.exp(-(theta + np.log(or_value))))
    return probs_from_cum(cum_ctrl), probs_from_cum(cum_trt)

def subgroup_nonpo_probs_threshold(p_dead, or_death):
    cum_ctrl = build_cum_control(p_dead)
    theta = np.log(cum_ctrl / (1 - cum_ctrl))
    betas = np.concatenate(([np.log(or_death)], np.zeros(len(theta) - 1, dtype=float)))
    cum_trt = 1 / (1 + np.exp(-(theta + betas)))
    return probs_from_cum(cum_ctrl), probs_from_cum(cum_trt)

def subgroup_nonpo_probs_rescaled_survivors(p_dead, or_death):
    # Alternative "death-only" DGM:
    # apply OR at death, then keep the relative composition of survivor levels fixed.
    odds0 = p_dead / (1 - p_dead)
    odds1 = odds0 * or_death
    p_dead_t = odds1 / (1 + odds1)

    control = probs_from_cum(build_cum_control(p_dead))
    survivor_ctrl = control[1:]
    survivor_weights = survivor_ctrl / survivor_ctrl.sum()

    treatment = np.empty_like(control)
    treatment[0] = p_dead_t
    treatment[1:] = (1 - p_dead_t) * survivor_weights
    return control, treatment

target_group = os.getenv("TARGET_GROUP", "D")
if target_group not in groups:
    raise ValueError(f"TARGET_GROUP must be one of {groups}, got {target_group}")

target_idx = groups.index(target_group)
target_p0 = baseline_dead[target_idx]
default_or = or_shrunk[target_idx]
target_or = float(os.getenv("TARGET_OR", f"{default_or}"))
output_tag = os.getenv("OUTPUT_TAG", "").strip()

nonpo_style = os.getenv("NONPO_STYLE", "rescaled_survivors").strip().lower()
if nonpo_style not in {"rescaled_survivors", "threshold_only"}:
    raise ValueError("NONPO_STYLE must be 'rescaled_survivors' or 'threshold_only'")

control_po, treat_po = subgroup_po_probs(target_p0, target_or)
if nonpo_style == "rescaled_survivors":
    control_do, treat_do = subgroup_nonpo_probs_rescaled_survivors(target_p0, target_or)
    nonpo_subtitle = "Death probability changes; non-death categories rescaled proportionally."
else:
    control_do, treat_do = subgroup_nonpo_probs_threshold(target_p0, target_or)
    nonpo_subtitle = "Death-threshold-only cumulative effect (current threshold-parameterized model)."

# Muted blue palette, light -> dark
colors = ["#e8eef2", "#cadbe6", "#a7c4d8", "#7ea8c4", "#4f88a8", "#1f5a7a"]

def draw_panel(ax, control, treatment, title, subtitle):
    y_positions = [1, 0]
    labels = ["Control", "Treatment"]

    for y, vals, lab in zip(y_positions, [control, treatment], labels):
        left = 0
        for v, c in zip(vals, colors):
            ax.barh(
                y, v, left=left, height=0.42,
                color=c, edgecolor="white", linewidth=1.2
            )
            if v >= 0.075:
                ax.text(
                    left + v / 2, y, f"{v*100:.0f}%",
                    ha="center", va="center",
                    color="white" if c in colors[-2:] else "#1f2933",
                    fontsize=10, fontweight="bold"
                )
            left += v

        ax.text(
            -0.01, y, lab, ha="right", va="center",
            fontsize=12, fontweight="bold", color="#2b2b2b"
        )

    # Connect corresponding categories
    c_cum = np.cumsum(control)
    t_cum = np.cumsum(treatment)
    c_left = np.concatenate([[0], c_cum[:-1]])
    t_left = np.concatenate([[0], t_cum[:-1]])
    c_centers = c_left + control / 2
    t_centers = t_left + treatment / 2

    for i in range(len(categories)):
        ax.plot(
            [c_centers[i], t_centers[i]], [1, 0],
            ls=(0, (3, 3)), lw=1, color="#9aa5b1", zorder=0
        )

    ax.set_xlim(-0.12, 1)
    ax.set_ylim(-0.65, 1.65)
    ax.set_yticks([])
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_xticklabels([f"{int(x*100)}" for x in np.linspace(0, 1, 6)], fontsize=10)
    ax.set_xlabel("Percentage of participants", fontsize=11, fontweight="bold", color="#2b2b2b")
    ax.set_title(title, fontsize=15, fontweight="bold", loc="left", pad=18, color="#1f2933")

    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#7b8794")
    ax.tick_params(axis="x", colors="#2b2b2b")
    ax.grid(False)

# ----- Figure layout -----
fig = plt.figure(figsize=(15.5, 4.8), facecolor="none")
gs = fig.add_gridspec(1, 2, wspace=0.18)

ax1 = fig.add_subplot(gs[0, 0], facecolor="none")
ax2 = fig.add_subplot(gs[0, 1], facecolor="none")

draw_panel(
    ax1,
    control_po,
    treat_po,
    "A. Proportional odds effect shift",
    ""
)

draw_panel(
    ax2,
    control_do,
    treat_do,
    "B. Non-proportional odds effect",
    ""
)

# Compact shared legend
legend_handles = [Patch(facecolor=col, edgecolor="none", label=cat) for cat, col in zip(categories, colors)]
fig.legend(
    handles=legend_handles,
    loc="lower center",
    ncol=len(categories),
    frameon=False,
    bbox_to_anchor=(0.5, -0.02),
    fontsize=10.5
)

fig.subplots_adjust(bottom=0.22)

plt.show()

# Save into paper/ regardless of execution cwd
out_dir = Path(__file__).resolve().parent
suffix = f"_{output_tag}" if output_tag else ""
fig.savefig(
    out_dir / f"ordinal_outcome_two_panel{suffix}.png",
    dpi=300,
    bbox_inches="tight",
    facecolor="none",
    transparent=True,
)
fig.savefig(
    out_dir / f"ordinal_outcome_two_panel{suffix}.pdf",
    bbox_inches="tight",
    facecolor="none",
    transparent=True,
)
