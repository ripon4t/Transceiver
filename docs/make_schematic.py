"""Draws the link block diagram used in the README."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OPT, ELE = "#2a6fbb", "#e07a1f"
fig, ax = plt.subplots(figsize=(11, 3.6)); ax.axis("off"); ax.set_xlim(0, 11); ax.set_ylim(0, 3.6)


def box(x, y, w, text, col):
    ax.add_patch(FancyBboxPatch((x, y), w, 0.8, boxstyle="round,pad=0.05", fc="white", ec=col, lw=2))
    ax.text(x + w / 2, y + 0.4, text, ha="center", va="center", fontsize=9.5)


def arrow(x0, y0, x1, y1, col):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=14, color=col, lw=2))


box(0.1, 2.3, 1.3, "PRBS\ngenerator", ELE); box(1.8, 2.3, 1.3, "NRZ pulse\ngenerator", ELE)
box(0.1, 0.5, 1.3, "CW laser", OPT); box(3.5, 1.4, 1.5, "Mach-Zehnder\nmodulator", OPT)
box(5.4, 1.4, 1.3, "Attenuator\n(channel)", OPT); box(7.1, 1.4, 1.3, "PIN\nphotodiode", OPT)
box(8.8, 1.4, 1.2, "Bessel\nfilter", ELE)
box(8.4, 2.7, 1.1, "Scope", ELE); box(9.7, 2.7, 1.2, "Eye diagram\nanalyzer", ELE)
arrow(1.4, 2.7, 1.8, 2.7, ELE); arrow(3.1, 2.7, 3.9, 2.2, ELE)
arrow(1.4, 0.9, 3.5, 1.7, OPT); arrow(5.0, 1.8, 5.4, 1.8, OPT); arrow(6.7, 1.8, 7.1, 1.8, OPT)
arrow(8.4, 1.8, 8.8, 1.8, OPT); arrow(9.4, 2.2, 8.95, 2.7, ELE); arrow(9.6, 2.2, 10.2, 2.7, ELE)
ax.plot([], [], color=OPT, lw=2, label="optical"); ax.plot([], [], color=ELE, lw=2, label="electrical")
ax.legend(loc="lower right", frameon=False)
ax.set_title("NRZ transceiver link simulated in INTERCONNECT")
fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "link_schematic.png"), dpi=200)
