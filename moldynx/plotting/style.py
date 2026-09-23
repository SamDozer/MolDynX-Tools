"""Global matplotlib style + palette (Nature-like: white bg, large fonts, no grid)."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PALETTE = {
    "primary": "#1f6f8b", "secondary": "#e07a5f", "accent": "#3d405b",
    "green": "#5f9e6e", "purple": "#8367c7", "muted": "#9aa0a6",
    "helix": "#c0392b", "sheet": "#2471a3", "coil": "#7f8c8d",
}
SEQ_CMAP = "viridis"
DIV_CMAP = "RdBu_r"
FEL_CMAP = "nipy_spectral"


def set_style() -> None:
    """Apply the global publication style. Call once per script/session."""
    plt.rcParams.update({
        "figure.figsize": (7.2, 4.8), "figure.dpi": 110, "savefig.dpi": 300,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "font.size": 13, "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica", "sans-serif"],
        "axes.titlesize": 15, "axes.labelsize": 14, "axes.titleweight": "bold",
        "xtick.labelsize": 12, "ytick.labelsize": 12,
        "legend.fontsize": 11, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 1.1, "axes.grid": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 5, "ytick.major.size": 5,
        "lines.linewidth": 1.8, "lines.solid_capstyle": "round",
        "legend.handlelength": 1.6,
    })


def new_axes(figsize: tuple[float, float] = (7.2, 4.8)):
    return plt.subplots(figsize=figsize)


def add_reference_line(ax, y, label=None, color=None):
    ax.axhline(y, ls="--", lw=1.2, color=color or PALETTE["accent"],
               alpha=0.8, label=label, zorder=1)
