"""Shared matplotlib style for publication-ready plots with LaTeX fonts."""

import matplotlib.pyplot as plt

# Colorblind-safe palette (Wong 2011)
CB_BLACK   = "#000000"
CB_ORANGE  = "#E69F00"
CB_SKY     = "#56B4E9"
CB_GREEN   = "#009E73"
CB_YELLOW  = "#F0E442"
CB_BLUE    = "#0072B2"
CB_VERMIL  = "#D55E00"
CB_PINK    = "#CC79A7"

COLORS = [CB_BLUE, CB_ORANGE, CB_GREEN, CB_VERMIL, CB_SKY, CB_PINK, CB_YELLOW, CB_BLACK]


def fmt_setting(name: str, label_map: dict[str, str] | None = None) -> str:
    r"""Format a setting name for display: use label_map if provided, else \texttt."""
    if label_map and name in label_map:
        return label_map[name]
    escaped = name.replace("_", r"\_")
    return rf"\texttt{{{escaped}}}"


def setup_latex_style():
    """Configure matplotlib for publication-quality LaTeX-rendered plots."""
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.grid": True,
        "grid.alpha": 0.3,
    })
