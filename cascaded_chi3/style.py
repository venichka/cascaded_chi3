"""Shared matplotlib styling for the cascaded-chi3 figures.

Call apply_style() once at the top of each fig_*.py (after the
matplotlib.use('Agg') backend selection but before importing pyplot from
elsewhere). The COLORS dict is the canonical palette: TCMT numerics in
blue, analytic curves in black, approximations in grey, failure modes
in red, special points in orange.
"""

from __future__ import annotations

import matplotlib as mpl


COLORS = {
    "tcmt":      "#1f77b4",   # numerical TCMT (blue, C0)
    "analytic":  "#000000",   # closed-form analytic (black)
    "adiabatic": "#666666",   # approximation / adiabatic (grey)
    "wrong":     "#d62728",   # failure mode / "wrong" (red, C3)
    "highlight": "#ff7f0e",   # special points / annotations (orange, C1)
    "linear":    "#888888",   # linear-cavity reference (light grey)
    "second":    "#2ca02c",   # second curve / partner (green, C2)
}


MARKERS = {
    "tcmt": dict(marker="o", ms=4.0, mfc="white", mew=1.2, ls="none"),
    "tcmt_filled": dict(marker="o", ms=4.5, mew=0.0, ls="none"),
    "wrong_marker": dict(marker="s", ms=4.0, mfc="white", mew=1.2, ls="none"),
}


def apply_style():
    mpl.rcParams.update({
        # text
        "font.size": 11.5,
        "axes.titlesize": 12.5,
        "axes.labelsize": 11.5,
        "legend.fontsize": 10,
        "legend.title_fontsize": 10,
        "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5,
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "axes.titleweight": "normal",
        "axes.titlepad": 8,
        # lines / markers
        "lines.linewidth": 1.8,
        "lines.markersize": 5,
        "lines.markeredgewidth": 1.0,
        # grid & spines
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.6,
        "grid.linestyle": "-",
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.linewidth": 0.9,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        # legend
        "legend.frameon": True,
        "legend.framealpha": 0.95,
        "legend.fancybox": False,
        "legend.edgecolor": "0.7",
        "legend.borderpad": 0.4,
        "legend.handlelength": 1.8,
        # output
        "figure.dpi": 110,
        "savefig.dpi": 160,
        "savefig.bbox": "tight",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })
