"""Shared matplotlib style for the paper figures (Optica single column).

Text width of the Optica template and of Supplement 1 is 5.25 in.  The main
text is set in Times, so figures use a Times face with STIX math at 8.5 pt;
fonts are embedded as TrueType (Type 42) and rasters at 300 dpi.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CODE = HERE.parents[1] / "cascaded_chi3"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

WIDTH = 5.25  # inches

BLUE = "#1f4e9c"
RED = "#b2352b"
GREEN = "#2e7d32"
ORANGE = "#d9822b"
GREY = "0.45"

RC = {
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "TeX Gyre Termes",
                   "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "lines.linewidth": 1.1,
    "lines.markersize": 3.2,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "axes.grid": False,
}


def apply():
    """(Re)apply the paper style; call after importing project modules that
    set their own matplotlib style."""
    plt.rcParams.update(RC)


apply()


def panel_label(ax, text, x=0.0, y=1.02, **kw):
    """Panel letter just above the top-left corner of the axes."""
    ax.text(x, y, text, transform=ax.transAxes, va="bottom", ha="left", **kw)


def placeholder(ax, text):
    """Empty framed panel for content that depends on the chosen design."""
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linestyle((0, (3, 2)))
        s.set_color("0.6")
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center",
            va="center", color="0.45", fontsize=7.5)


def save(fig, name):
    out = HERE / name
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"wrote {out}.pdf/.png")
