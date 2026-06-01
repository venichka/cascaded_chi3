"""2D design map of cascaded phase-per-pump in (Q_1, Q_2) space.

Companion to the 1D `fig_phase_shift_per_pump`. The cascade has two
free Q-knobs; the native architecture has only Q_1. Plotting in the
(Q_1, Q_2) plane shows the whole design landscape at once:

  - colour: log10(phi_NL / |s_+|^2) for the cascade scheme. Higher = better.
    phi_cas ~ Q_1^2 * Q_2 * (host structural factor).
  - white dashed contours: constant operating bandwidth
    Delta f = min(gamma_1, gamma_2) / (2 pi). The contours are
    L-shaped: vertical where gamma_2 (Q_2) is the bottleneck,
    horizontal where gamma_1 (Q_1) is the bottleneck, meeting at the
    diagonal Q_1 = Q_2/2 where the two losses balance.
  - red horizontal: Q_2 = Q_2_crossover, the threshold above which
    the cascade beats the native architecture at every Q_1.
  - native scheme is the Q_2 -> 0 limit: phi_nat depends only on Q_1.
    On this map, native iso-phase contours are vertical lines; we
    overlay them as thin grey verticals so the reader can compare
    native and cascade at any (Q_1, Q_2) point directly.

We do two panels: LiNbO_3 (cascade wins easily) and AlGaAs
(native wins until Q_2 is heroic).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np
from matplotlib.colors import LogNorm

from params import ETA_GEOM, HOSTS, OMEGA1_REAL

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def enhancement(Q2, host, eta=ETA_GEOM):
    """R = |chi3_casc|/|chi3_native| at Delta omega = 2 gamma_2."""
    return (Q2 / (24.0 * host["n"] ** 2)) * eta * host["chi2"] ** 2 / host["chi3"]


def Q2_crossover(host, eta=ETA_GEOM):
    """Q_2 at which R = 1, at Delta omega = 2 gamma_2."""
    return 24.0 * host["n"] ** 2 * host["chi3"] / (eta * host["chi2"] ** 2)


def phi_cascade_grid(Q1_grid, Q2_grid, host):
    """phi/pump for cascade, in the same arbitrary units as the 1D fig
    (chi^(3)/n^4 prefactor times Q_1^2 times R(Q_2))."""
    return Q1_grid ** 2 * enhancement(Q2_grid, host) * host["chi3"] / host["n"] ** 4


def phi_native_1d(Q1, host):
    return Q1 ** 2 * host["chi3"] / host["n"] ** 4


def main():
    Q1_axis = np.geomspace(10, 1e5, 220)
    Q2_axis = np.geomspace(100, 1e5, 220)
    Q1_grid, Q2_grid = np.meshgrid(Q1_axis, Q2_axis)

    bw = np.minimum(OMEGA1_REAL / Q1_grid,
                    2 * OMEGA1_REAL / Q2_grid) / (2 * np.pi)

    fig, axs = plt.subplots(1, 2, figsize=(13.0, 5.6), sharey=True,
                            gridspec_kw={"wspace": 0.08})

    # Use the same colour scale across both panels so the materials are
    # directly comparable.
    all_phis = []
    for host_name in ("LiNbO3", "AlGaAs"):
        all_phis.append(phi_cascade_grid(Q1_grid, Q2_grid, HOSTS[host_name]))
    vmin = max(min(arr.min() for arr in all_phis), 1e-25)
    vmax = max(arr.max() for arr in all_phis)

    bw_levels_hz = [1e9, 1e10, 1e11, 1e12, 1e13]
    bw_labels = {1e9: r"1 GHz", 1e10: r"10 GHz", 1e11: r"100 GHz",
                 1e12: r"1 THz", 1e13: r"10 THz"}

    for ax, host_name in zip(axs, ("LiNbO3", "AlGaAs")):
        host = HOSTS[host_name]
        phi_cas = phi_cascade_grid(Q1_grid, Q2_grid, host)

        pcm = ax.pcolormesh(Q1_axis, Q2_axis, phi_cas,
                            norm=LogNorm(vmin=vmin, vmax=vmax),
                            cmap="viridis", shading="gouraud")

        # Iso-bandwidth contours (white dashed L-shapes)
        cs_bw = ax.contour(Q1_axis, Q2_axis, bw, levels=bw_levels_hz,
                           colors="white", linewidths=1.0, linestyles="--",
                           alpha=0.85)
        ax.clabel(cs_bw, fmt={lv: bw_labels[lv] for lv in bw_levels_hz},
                  fontsize=8.5, inline=True, inline_spacing=4)

        # Diagonal Q_1 = Q_2/2: boundary between FM- and SH-bottleneck
        ax.loglog(Q1_axis, 2 * Q1_axis, "-", color="0.8", lw=0.8,
                  alpha=0.7)

        # Cascade-vs-native horizontal threshold
        Q2x = Q2_crossover(host)
        if Q2_axis[0] <= Q2x <= Q2_axis[-1]:
            ax.axhline(Q2x, color="C3", ls="-", lw=2.0)
            ax.text(0.97 * Q1_axis[-1], Q2x * 1.18,
                    rf"$Q_2^{{\rm crossover}} = {Q2x:.0f}$  "
                    rf"(cascade $\geq$ native above)",
                    color="C3", fontsize=9, ha="right",
                    bbox=dict(facecolor="white", edgecolor="C3",
                              boxstyle="round,pad=0.25", alpha=0.9))

        # Mark the LiNbO3 design point Q_1 = 1000, Q_2 = 4000
        ax.plot([1000], [4000], "o", color="white", mec="k", mew=1.2, ms=10)
        ax.annotate("design pt.\n$(Q_1, Q_2) = (10^3, 4\\!\\times\\!10^3)$",
                    xy=(1000, 4000), xytext=(1.5e3, 1.4e4),
                    fontsize=8.5, color="k",
                    bbox=dict(facecolor="white", edgecolor="0.5",
                              boxstyle="round,pad=0.25"),
                    arrowprops=dict(arrowstyle="-", color="0.4", lw=0.7))

        R_at_design = enhancement(4000, host)
        ax.set_title(f"{host_name}  ($R(Q_2{{=}}4000) \\approx {R_at_design:.1f}$)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"FM-cavity quality factor  $Q_1$")
        if ax is axs[0]:
            ax.set_ylabel(r"SH-cavity quality factor  $Q_2$")
        ax.set_xlim(Q1_axis[0], Q1_axis[-1])
        ax.set_ylim(Q2_axis[0], Q2_axis[-1])

    cbar = fig.colorbar(pcm, ax=axs, location="right", pad=0.015,
                        fraction=0.04, aspect=30)
    cbar.set_label(r"cascade $\phi_{\rm NL}/|s_+|^2$  (arb., log scale)")

    fig.suptitle(
        "Design landscape: cascade $\\phi_{\\rm NL}/|s_+|^2$ on the "
        r"$(Q_1, Q_2)$ plane  -- "
        r"shared colour scale; iso-bandwidth contours (white dashed); "
        r"cascade-vs-native crossover (red)",
        fontsize=11, y=1.00,
    )

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "phase_shift_per_pump_2d"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
