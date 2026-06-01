"""Material survey: cascaded enhancement ratio R vs Q_2, by host.

Following the table in Section 11 of the main note at the conservative
operating point Delta omega = 2 gamma_2 (= 4 omega_1 / Q_2, since
gamma_2 = omega_2/Q_2 = 2 omega_1/Q_2). For materials with chi^(2) != 0
(non-centrosymmetric, cascade available) the simple-form ratio is

    R(Q_2) = (Q_2 / (24 n_2^2)) * eta_geom * |chi^(2)|^2 / chi^(3)_native

grows linearly with Q_2 and crosses unity at

    Q_2^crossover = (24 n_2^2 chi^(3)_native) / (|chi^(2)|^2 eta_geom).

The factor 1/24 comes from omega_1/Delta omega = Q_2/4 at this operating
point. (The often-quoted factor 1/12 corresponds to Delta omega ~ gamma_2,
i.e. the edge of adiabatic validity, doubling R but coinciding with the
peak of the cascaded TPA Lorentzian -- not a safe operating point.)

For inversion-symmetric hosts (Si, Si_3 N_4) the cascade route is
structurally unavailable; we annotate them rather than draw a curve.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from params import ETA_GEOM, HOSTS

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def enhancement(Q2: float, host: dict, eta: float = ETA_GEOM) -> float:
    """R at Delta omega = 2 gamma_2 (factor 24 in the prefactor)."""
    return (Q2 / (24.0 * host["n"] ** 2)) * eta * host["chi2"] ** 2 / host["chi3"]


def Q2_crossover(host: dict, eta: float = ETA_GEOM) -> float:
    """Q_2 at which R = 1, at the same conservative operating point."""
    return 24.0 * host["n"] ** 2 * host["chi3"] / (eta * host["chi2"] ** 2)


def main():
    Q2_grid = np.geomspace(10, 1e6, 400)

    fig, ax = plt.subplots(figsize=(9.5, 5.6))

    # Achievable-Q_2 band (rough qBIC literature range)
    ax.axvspan(1e3, 1e5, alpha=0.07, color="grey", zorder=0,
               label=r"qBIC range, $Q_2 \in [10^3, 10^5]$")

    ax.axhline(1.0, color="k", lw=0.8, ls="--",
               label=r"$R = 1$: cascade $=$ native")

    cascade_hosts = [(n, h) for n, h in HOSTS.items() if h["cascade_available"]]
    inv_sym_hosts = [(n, h) for n, h in HOSTS.items() if not h["cascade_available"]]

    for name, host in cascade_hosts:
        R = enhancement(Q2_grid, host)
        Q2_x = Q2_crossover(host)
        label = (f"{name}: $|\\chi^{{(2)}}|={host['chi2']*1e12:.0f}$ pm/V, "
                 f"$\\chi^{{(3)}}_{{\\rm nat}}={host['chi3']*1e21:.2g}\\times 10^{{-21}}$")
        ax.loglog(Q2_grid, R, color=host["color"], lw=2.0, label=label)

        if Q2_grid[0] <= Q2_x <= Q2_grid[-1]:
            ax.plot([Q2_x], [1.0], "v", color=host["color"], ms=9, mec="k",
                    mew=0.8)
            x_text = Q2_x * 1.25 if Q2_x < 1e5 else Q2_x * 0.06
            ax.annotate(rf"$Q_2^{{\rm x}} \approx {Q2_x:.0f}$",
                        xy=(Q2_x, 1.0), xytext=(x_text, 0.2),
                        fontsize=9, color=host["color"])

    # Inversion-symmetric: annotation in the lower-right (data-free zone)
    if inv_sym_hosts:
        text = (r"Inversion-symmetric ($\chi^{(2)} \equiv 0$):" + "\n  "
                + ", ".join(n for n, _ in inv_sym_hosts)
                + "\n  cascade unavailable;" + "\n  native-Kerr only.")
        ax.text(0.97, 0.04, text, transform=ax.transAxes,
                fontsize=9, color="0.25", ha="right", va="bottom",
                bbox=dict(facecolor="white", edgecolor="0.7",
                          boxstyle="round,pad=0.35"))

    ax.set_xlabel(r"SH-cavity quality factor  $Q_2$")
    ax.set_ylabel(r"$R = |\chi^{(3)}_{\rm eff,casc}| \,/\, |\chi^{(3)}_{\rm native}|$")
    ax.set_title(r"Cascade-vs-native by host (note \S 11, table); $\Delta\omega = 2\gamma_2$")
    ax.set_xlim(Q2_grid[0], Q2_grid[-1])
    ax.set_ylim(1e-2, 1e5)
    ax.grid(True, which="both", alpha=0.3)
    # Move legend below the plot to free up the upper-left curves
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2,
              fontsize=8.5, framealpha=0.95)

    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "cascaded_vs_native"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))

    print(f"wrote {out}.{{png,pdf}}")
    print(f"  R at Q_2 = 4000 (note's design point):")
    for name, host in cascade_hosts:
        R4k = enhancement(4_000, host)
        Q2x = Q2_crossover(host)
        print(f"    {name:8s}: R = {R4k:7.2f}, Q_2_crossover = {Q2x:>9.0f}")
    plt.close(fig)


if __name__ == "__main__":
    main()
