"""Cascaded enhancement vs the dimensionless tensor overlap eta_geom.

From eq. (chi3casc_explicit) of the note:

    |chi^(3)_casc| ~ (omega_1 / (6 n_2^2 |Domega|)) * |chi^(2)|^2 * eta_geom

so the enhancement over the native chi^(3) scales linearly with eta_geom.
At fixed Q_2 and design wavelength,

    enhancement = (Q_2 / (12 n_2^2)) * eta_geom * |chi^(2)|^2 / |chi^(3)|_native.

Section 10(5) of the note: random mode pairings give eta ~ 1e-3, while a
properly engineered qBIC pair gives ~0.3. The geometric overlap is the
binding constraint, not Q_2.

Plot enhancement vs eta_geom on log-log for two Q_2 values, with
horizontal markers for "= native chi^(3)" and shaded regimes.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from params import CHI2_REAL, N2_REAL

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"

CHI3_NATIVE = 2e-21


def enhancement(Q2, eta):
    return (Q2 / (12.0 * N2_REAL ** 2)) * eta * CHI2_REAL ** 2 / CHI3_NATIVE


def main():
    eta_grid = np.geomspace(1e-4, 1.0, 400)

    fig, ax = plt.subplots(figsize=(8.5, 5.2))

    for Q2, color in [(4_000, "C0"), (10_000, "C3")]:
        enh = enhancement(Q2, eta_grid)
        ax.loglog(eta_grid, enh, color=color, lw=2.0, label=f"$Q_2 = {Q2}$")

    # Mark regimes (Section 10(5) of the note).
    ax.axvspan(1e-4, 1e-2, color="C3", alpha=0.07,
               label=r"random mode pairing ($\eta \sim 10^{-3}$)")
    ax.axvspan(0.1, 0.6, color="C2", alpha=0.10,
               label=r"engineered qBIC pair ($\eta \sim 0.3$)")
    ax.axvline(1.0, color="k", lw=0.6, ls=":")
    ax.text(1.02, 4e-1, "perfect overlap", fontsize=9, color="k")

    ax.axhline(1.0, color="k", lw=0.7, ls="--")
    ax.text(1.5e-4, 1.2,
            r"$|\chi^{(3)}_{\rm casc}| = |\chi^{(3)}_{\rm native}|$",
            fontsize=9, color="k")

    ax.set_xlabel(r"tensor overlap  $\eta_{\rm geom}$")
    ax.set_ylabel(r"$|\chi^{(3)}_{\rm casc}|/|\chi^{(3)}_{\rm native}|$")
    ax.set_title(r"Geometric overlap is the binding constraint, not $Q_2$")
    ax.set_xlim(1e-4, 1.5)
    ax.set_ylim(1e-2, 5e2)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right", fontsize=9.5)

    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "overlap_eta"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    print(f"  Q2=4000, eta=0.3: enh = {enhancement(4000, 0.3):.1f}")
    print(f"  Q2=4000, eta=1e-3: enh = {enhancement(4000, 1e-3):.3f}")
    print(f"  Q2=10000, eta=0.3: enh = {enhancement(10000, 0.3):.1f}")
    plt.close(fig)


if __name__ == "__main__":
    main()
