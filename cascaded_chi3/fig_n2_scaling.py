"""1/n_2^2 scaling vs the naive 1/n^3.

The note's "Explicit n-dependence" paragraph (Section 5):

    correct (the note's eq. chi3casc_explicit):
        |chi^(3)_casc| ~ omega_1 / (6 n_2^2 |Domega|) * |chi^(2)|^2 * eta_geom

    naive (the original-derivation form, single bulk index n):
        |chi^(3)_casc|_naive ~ omega_1 / (6 n^3 |Domega|) * |chi^(2)|^2 * eta_geom

The four powers of n_1 in the numerator and denominator cancel exactly --
only n_2 survives, with two powers (the SH-mode normalisation). Setting
n = (n_1+n_2)/2 in the naive form, the discrepancy is
correct/naive = n^3 / n_2^2, typically ~2 for LiNbO3 / GaAs.

Plot |chi^(3)_casc| vs n_2 (with n_1 frozen at 2.21) for both formulas.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from params import (CHI2_REAL, ETA_GEOM, N1_REAL, OMEGA1_REAL, Q2_REF)

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def chi3_correct(n2, omega1, Domega, eta):
    """note eq. (chi3casc_explicit). Only n_2 enters."""
    return omega1 / (6.0 * n2 ** 2 * abs(Domega)) * CHI2_REAL ** 2 * eta


def chi3_naive(n_avg, omega1, Domega, eta):
    """Single-index plane-wave-style formula: 1/n^3."""
    return omega1 / (6.0 * n_avg ** 3 * abs(Domega)) * CHI2_REAL ** 2 * eta


def main():
    omega1 = OMEGA1_REAL
    gamma2 = 2 * omega1 / Q2_REF
    Domega = 2 * gamma2

    n2_grid = np.linspace(1.5, 4.0, 400)
    n_avg_grid = 0.5 * (N1_REAL + n2_grid)

    chi_corr = chi3_correct(n2_grid, omega1, Domega, ETA_GEOM)
    chi_naiv = chi3_naive(n_avg_grid, omega1, Domega, ETA_GEOM)

    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.5))

    axs[0].plot(n2_grid, chi_corr * 1e21, "-", color="C0", lw=2.0,
                label=r"correct (note): $\propto 1/n_2^2$")
    axs[0].plot(n2_grid, chi_naiv * 1e21, "--", color="C3", lw=2.0,
                label=r"naive plane-wave: $\propto 1/\bar{n}^3$, $\bar{n}=(n_1+n_2)/2$")
    axs[0].axvline(2.26, color="grey", lw=0.7, ls=":")
    axs[0].annotate(r"LiNbO$_3$, $n_2 = 2.26$",
                    xy=(2.26, 24), xytext=(2.95, 25),
                    fontsize=9.5, color="0.3",
                    arrowprops=dict(arrowstyle="-", color="0.5", lw=0.7))
    axs[0].set_xlabel(r"SH-mode effective index  $n_2$")
    axs[0].set_ylabel(r"$|\chi^{(3)}_{\rm casc}|$  $[10^{-21}\,\mathrm{m}^2/\mathrm{V}^2]$")
    axs[0].set_title(r"(a) Cascaded $\chi^{(3)}$ vs SH-mode index $n_2$")
    axs[0].grid(True, alpha=0.3)
    axs[0].legend(fontsize=9, loc="upper right")

    axs[1].plot(n2_grid, chi_corr / chi_naiv, "-", color="k", lw=2.0)
    axs[1].axvline(2.26, color="grey", lw=0.7, ls=":")
    ratio_at_LN = chi3_correct(2.26, omega1, Domega, ETA_GEOM) / chi3_naive(0.5*(N1_REAL+2.26), omega1, Domega, ETA_GEOM)
    axs[1].plot([2.26], [ratio_at_LN], "ko", ms=8)
    axs[1].annotate(f"LiNbO$_3$: {ratio_at_LN:.2f}$\\times$",
                    xy=(2.26, ratio_at_LN), xytext=(2.5, ratio_at_LN + 0.1),
                    fontsize=10)
    axs[1].axhline(1.0, color="grey", lw=0.6, ls=":")
    axs[1].set_xlabel(r"SH-mode index  $n_2$")
    axs[1].set_ylabel(r"correct / naive")
    axs[1].set_title(r"(b) Naive $1/n^3$ underestimates by $\sim$2$\times$ at telecom indices")
    axs[1].grid(True, alpha=0.3)

    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "n2_scaling"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    print(f"  LiNbO3 ratio (correct/naive) at n_2=2.26: {ratio_at_LN:.3f}")
    plt.close(fig)


if __name__ == "__main__":
    main()
