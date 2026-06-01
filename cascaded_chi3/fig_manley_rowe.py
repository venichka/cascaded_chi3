"""Manley-Rowe / energy-conservation check.

Lossless TCMT (gamma_1 = gamma_2 = 0, no input port). Energy
|a_1|^2 + |a_2|^2 is conserved if and only if the SHG and DFG couplings
carry the same constant -- the factor-of-2 bookkeeping in the DFG
polarisation P_i^(2)(omega_1) = 2 epsilon_0 chi^(2) E(2omega_1) E*(omega_1)
that the note rederives in Section 2.

We compare two RHS:
  correct : da1/dt = +i conj(a1) a2,    da2/dt = +i a1^2
  wrong   : da1/dt = +i (1/2) conj(a1) a2,  da2/dt = +i a1^2
            (the original note's "kappa vs kappa/2" mismatch)

Initial state (a1, a2) = (1+0j, 0+0j); Domega = 0 (degenerate, on-resonance
SHG/DFG, fastest energy exchange). Plot |a1|^2, |a2|^2, and the total
vs time. Correct: total flat. Wrong: total drifts.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def make_rhs(dfg_factor: float):
    """Lossless RHS with adjustable DFG coupling. dfg_factor=1.0 is correct."""
    def rhs(t, y):
        a1 = y[0] + 1j * y[2]
        a2 = y[1] + 1j * y[3]
        da1 = 1j * dfg_factor * np.conj(a1) * a2  # DFG only, no losses
        da2 = 1j * a1 * a1                          # SHG, no losses
        return [da1.real, da2.real, da1.imag, da2.imag]
    return rhs


def integrate_lossless(dfg_factor, y0, t_eval):
    sol = solve_ivp(
        make_rhs(dfg_factor),
        (t_eval[0], t_eval[-1]), y0, t_eval=t_eval,
        method="DOP853", rtol=1e-12, atol=1e-14, max_step=0.05,
    )
    a1 = sol.y[0] + 1j * sol.y[2]
    a2 = sol.y[1] + 1j * sol.y[3]
    return a1, a2


def main():
    t_eval = np.linspace(0, 10.0, 4000)
    y0 = [1.0, 0.0, 0.0, 0.0]  # a1=1, a2=0

    a1_c, a2_c = integrate_lossless(1.0, y0, t_eval)
    a1_w, a2_w = integrate_lossless(0.5, y0, t_eval)

    n1_c, n2_c = np.abs(a1_c) ** 2, np.abs(a2_c) ** 2
    n1_w, n2_w = np.abs(a1_w) ** 2, np.abs(a2_w) ** 2

    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.4))

    axs[0].plot(t_eval, n1_c, color="C0", lw=1.5, label=r"$|a_1|^2$")
    axs[0].plot(t_eval, n2_c, color="C3", lw=1.5, label=r"$|a_2|^2$")
    axs[0].plot(t_eval, n1_c + n2_c, color="k", lw=2.0, ls="--",
                label=r"$|a_1|^2 + |a_2|^2$ (total)")
    axs[0].axhline(1.0, color="grey", lw=0.7, ls=":")
    axs[0].set_xlabel(r"$t \cdot \beta$")
    axs[0].set_ylabel("mode energy")
    axs[0].set_title(r"(a) Correct DFG coupling — total energy conserved")
    axs[0].legend(fontsize=9)
    axs[0].grid(True, alpha=0.3)
    axs[0].set_ylim(-0.05, 1.15)

    axs[1].plot(t_eval, n1_w, color="C0", lw=1.5, label=r"$|a_1|^2$")
    axs[1].plot(t_eval, n2_w, color="C3", lw=1.5, label=r"$|a_2|^2$")
    axs[1].plot(t_eval, n1_w + n2_w, color="k", lw=2.0, ls="--",
                label=r"$|a_1|^2 + |a_2|^2$ (total)")
    axs[1].axhline(1.0, color="grey", lw=0.7, ls=":")
    axs[1].set_xlabel(r"$t \cdot \beta$")
    axs[1].set_ylabel("mode energy")
    axs[1].set_title(r"(b) DFG halved — Manley$-$Rowe violated, total drifts")
    axs[1].legend(fontsize=9)
    axs[1].grid(True, alpha=0.3)
    axs[1].set_ylim(-0.05, 1.6)

    drift_correct = (n1_c + n2_c).max() - (n1_c + n2_c).min()
    drift_wrong = (n1_w + n2_w).max() - (n1_w + n2_w).min()
    print(f"  total-energy drift, correct DFG : {drift_correct:.2e}")
    print(f"  total-energy drift, halved DFG  : {drift_wrong:.2e}")

    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "manley_rowe"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
