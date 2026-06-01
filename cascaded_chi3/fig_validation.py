"""Validation: numerical effective coefficients from steady-state TCMT vs analytic.

Pump on FM resonance (delta1=0), small |s+|. The cascaded coupling on the
FM line, in steady state, factorises as

    i * conj(a1) * a2 = (i*alpha_eff_casc + r) * |a1|^2 * a1
        ->  C := i * conj(a1) * a2 / (a1 * |a1|^2)
            alpha_eff_casc = Im(C),       (Kerr part)
            gamma_2pa_eff  = -Re(C)       (cascaded TPA, positive)

Two panels:
  (a) Cascaded Kerr alpha_eff_casc(Domega): odd in Domega, sign-tunable.
      TCMT numerics overlaid with the full analytic and the simple
      |Domega| >> gamma2/2 form -1/Domega; the simple form deviates near
      Domega = gamma2/2 (the adiabatic-elimination breakdown).
  (b) Cascaded two-photon-like loss gamma_2pa_eff(Domega): even in
      Domega, Lorentzian peak at Domega = 0 with FWHM gamma2 (the
      analogue of two-photon absorption, set by gamma2 rather than
      Im chi^(3); see Section 5(ii) of the note).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from analytic import (alpha3_eff_casc_full, alpha3_eff_casc_simple,
                      gamma_2pa_eff)
from params import TCMTParams
from tcmt import steady_state

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def main():
    # Sweep both signs to show sign tunability of the Kerr part and the
    # even symmetry of the TPA part.
    Dw_neg = -np.geomspace(0.6, 30.0, 35)[::-1]
    Dw_pos = np.geomspace(0.6, 30.0, 35)
    Dw_arr = np.concatenate([Dw_neg, Dw_pos])

    s_plus = 0.01
    base = dict(gamma1=10.0, gamma1_rad=5.0, alpha3=0.0, delta1=0.0)

    alpha_num = np.empty_like(Dw_arr)
    gamma_num = np.empty_like(Dw_arr)
    for i, Dw in enumerate(Dw_arr):
        p = TCMTParams(Domega=float(Dw), **base)
        a1, a2 = steady_state(p, s_plus, t_settle=200.0 / min(p.gamma1, 1.0))
        C = 1j * np.conj(a1) * a2 / (a1 * abs(a1) ** 2)
        alpha_num[i] = C.imag
        gamma_num[i] = -C.real

    Dw_fine_neg = np.linspace(-30, -0.05, 600)
    Dw_fine_pos = np.linspace(0.05, 30, 600)
    Dw_fine_full = np.linspace(-30, 30, 1200)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    # Panel (a): Kerr coefficient
    ax1.plot(Dw_fine_full, alpha3_eff_casc_full(Dw_fine_full), "-", color="k",
             lw=1.6, label=r"full: $-\Delta\omega/(\Delta\omega^2+1/4)$")
    for seg in (Dw_fine_neg, Dw_fine_pos):
        ax1.plot(seg, alpha3_eff_casc_simple(seg), "--", color="grey", lw=1.2,
                 label=(r"simple: $-1/\Delta\omega$" if seg is Dw_fine_neg else None))
    ax1.plot(Dw_arr, alpha_num, "o", color="C0", ms=4, mfc="white", mew=1.2,
             label="TCMT (small-signal)")
    ax1.axhline(0, color="k", lw=0.5)
    ax1.axvline(0, color="k", lw=0.5)
    ax1.set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax1.set_ylabel(r"$\widetilde{\alpha}_3^{\rm eff,casc}$")
    ax1.set_title("(a) Cascaded Kerr (real part) — sign-tunable through $\\Delta\\omega = 0$")
    ax1.set_ylim(-1.3, 1.3)
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Panel (b): TPA-like loss
    ax2.plot(Dw_fine_full, gamma_2pa_eff(Dw_fine_full), "-", color="k", lw=1.6,
             label=r"analytic: $(1/2)/(\Delta\omega^2+1/4)$")
    ax2.plot(Dw_arr, gamma_num, "s", color="C3", ms=4, mfc="white", mew=1.2,
             label="TCMT (small-signal)")
    ax2.axhline(0, color="k", lw=0.5)
    ax2.axvline(0, color="k", lw=0.5)
    ax2.set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax2.set_ylabel(r"$\widetilde{\gamma}_{2\rm PA}^{\rm eff}$")
    ax2.set_title(r"(b) Cascaded two-photon-like loss — Lorentzian, peak at $\Delta\omega = 0$")
    ax2.legend(loc="upper right", fontsize=9)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "validation"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
