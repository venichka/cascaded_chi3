"""Parametric-oscillation onset above the saturation threshold.

Two complementary panels:

(a) Linear-stability calculation (derived, not asserted).
    For the trivial fixed point a_2 = i*a_1^2/(1/2 - i*Delta_SH) we
    build the 4x4 real Jacobian of the TCMT and read off the largest
    real part of its eigenvalues, vs the saturation parameter
    sat = |beta a_1|/|Delta_SH|. The crossing of zero gives the
    parametric-oscillation threshold sat_threshold. We do this for
    three FM-damping levels gamma_1/gamma_2 in {2, 5, 10} so the
    realistic LiNbO3 design point (gamma_1/gamma_2 = 10) is included
    alongside the lighter-damped cases.

(b) Time-domain validation. Drive CW at four input levels for the
    lightest-damped case (gamma_1/gamma_2 = 2, where the threshold
    sits at sat ~ 1.24); seed a 1e-3j perturbation off the trivial FP
    and integrate. Below threshold: decay back. Above threshold:
    exponential growth into a limit cycle.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from analytic import (a2_adiabatic, alpha3_eff_casc_full,
                      cw_kerr_only_intracavity, gamma_2pa_eff)
from linear_stability import find_threshold, stability_vs_sat
from params import TCMTParams
from tcmt import integrate

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def s_plus_for_sat(p: TCMTParams, target_sat: float) -> float:
    """Solve the Kerr-only cubic for |s+| that yields the requested sat."""
    target_n1 = (target_sat * abs(p.Domega)) ** 2
    a_e = alpha3_eff_casc_full(p.Domega)
    g_e = gamma_2pa_eff(p.Domega)
    g1 = 0.5 * p.gamma1
    bracket = (g1 + g_e * target_n1) ** 2 + (a_e * target_n1) ** 2
    return float(np.sqrt(target_n1 * bracket / p.gamma1_rad))


def main():
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 4.7))

    # ----- panel (a): linear-stability eigenvalue vs sat -----------------
    sat_grid = np.linspace(0.4, 3.0, 600)
    cases_lin = [
        (2.0,  "C3", r"$\gamma_1/\gamma_2 = 2$  (lighter damped)"),
        (5.0,  "C2", r"$\gamma_1/\gamma_2 = 5$  (intermediate)"),
        (10.0, "C0", r"$\gamma_1/\gamma_2 = 10$ (LiNbO$_3$ design pt.)"),
    ]
    thresholds = []
    for gamma1, color, label in cases_lin:
        p = TCMTParams(Domega=2.0, gamma1=gamma1,
                       gamma1_rad=0.5 * gamma1, alpha3=0.0, delta1=0.0)
        re_max = stability_vs_sat(p, sat_grid)
        axs[0].plot(sat_grid, re_max, color=color, lw=1.8, label=label)
        thr, _, _ = find_threshold(p, 0.4, 3.0, 600)
        if thr is not None:
            axs[0].axvline(thr, color=color, lw=0.8, ls=":", alpha=0.7)
            thresholds.append((thr, color, gamma1))

    # Stack threshold annotations vertically in a clear block (above the
    # axis) so they don't crowd the zero crossings.
    for k, (thr, color, gamma1) in enumerate(thresholds):
        axs[0].annotate(rf"$\gamma_1/\gamma_2={gamma1:.0f}$: "
                        rf"$\mathrm{{sat}}_\mathrm{{thr}} = {thr:.2f}$",
                        xy=(0.97, 0.97 - 0.07 * k),
                        xycoords="axes fraction", ha="right", va="top",
                        color=color, fontsize=9)

    axs[0].axhline(0.0, color="k", lw=0.6)
    axs[0].axvline(1.0, color="grey", lw=0.6, ls="--", alpha=0.6)
    axs[0].annotate(r"naive $|\beta a_1| = |\Delta_{\mathrm{SH}}|$",
                    xy=(1.0, -2.0), xytext=(1.5, -2.2),
                    color="0.45", fontsize=9,
                    arrowprops=dict(arrowstyle="-", color="0.6", lw=0.6))
    axs[0].set_xlabel(r"saturation parameter  $|\beta a_1|/|\Delta_{\mathrm{SH}}|$")
    axs[0].set_ylabel(r"$\max\,\mathrm{Re}(\lambda)$ of the 4$\times$4 Jacobian")
    axs[0].set_title("(a) Linear-stability threshold vs FM damping")
    axs[0].set_xlim(0.4, 3.0)
    axs[0].set_ylim(-2.5, 2.5)
    axs[0].grid(True, alpha=0.3)
    axs[0].legend(fontsize=9, loc="lower right")

    # ----- panel (b): time-domain validation, gamma_1/gamma_2 = 2 ---------
    p = TCMTParams(Domega=2.0, gamma1=2.0, gamma1_rad=1.0, alpha3=0.0,
                   delta1=0.0)
    thr_lin, _, _ = find_threshold(p, 0.4, 3.0, 600)

    cases_td = [
        (0.7, "C0", r"sat $= 0.7$ (below)"),
        (1.0, "C2", r"sat $= 1.0$ (at naive)"),
        (1.3, "C3", rf"sat $= 1.3$ (above stab. thr. $\approx {thr_lin:.2f}$)"),
        (1.7, "C1", r"sat $= 1.7$ (deep above)"),
    ]
    t_eval = np.linspace(0, 100.0, 4000)
    for target_sat, color, lab in cases_td:
        sp = s_plus_for_sat(p, target_sat)
        n1_triv = cw_kerr_only_intracavity(p, sp)[0]
        a1_triv = np.sqrt(n1_triv) + 0j
        a2_triv = a2_adiabatic(a1_triv, p.Domega)
        seed = (a1_triv, a2_triv + 1e-3j)
        sol = integrate(p, sp, (0, t_eval[-1]), y0=seed, t_eval=t_eval,
                        max_step=0.05)
        n2_t = np.abs(sol.a[1]) ** 2
        axs[1].semilogy(t_eval, n2_t, color=color, lw=1.5, label=lab)
        axs[1].axhline(abs(a2_triv) ** 2, color=color, lw=0.6, ls=":",
                       alpha=0.7)

    axs[1].set_xlabel(r"$t \cdot \gamma_2$")
    axs[1].set_ylabel(r"$|a_2|^2$ (intracavity SH energy)")
    axs[1].set_title(rf"(b) Time-domain validation, $\gamma_1/\gamma_2 = 2$")
    # Legend placed outside upper-right to avoid the limit-cycle data
    axs[1].legend(loc="upper left", bbox_to_anchor=(0.01, 0.98),
                  fontsize=9, framealpha=0.95)
    axs[1].grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "parametric_onset"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    print(f"  linear-stability thresholds:")
    for gamma1, _, _ in cases_lin:
        p = TCMTParams(Domega=2.0, gamma1=gamma1,
                       gamma1_rad=0.5 * gamma1, alpha3=0.0, delta1=0.0)
        thr, _, _ = find_threshold(p, 0.4, 3.0, 600)
        print(f"    gamma1/gamma2 = {gamma1}: sat_threshold = {thr:.3f}")
    plt.close(fig)


if __name__ == "__main__":
    main()
