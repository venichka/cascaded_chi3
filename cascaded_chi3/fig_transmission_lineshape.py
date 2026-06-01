"""Power-resolved complex transmission t(omega_p) at the FM port.

This is the most directly measurable signature of the cascaded chi^(3):
sweep the pump frequency across the FM resonance at several input
powers and watch the lineshape distort.

Output is the FM-port transmission |s_-/s_+|^2 with
  s_- = -s_+ + sqrt(gamma_1_rad) * a_1     (Suh-Wang-Fan).

Realistic LiNbO3 design point (Section 9 of the note):
  Domega/gamma_2 = 2,  gamma_1/gamma_2 = 10,  gamma_1_rad/gamma_1 = 0.5.
With Domega > 0 the cascaded Kerr is negative (alpha_eff_casc < 0): the
effective FM resonance shifts DOWN with intensity, so the dip walks to
LOWER pump frequency and develops a steep edge on the red side.

At low power: clean Lorentzian dip centred at delta_1 = 0.
At intermediate power: dip pulled red, asymmetric.
Past the cusp: hysteresis -- the up-sweep (delta_1 from blue to red) and
down-sweep traces differ. The figure shows both directions for the
highest-power curve.

We use the analytic Kerr+2PA cubic root for the up- and down-sweep
selection (continuation), bypassing the integrator's slow drift through
the saddle-node region while staying physically equivalent to the TCMT.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from analytic import cw_kerr_only_intracavity
from params import TCMTParams

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def transmission_branches(p: TCMTParams, s_plus: float):
    """Return (low, mid, high) intracavity-energy branches for given delta_1.
    Each is np.nan if the corresponding root doesn't exist."""
    roots = cw_kerr_only_intracavity(p, s_plus)
    if roots.size == 1:
        return roots[0], np.nan, np.nan
    if roots.size >= 3:
        rs = sorted(roots)
        return rs[0], rs[1], rs[2]
    if roots.size == 2:
        rs = sorted(roots)
        return rs[0], np.nan, rs[1]
    return np.nan, np.nan, np.nan


def transmission_from_n1(p, s_plus, n1):
    """|s-/s+|^2 given a chosen intracavity-energy branch root n1.

    From the FM steady state: a_1 (gamma1/2 + gamma_eff_2pa n1 - i (delta1 + alpha_eff n1))
                            = sqrt(gamma1_rad) s_plus
    -> a_1 = sqrt(gamma1_rad) s_plus / [gamma1/2 + g_e n1 - i (delta1 + a_e n1)]
    """
    from analytic import alpha3_eff_casc_full, gamma_2pa_eff
    Delta_SH = 2 * p.delta1 + p.Domega
    a_e = p.alpha3 + alpha3_eff_casc_full(Delta_SH)
    g_e = gamma_2pa_eff(Delta_SH)
    g1 = 0.5 * p.gamma1
    if not np.isfinite(n1):
        return np.nan
    denom = (g1 + g_e * n1) - 1j * (p.delta1 + a_e * n1)
    a1 = np.sqrt(p.gamma1_rad) * s_plus / denom
    s_minus = -s_plus + np.sqrt(p.gamma1_rad) * a1
    return abs(s_minus / s_plus) ** 2


def main():
    base = dict(Domega=2.0, gamma1=10.0, gamma1_rad=5.0, alpha3=0.0)
    delta1_grid = np.linspace(-15.0, 15.0, 800)

    # Three input-power levels.
    s2_levels = [(0.5, "low (small-signal)", "C0"),
                 (10.0, "moderate", "C2"),
                 (60.0, "high (near cusp)", "C3")]

    fig, ax = plt.subplots(figsize=(8.5, 5.0))

    for s2, label, color in s2_levels:
        s_plus = np.sqrt(s2)
        T_low = np.empty_like(delta1_grid)
        T_high = np.empty_like(delta1_grid)
        for i, d1 in enumerate(delta1_grid):
            p = TCMTParams(delta1=float(d1), **base)
            n1_lo, n1_mid, n1_hi = transmission_branches(p, s_plus)
            T_low[i] = transmission_from_n1(p, s_plus, n1_lo)
            T_high[i] = transmission_from_n1(p, s_plus, n1_hi)

        # If the high branch exists somewhere, draw it as the "down sweep"
        # (you arrive there by coming from large amplitude at high blue
        # detuning); the low branch is the "up sweep".
        ax.plot(delta1_grid, T_low, color=color, lw=1.7, label=label)
        if np.any(np.isfinite(T_high)):
            mask = np.isfinite(T_high)
            ax.plot(delta1_grid[mask], T_high[mask], color=color, lw=1.7,
                    ls="--", alpha=0.8)

    ax.axhline(1.0, color="k", lw=0.5, ls=":")
    ax.axvline(0.0, color="grey", lw=0.7, ls=":")
    ax.set_xlabel(r"pump detuning  $\delta_1 / \gamma_2 = (\omega_p - \omega_1)/\gamma_2$")
    ax.set_ylabel(r"transmission  $|s_-/s_+|^2$")
    ax.set_title(
        "FM-port transmission vs pump frequency\n"
        r"(solid = lower-branch root; dashed = upper-branch root, only inside the bistable region)",
        fontsize=11,
    )
    ax.set_ylim(0.0, 1.05)
    ax.set_xlim(-15, 15)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", title="input power")

    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "transmission_lineshape"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
