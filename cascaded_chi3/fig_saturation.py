"""Saturation breakdown: as |beta*a1| approaches |Domega| the SH energy
becomes comparable to the FM energy and the small-signal "weakly cascaded"
description breaks down -- pump depletion appears.

The full TCMT and the Kerr-only adiabatic cubic agree exactly at any
deterministic steady state (the substitution a_2 = i*a_1^2/(1/2 - i*Domega)
is algebraically exact). HOWEVER, that trivial fixed point becomes
PARAMETRICALLY UNSTABLE above the parametric-oscillation threshold (when
|beta*a_1| ~ |Domega|, equivalently sat parameter ~ 1). Above threshold
the deterministic ODE enters a limit cycle and there is no longer a
single steady-state value to plot. We therefore stop the sweep just past
threshold (sat parameter ~ 1.1, where TCMT is still convergent) and
overlay the Kerr-only cubic to verify exact agreement everywhere shown.

Drive on FM resonance (delta1 = 0). Two panels:
  left  : |a1|^2 vs |s+|^2 -- full TCMT vs Kerr-only cubic vs a linear
          (no-cascading) cavity. The TCMT curve rolls below the linear
          one as the cascaded 2PA loss drains FM energy into the SH;
          it tracks the cubic to machine precision.
  right : intracavity energy fraction |a2|^2 / (|a1|^2 + |a2|^2) plus
          saturation parameter |beta*a1|/|Domega| on a twin axis.
          The fraction climbs toward 0.5 as sat parameter approaches 1.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from analytic import cw_kerr_only_intracavity, a2_adiabatic
from params import TCMTParams
from tcmt import integrate

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def main():
    # Realistic LiNbO3 qBIC metasurface design point from Section 9 of the note:
    #   Domega/gamma_2 = 2     (note's chosen operating point)
    #   gamma_1/gamma_2 = 10   (Q_1 ~ Q_2/10 = 400 -- broader, in-coupled FM)
    #   gamma_1_rad/gamma_1 = 0.5  (qBIC with controlled radiative coupling)
    p = TCMTParams(Domega=2.0, gamma1=10.0, gamma1_rad=5.0, alpha3=0.0,
                   delta1=0.0)

    # Linear cavity reference (cascading turned off): same gamma1 etc., no SHG.
    # |a1_lin|^2 = gamma1_rad * |s+|^2 / (gamma1/2)^2.
    g1 = 0.5 * p.gamma1

    # Stop just past the parametric-oscillation threshold. Below threshold
    # the TCMT trivial fixed point matches the Kerr-only cubic exactly
    # (we verify this on the plot). Above threshold the trivial fixed
    # point becomes unstable and the deterministic ODE enters a limit
    # cycle -- the integrator returns a time-averaged value that is NOT
    # the steady state and the comparison stops being meaningful.
    # Cap input below the parametric threshold (|a1|^2 = |Domega|^2 = 4
    # at sat parameter = 1). |s+|^2 = 20 gives peak sat ~ 0.9 with the
    # rolled-over shape clearly visible on a linear x-axis.
    s2 = np.geomspace(1e-3, 20.0, 70)

    n1_full = np.empty_like(s2)
    n2_full = np.empty_like(s2)
    n1_lin = p.gamma1_rad * s2 / g1**2

    state = (0 + 0j, 0 + 0j)
    t_settle = 1500.0 / min(p.gamma1, 1.0)
    for i, sq in enumerate(s2):
        sp = np.sqrt(sq)
        sol = integrate(p, sp, (0.0, t_settle), y0=state, max_step=2.0)
        a1, a2 = sol.a[0, -1], sol.a[1, -1]
        n1_full[i] = abs(a1) ** 2
        n2_full[i] = abs(a2) ** 2
        state = (a1, a2)

    # Analytic Kerr-only cubic: smallest positive root of the cubic in |a1|^2.
    n1_cubic = np.empty_like(s2)
    n2_cubic = np.empty_like(s2)
    for i, sq in enumerate(s2):
        roots = cw_kerr_only_intracavity(p, np.sqrt(sq))
        n1_cubic[i] = roots[0] if roots.size else np.nan
        # corresponding SH from the algebraic relation a2 = i*a1^2/(1/2-i*Domega)
        a1c = np.sqrt(n1_cubic[i])
        n2_cubic[i] = abs(a2_adiabatic(a1c, p.Domega)) ** 2

    sat = np.sqrt(n1_full) / abs(p.Domega)
    e_frac = n2_full / (n1_full + n2_full)

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))

    axs[0].loglog(s2, n1_lin, "-", color="grey", lw=2,
                  label="linear cavity (no cascading)")
    axs[0].loglog(s2, n1_cubic, "-", color="k", lw=1.2, alpha=0.6,
                  label="Kerr-only cubic (analytic SS)")
    axs[0].loglog(s2, n2_cubic, "-", color="k", lw=1.2, alpha=0.3)
    axs[0].loglog(s2, n1_full, "o", color="C0", ms=3.5, mfc="none", mew=1.0,
                  label=r"full TCMT $|a_1|^2$")
    axs[0].loglog(s2, n2_full, "s", color="C3", ms=3.5, mfc="none", mew=1.0,
                  label=r"full TCMT $|a_2|^2$")
    axs[0].axhline(p.Domega**2, color="C2", lw=0.8, ls=":",
                   label=r"$|\beta a_1|^2 = |\Delta\omega|^2$ (param. threshold)")
    axs[0].set_xlabel(r"$|s_+|^2$  (input power)")
    axs[0].set_ylabel(r"intracavity energy")
    axs[0].set_title("(a) Pump depletion: FM energy diverted into SH")
    axs[0].grid(True, which="both", alpha=0.3)
    axs[0].legend(fontsize=9, loc="lower right")

    line_e, = axs[1].plot(s2, e_frac, "o-", color="C3", ms=3.5,
                           label=r"SH energy fraction $|a_2|^2/(|a_1|^2+|a_2|^2)$")
    axs[1].axhline(0.5, color="C3", lw=1.0, ls="--", alpha=0.5)
    # Linear x-axis (as in the reference) makes the saturating shape visible;
    # log x makes the same curve look like a steady climb across decades.
    axs[1].set_xlabel(r"$|s_+|^2$  (input power)")
    axs[1].set_ylabel("SH energy fraction", color="C3")
    axs[1].tick_params(axis="y", labelcolor="C3")
    axs[1].set_ylim(0, 0.6)
    axs[1].grid(True, which="both", alpha=0.3)
    axs[1].set_title("(b) SH grows; saturation parameter approaches unity")

    ax2 = axs[1].twinx()
    line_s, = ax2.plot(s2, sat, "^-", color="C2", ms=3.5,
                       label=r"$|\beta a_1|/|\Delta\omega|$")
    ax2.axhline(1.0, color="C2", lw=0.8, ls=":", alpha=0.7)
    ax2.set_ylabel(r"saturation parameter  $|\beta a_1|/|\Delta\omega|$",
                   color="C2")
    ax2.tick_params(axis="y", labelcolor="C2")
    ax2.set_ylim(0, 1.3)

    axs[1].legend([line_e, line_s],
                  [line_e.get_label(), line_s.get_label()],
                  loc="upper left", fontsize=9)

    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "saturation"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
