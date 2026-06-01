"""Cascade as a passive optical limiter.

Tier 2, item 7 of the note's Section 13. The cascade saturation
|beta a_1| -> |Domega| is a hard, passive intensity ceiling at
intracavity energy |a_1|^2_cap = (Domega/beta)^2. Above this threshold
the trivial fixed point destabilises and the FM clamps. In an
integrated metasurface this is an all-optical power limiter, with the
twist that the ceiling is tunable via Domega.

Native-Kerr resonance has no such mechanism: intracavity intensity
scales linearly with input up to material damage.

Two panels:
  (a) Transmitted power |s_-|^2 vs incident |s_+|^2 for cascade
      (three Domega values) and a linear-cavity reference.
  (b) The saturated transmission level |s_-|^2_cap vs Domega/gamma_2:
      decreasing Domega lowers the cap (tighter limiter).
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
from tcmt import attractor


HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def transmitted_power_full(p: TCMTParams, s_p: float):
    """Return (<|a_1|^2>, <|s_-/s_+|^2>, is_limit_cycle) from full-TCMT
    attractor-averaging.  Works for stable fixed point and limit cycle alike."""
    res = attractor(p, complex(s_p), t_settle=200.0, t_avg=200.0,
                    n_samples=800)
    return res["n1_mean"], res["T_mean"], res["is_limit_cycle"]


def main():
    # Same FM cavity in all cases; sweep Domega for the cascade
    base = dict(gamma1=2.0, gamma1_rad=1.0, alpha3=0.0, delta1=0.0)

    Domega_values = [1.5, 2.5, 4.0]
    colors = ["C3", "C0", "C2"]
    s2_grid = np.geomspace(1e-3, 100.0, 32)

    fig, axs = plt.subplots(1, 2, figsize=(13.0, 5.0),
                            gridspec_kw={"wspace": 0.32})

    # Linear cavity reference (no NL): |a_1|^2 = gamma_rad |s_+|^2 / (gamma_1/2)^2
    g1 = 0.5 * base["gamma1"]
    n1_linear = base["gamma1_rad"] * s2_grid / g1 ** 2

    # ---- Panel (a): intracavity energy vs incident, three Domega values ----
    # Full TCMT attractor-averaging at each pump.  Points where |beta a_1| >
    # |Delta omega| are past the parametric-stability boundary and are shown
    # as open markers (limit-cycle time average).
    for Domega, color in zip(Domega_values, colors):
        p = TCMTParams(Domega=Domega, beta=1.0, **base)
        n1 = np.empty_like(s2_grid)
        for i, s2 in enumerate(s2_grid):
            n1[i], _, _ = transmitted_power_full(p, np.sqrt(s2))

        cap = Domega ** 2  # |beta a_1|^2 = |Delta omega|^2 ceiling
        unstable = n1 > cap   # |beta a_1|^2 > |Delta omega|^2 marks past-cap
        # Stable points (filled, solid)
        axs[0].loglog(s2_grid[~unstable], n1[~unstable], "o-",
                      color=color, ms=5, mfc=color, mew=0, lw=1.7,
                      label=rf"cascade, $\Delta\omega/\gamma_2 = {Domega:g}$"
                            rf"  $\Rightarrow$ cap $= {cap:g}$")
        # Past-cap points (open, dashed)
        if unstable.any():
            axs[0].loglog(s2_grid[unstable], n1[unstable], "s--",
                          color=color, ms=5, mfc="white", mew=1.2,
                          lw=1.0, alpha=0.7)
        axs[0].axhline(cap, color=color, lw=0.8, ls=":", alpha=0.5)

    # Linear reference
    axs[0].loglog(s2_grid, n1_linear, "-", color="0.55", lw=2.0,
                  label=r"linear cavity (no NL): $\propto |s_+|^2$")

    axs[0].set_xlabel(r"incident power  $|s_+|^2$")
    axs[0].set_ylabel(r"intracavity FM energy  $\langle|a_1|^2\rangle_t$")
    axs[0].set_title(r"(a) Cascade clamps $|a_1|^2$; cap tunable via $\Delta\omega$")
    axs[0].grid(True, which="both", alpha=0.3)
    axs[0].legend(loc="lower right", fontsize=9)
    axs[0].text(0.03, 0.97,
                "filled: stable fixed point\n"
                r"open: past $|\beta a_1|>|\Delta\omega|$" + "\n"
                "   limit-cycle time average",
                transform=axs[0].transAxes,
                fontsize=8.5, color="0.20", ha="left", va="top",
                bbox=dict(facecolor="white", edgecolor="0.6",
                          boxstyle="round,pad=0.35"))

    # ---- Panel (b): saturated transmission vs Domega ---------------------
    # Three curves:
    #   (i)  Full-TCMT attractor-averaged <|a_1|^2> at |s+|^2 = 50 (markers)
    #   (ii) Self-consistent cubic (cw_kerr_only_intracavity) at the same pump:
    #        the lower stable root.  This is the right analytic prediction
    #        for the steady state below parametric threshold.
    #   (iii) Cap |Delta omega/beta|^2: the parametric-stability THRESHOLD,
    #        not the actual saturated value.  Asymptote / upper bound only.
    Domega_axis = np.linspace(0.5, 8.0, 60)
    n1_cap = np.zeros(Domega_axis.size)        # TCMT, attractor-averaged
    is_lc = np.zeros_like(Domega_axis, dtype=bool)
    n1_cubic = np.full_like(Domega_axis, np.nan)   # analytic cubic (lower root)
    s_input_high = 50.0
    for i, D in enumerate(Domega_axis):
        p = TCMTParams(Domega=D, beta=1.0, **base)
        n1_cap[i], _, is_lc[i] = transmitted_power_full(p, np.sqrt(s_input_high))
        roots = cw_kerr_only_intracavity(p, np.sqrt(s_input_high))
        if roots.size > 0:
            n1_cubic[i] = roots[0]   # lower stable root

    # Use the cubic to decide which TCMT points are "stable" (full TCMT
    # tracks the cubic) vs "past parametric threshold" (TCMT diverges,
    # limit-cycle time-average instead).
    stable = n1_cap < Domega_axis ** 2

    # Cubic analytic line (the right prediction below threshold)
    axs[1].plot(Domega_axis, n1_cubic, "-", color="k", lw=1.8,
                label=r"analytic cubic  $|a_1|^2_{\rm cubic}(|s_+|^2{=}50)$")
    # Cap asymptote (light dotted, in the background)
    axs[1].plot(Domega_axis, Domega_axis ** 2, ":", color="C3", lw=1.4,
                label=r"cap  $|\Delta\omega/\beta|^2$ (parametric threshold)")
    # TCMT markers
    axs[1].plot(Domega_axis[stable], n1_cap[stable], "o", color="C0",
                ms=6, mfc="C0", mew=0,
                label=r"TCMT $\langle|a_1|^2\rangle_t$  (stable FP)")
    if (~stable).any():
        axs[1].plot(Domega_axis[~stable], n1_cap[~stable], "s",
                    color="C0", ms=6, mfc="white", mew=1.3, alpha=0.7,
                    label=r"TCMT past cap (limit-cycle avg)")

    axs[1].set_xlabel(r"SH detuning  $\Delta\omega/\gamma_2$")
    axs[1].set_ylabel(r"saturated  $\langle|a_1|^2\rangle_t$  at  $|s_+|^2 = 50$")
    axs[1].set_title(r"(b) Saturated FM energy vs $\Delta\omega$ at high pump")
    axs[1].grid(True, which="both", alpha=0.3)
    axs[1].legend(loc="upper left", fontsize=8.5)

    fig.suptitle(
        r"Cascade as a passive optical limiter: intensity-clamped at "
        r"$|a_1|^2_{\rm cap} = |\Delta\omega/\beta|^2$",
        fontsize=11.5, y=1.00,
    )
    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "optical_limiter"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
