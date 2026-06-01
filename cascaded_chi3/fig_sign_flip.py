"""Sign tunability: sweep Domega through 0, watch alpha3_eff change sign.

Two curves:
  - alpha3_bare = 0 (pure cascading): zero crossing exactly at Domega = 0.
  - alpha3_bare > 0 (e.g. 0.05): zero shifts to Domega = 1/alpha3_bare.

The grey shaded band |Domega| < 1 marks the regime where adiabatic
elimination starts to break down (flagged by the validation figure).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from analytic import alpha3_eff_casc_full
from params import TCMTParams
from tcmt import steady_state

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def alpha_eff_numerical(p, s_plus=0.01):
    """Total effective Kerr (cascaded + bare) extracted from steady state.

    The cascaded part comes from i*conj(a1)*a2 / (a1 * |a1|^2) (the
    imaginary part). a_2 has no alpha_3 dependence, so the bare Kerr
    must be added explicitly to compare with the analytic full curve.
    """
    a1, a2 = steady_state(p, s_plus, t_settle=200.0 / min(p.gamma1, 1.0))
    C = 1j * np.conj(a1) * a2 / (a1 * abs(a1) ** 2)
    return p.alpha3 + C.imag


def main():
    # Plot the full Lorentzian form straight through Domega = 0;
    # the cascaded Kerr passes smoothly through zero (it does NOT
    # diverge -- that's the "$-1/Domega$" approximation, not the full
    # form). The cost of operating near zero is captured by the
    # cascaded TPA gamma_2pa_eff = (gamma_2/2)/(Domega^2 + (gamma_2/2)^2),
    # which peaks at Domega = 0; we annotate that fact.
    Dw = np.linspace(-30, 30, 600)

    fig, ax = plt.subplots(figsize=(8, 4.6))

    for alpha3, color, label in [
        (0.00, "C0", r"pure cascading ($\widetilde{\alpha}_3=0$)"),
        (0.05, "C3", r"with bare Kerr $\widetilde{\alpha}_3=+0.05$"),
    ]:
        # analytic: bare + cascaded
        alpha_an = alpha3 + alpha3_eff_casc_full(Dw)
        ax.plot(Dw, alpha_an, "-", color=color, lw=2, label=label + " (analytic)")

        # numerical: thin out for speed; skip a tiny window |Domega| < 0.6
        # where the FM steady state is dominated by the cascaded TPA and
        # the small-signal extraction has poor SNR.
        Dw_num_full = Dw[::20]
        Dw_num = Dw_num_full[np.abs(Dw_num_full) > 0.6]
        alpha_num = np.array([
            alpha_eff_numerical(TCMTParams(Domega=d, gamma1=10.0,
                                           gamma1_rad=5.0, alpha3=alpha3))
            for d in Dw_num
        ])
        ax.plot(Dw_num, alpha_num, "o",
                color=color, ms=5, mfc="white", mew=1.3,
                label=label + " (TCMT)")

        # mark the zero crossing of the analytic curve
        if alpha3 != 0:
            zero_at = 1.0 / alpha3
            if min(Dw) < zero_at < max(Dw):
                ax.plot([zero_at], [0], "v", color=color, ms=8, mfc=color)
                ax.annotate(rf"zero @ $\Delta\omega/\gamma_2={zero_at:.0f}$",
                            xy=(zero_at, 0), xytext=(zero_at - 8, 0.04),
                            fontsize=9, color=color)

    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(0, color="k", lw=0.5)
    # Annotate the TPA-dominated regime: the Kerr passes through 0,
    # but the imaginary part (cascaded TPA) peaks here.
    ax.annotate(
        "cascaded Kerr $\\to 0$,\n cascaded TPA peaks here",
        xy=(0, 0), xytext=(3, -0.55),
        fontsize=9, color="0.35",
        arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8),
    )
    ax.set_xlabel(r"$\Delta\omega/\gamma_2 = (2\omega_1-\omega_2)/\gamma_2$")
    ax.set_ylabel(r"$\widetilde{\alpha}_3^{\rm eff}$")
    ax.set_title("Sign tunability: walking the SH resonance through twice the FM")
    ax.set_ylim(-1.25, 1.25)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.95)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "sign_flip"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
