"""Loss budget: cascade-specific TPA penalty (Caveat #1 of Section 12).

The cascaded coefficient on |a_1|^2 a_1 in the FM equation has both
real (Kerr) and imaginary (loss) parts:

    alpha_eff_casc = -|beta|^2 Domega / (Domega^2 + (gamma_2/2)^2)     [Kerr]
    gamma_eff_2pa  =  |beta|^2 (gamma_2/2) / (Domega^2 + (gamma_2/2)^2) [loss]

Ratio of Kerr to cascaded TPA:
    |alpha_eff_kerr| / gamma_eff_2pa  =  2 Domega / gamma_2.

So at the bandwidth-optimal point Domega = gamma_2/2 the Kerr-to-loss
ratio is unity (terrible contrast). Even at Domega = 2 gamma_2 (note's
nominal design point) it's only 4. Native-Kerr in a transparent host
like LiNbO_3 has essentially infinite contrast at telecom (material
absorption ~ 10^4 below the Kerr coefficient). Cascade therefore
carries an unavoidable contrast ceiling, even before considering
saturation. This is the central "price" of cascading.

Two panels:
  (a) |alpha_eff_kerr|/|beta|^2 and gamma_eff_2pa/|beta|^2 vs Domega/gamma_2;
      Kerr falls as 1/Domega for large Domega, TPA falls as 1/Domega^2.
  (b) Switching contrast |alpha_eff_kerr| / gamma_eff_2pa vs Domega/gamma_2
      for cascade; native horizontal line at material-absorption limit.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np


HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"

# Material-absorption-limited contrast for transparent native host
# (LiNbO_3 at 1550 nm: |chi^(3)|/Im chi^(3) ~ 10^4 typical).
NATIVE_CONTRAST_FLOOR = 1.0e4


def alpha_kerr(Dw):
    return Dw / (Dw ** 2 + 0.25)


def gamma_tpa(Dw):
    return 0.5 / (Dw ** 2 + 0.25)


def main():
    Dw = np.geomspace(0.1, 100, 600)

    fig, axs = plt.subplots(1, 2, figsize=(12.0, 4.6))

    # ---- Panel (a): real (Kerr) and imaginary (TPA) coefficients --------
    axs[0].loglog(Dw, alpha_kerr(Dw), color="C0", lw=2.0,
                  label=r"$|\alpha_{\rm eff,Kerr}|/|\beta|^2$ (real)")
    axs[0].loglog(Dw, gamma_tpa(Dw), color="C3", lw=2.0,
                  label=r"$\gamma_{\rm eff,2PA}/|\beta|^2$ (loss)")

    # Slope guides (placed off the lines to avoid label/curve overlap)
    Dw_far = Dw[Dw > 3]
    axs[0].loglog(Dw_far, 1.0 / Dw_far, ":", color="C0", lw=0.9, alpha=0.5)
    axs[0].text(50, 0.025, r"$\propto 1/\Delta\omega$", color="C0", fontsize=10)
    axs[0].loglog(Dw_far, 0.5 / Dw_far ** 2, ":", color="C3", lw=0.9, alpha=0.5)
    axs[0].text(50, 1.4e-4, r"$\propto 1/\Delta\omega^2$", color="C3",
                fontsize=10)

    axs[0].axvline(2.0, color="0.4", ls="-.", lw=0.9)
    axs[0].text(2.1, 3.5, r"design pt. $\Delta\omega = 2\gamma_2$",
                color="0.3", fontsize=9)
    axs[0].set_xlabel(r"$\Delta\omega / \gamma_2$")
    axs[0].set_ylabel(r"cascaded coefficient / $|\beta|^2$")
    axs[0].set_title("(a) Kerr falls as $1/\\Delta\\omega$; TPA falls faster, as $1/\\Delta\\omega^2$")
    axs[0].grid(True, which="both", alpha=0.3)
    axs[0].legend(loc="lower left", fontsize=9.5)
    axs[0].set_ylim(1e-4, 5)

    # ---- Panel (b): contrast ratio Kerr / TPA ----------------------------
    contrast_casc = alpha_kerr(Dw) / gamma_tpa(Dw)  # = 2 Dw
    axs[1].loglog(Dw, contrast_casc, color="C0", lw=2.2,
                  label=r"cascade: $|\alpha_{\rm Kerr}|/\gamma_{\rm 2PA} = 2\Delta\omega/\gamma_2$")
    axs[1].axhline(NATIVE_CONTRAST_FLOOR, color="C2", lw=2.0, ls="--",
                   label=r"native LiNbO$_3$ at telecom $\sim 10^4$")
    axs[1].axhline(1.0, color="0.6", lw=0.8, ls=":")
    axs[1].text(0.12, 1.3, r"$|\alpha_{\rm Kerr}| = \gamma_{\rm 2PA}$",
                color="0.4", fontsize=9)

    axs[1].axvline(2.0, color="0.4", ls="-.", lw=0.9)
    axs[1].axvline(0.5, color="0.4", ls="-.", lw=0.9, alpha=0.5)
    axs[1].text(0.52, 12, "BW-optimal\n($\\Delta\\omega = \\gamma_2/2$)",
                color="0.3", fontsize=8.5)
    axs[1].text(2.1, 12, "design pt.\n($\\Delta\\omega = 2\\gamma_2$)",
                color="0.3", fontsize=8.5)

    # Mark Kerr/loss values at the two key operating points
    for Dw_mark, label in [(0.5, "1"), (2.0, "4")]:
        c = 2 * Dw_mark
        axs[1].plot([Dw_mark], [c], "o", color="C0", ms=8, mec="k", mew=0.8)
        axs[1].annotate(f"$={label}$", xy=(Dw_mark, c),
                        xytext=(Dw_mark * 1.4, c * 0.6),
                        color="C0", fontsize=9)

    axs[1].set_xlabel(r"$\Delta\omega / \gamma_2$")
    axs[1].set_ylabel(r"switching contrast  $|\alpha_{\rm Kerr}|/\gamma_{\rm 2PA}$")
    axs[1].set_title("(b) Cascade is contrast-limited; native is essentially lossless")
    axs[1].grid(True, which="both", alpha=0.3)
    axs[1].legend(loc="lower right", fontsize=9.5)
    axs[1].set_ylim(0.5, 5e4)

    fig.suptitle("Loss budget: cascade carries unavoidable TPA-like absorption",
                 fontsize=11.5, y=1.00)
    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "loss_budget"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    print(f"  Kerr/TPA contrast at Domega = 0.5 (BW-optimal): {2*0.5:.1f}")
    print(f"  Kerr/TPA contrast at Domega = 2.0 (design pt):  {2*2.0:.1f}")
    print(f"  Kerr/TPA contrast at Domega = 10 (far-detuned): {2*10:.1f}")
    print(f"  Native contrast floor (LiNbO_3 telecom):        {NATIVE_CONTRAST_FLOOR:.0e}")
    plt.close(fig)


if __name__ == "__main__":
    main()
