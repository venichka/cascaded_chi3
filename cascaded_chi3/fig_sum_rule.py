"""Bandwidth-enhancement sum rule (CW).

The cascaded Kerr response has a complex transfer function for pump
amplitude modulation at frequency Omega:

    kappa_casc(Omega) = -|beta|^2 / (gamma2/2 - i*(Omega + Delta_omega))

To make the sum rule visible directly, we fix the SH detuning
Delta_omega in physical units (so the lineshape of every curve is
centred at the same Omega) and vary gamma2 = omega2 / Q2. Then:

  - peak |kappa_casc|         = |beta|^2 / (gamma2/2)         ~  Q2
  - FWHM of |kappa_casc(Omega)|     = sqrt(3) * gamma2          ~  1/Q2
  - peak * FWHM               = sqrt(3) * |beta|^2 / (1/2) = 2*sqrt(3)
                                                              (constant)

This is exactly the sum rule of Section 5: resonant enhancement
redistributes nonlinear action; it does not increase the area.

All three Q2 curves are evaluated on the SAME Omega grid so the
peaks line up at Omega = -Delta_omega.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from params import OMEGA1_REAL

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def kappa(Omega, gamma2, Domega):
    return -1.0 / (0.5 * gamma2 - 1j * (Omega + Domega))


def fwhm(x, y):
    """FWHM of a single-peaked nonneg array y(x) sampled on uniform x;
    linear interp on each side."""
    y_max = y.max()
    half = 0.5 * y_max
    above = y >= half
    if not above.any():
        return np.nan
    i_lo = above.argmax()
    i_hi = len(y) - 1 - above[::-1].argmax()
    if i_lo == 0 or i_hi == len(y) - 1:
        return np.nan
    x_l = np.interp(half, [y[i_lo - 1], y[i_lo]], [x[i_lo - 1], x[i_lo]])
    x_r = np.interp(half, [y[i_hi + 1], y[i_hi]], [x[i_hi + 1], x[i_hi]])
    return x_r - x_l


def main():
    omega2 = 2 * OMEGA1_REAL
    Q2_values = [1_000, 4_000, 16_000]
    colors = ["C0", "C1", "C3"]

    # Fix Delta_omega in physical units. Choosing it = (gamma2 at the
    # mid Q2)/2 puts the mid-Q2 curve on the bandwidth-optimal point.
    gamma2_mid = omega2 / Q2_values[1]
    Domega = 0.5 * gamma2_mid

    # One shared Omega grid: wide enough for the broadest curve
    # (lowest Q2 -> largest gamma2), centred on Omega = -Delta_omega so
    # all peaks land at the same point on the axis.
    gamma2_max = omega2 / min(Q2_values)
    Omega = np.linspace(-Domega - 6 * gamma2_max,
                        -Domega + 6 * gamma2_max, 8000)

    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.6))

    peak_record = []
    fwhm_record = []
    Omega_GHz = Omega / (2 * np.pi) / 1e9
    Omega_GHz_offset = (Omega + Domega) / (2 * np.pi) / 1e9

    for Q2, c in zip(Q2_values, colors):
        gamma2 = omega2 / Q2
        K = np.abs(kappa(Omega, gamma2, Domega))
        peak_record.append(K.max())
        fwhm_record.append(fwhm(Omega, K))
        axs[0].plot(Omega_GHz_offset, K / peak_record[0], color=c,
                    label=f"$Q_2 = {Q2}$  ($\\gamma_2/2\\pi = "
                          f"{gamma2/(2*np.pi)/1e9:.1f}$ GHz)")

    axs[0].set_xlabel(r"$(\Omega - \Omega_{\rm peak})/(2\pi)$  [GHz]")
    axs[0].set_ylabel(r"$|\kappa_{\rm casc}(\Omega)|$  (normalised to lowest $Q_2$)")
    axs[0].set_title("(a) Frequency response — peaks aligned")
    axs[0].legend(fontsize=9)
    axs[0].grid(True, alpha=0.3)
    axs[0].set_xlim(Omega_GHz_offset[0], Omega_GHz_offset[-1])

    Q2_arr = np.array(Q2_values, float)
    peak_arr = np.array(peak_record)
    fwhm_arr = np.array(fwhm_record)
    product = peak_arr * fwhm_arr

    axs[1].loglog(Q2_arr, peak_arr / peak_arr[0], "o-", color="C0",
                  label=r"peak  $\propto Q_2$")
    axs[1].loglog(Q2_arr, fwhm_arr / fwhm_arr[0], "s-", color="C3",
                  label=r"FWHM  $\propto 1/Q_2$")
    axs[1].loglog(Q2_arr, product / product[0], "D-", color="k",
                  label=r"peak $\times$ FWHM  (constant)")
    axs[1].loglog(Q2_arr, Q2_arr / Q2_arr[0], ":", color="C0", alpha=0.5)
    axs[1].loglog(Q2_arr, Q2_arr[0] / Q2_arr, ":", color="C3", alpha=0.5)
    axs[1].set_xlabel(r"$Q_2$")
    axs[1].set_ylabel("normalised value")
    axs[1].set_title("(b) Bandwidth-enhancement sum rule")
    axs[1].legend(fontsize=9)
    axs[1].grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "sum_rule"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    print(f"  peak ratios   : {peak_arr / peak_arr[0]}")
    print(f"  FWHM ratios   : {fwhm_arr / fwhm_arr[0]}")
    print(f"  product ratios: {product / product[0]}")
    plt.close(fig)


if __name__ == "__main__":
    main()
