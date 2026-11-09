"""Fig. 4 (draft): singly resonant metasurface and the proposed Z-scan.

(a) Placeholder for the unit cell and second-harmonic mode of the chosen
    design.
(b) Real (solid) and imaginary (dashed) parts of chi_casc / chi_nat versus
    Domega / gamma2 for a design with R = 3 at |Domega| = 2 gamma2; the
    dotted line is the native susceptibility.
(c) Closed-aperture Z-scan transmission change of a thin sample,
    Delta T = [4 x DPhi - 2 (x^2+3) DPsi] / [(x^2+9)(x^2+1)], x = z/z0,
    with the refractive (DPhi) and absorptive (DPsi) phases set by
    chi_nat + chi_casc at Domega = +2 gamma2 (solid), -2 gamma2 (dashed) and
    by chi_nat alone (dotted), in units of the native phase DPhi_nat.
(d) Z-scan (I^2-weighted) nonlinear phase for Gaussian pulses relative to
    its continuous-wave value versus the pulse duration (intensity FWHM),
    at Domega = 2 gamma2 for Q2 = 1000 (solid) and 4000 (dashed), pump at
    1.55 um.  Model: cascaded_chi3/sample_pulsed.py.
"""

import numpy as np
import matplotlib.pyplot as plt

from paperstyle import (apply, WIDTH, BLUE, RED, GREY, panel_label, placeholder,
                        save)
from sample_pulsed import reduction_factors

R0 = 3.0                      # draft value of R at |Domega| = 2 gamma2
C0 = 2.99792458e8
W1 = 2 * np.pi * C0 / 1.55e-6


def chi_ratio(x):
    """chi_casc / chi_nat for R(2 gamma2) = R0, x = Domega/gamma2."""
    shape = -1.0 / (x + 0.5j)
    return shape * R0 / abs((-1.0 / (2 + 0.5j)).real)


def zscan(x, dphi, dpsi):
    return (4 * x * dphi - 2 * (x ** 2 + 3) * dpsi) / ((x ** 2 + 9) * (x ** 2 + 1))


def main():
    apply()
    fig, axs = plt.subplots(2, 2, figsize=(WIDTH, 4.0))
    (ax_a, ax_b), (ax_c, ax_d) = axs

    placeholder(ax_a, "unit cell and second-harmonic\nmode of the chosen design")
    panel_label(ax_a, "(a)")

    # (b)
    x = np.linspace(-6, 6, 601)
    c = chi_ratio(x)
    ax_b.axhline(0, color="0.75", lw=0.5, zorder=0)
    ax_b.plot(x, c.real, "-", color=BLUE)
    ax_b.plot(x, c.imag, "--", color=RED)
    ax_b.axhline(1.0, color=GREY, ls=":", lw=1.0)
    ax_b.set_xlim(-6, 6)
    ax_b.set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax_b.set_ylabel(r"$\chi^{(3)}_{\mathrm{casc}}/\chi^{(3)}_{\mathrm{nat}}$")
    panel_label(ax_b, "(b)")

    # (c)
    z = np.linspace(-6, 6, 601)
    for xd, ls, col in ((2.0, "-", BLUE), (-2.0, "--", RED)):
        cr = chi_ratio(np.array([xd]))[0]
        dphi = 1.0 + cr.real
        dpsi = cr.imag
        ax_c.plot(z, zscan(z, dphi, dpsi), ls, color=col)
    ax_c.plot(z, zscan(z, 1.0, 0.0), ":", color=GREY, lw=1.0)
    ax_c.axhline(0, color="0.75", lw=0.5, zorder=0)
    ax_c.set_xlim(-6, 6)
    ax_c.set_xlabel(r"$z/z_0$")
    ax_c.set_ylabel(r"$\Delta T/\Delta\Phi_{\mathrm{nat}}$")
    panel_label(ax_c, "(c)")

    # (d)
    tau = np.logspace(-1, 1.3, 40)               # ps
    for Q2, ls in ((1000, "-"), (4000, "--")):
        g2 = 2 * W1 / Q2
        r = [reduction_factors(2.0, g2 * t * 1e-12)[0] for t in tau]
        ax_d.semilogx(tau, r, ls, color=BLUE)
    ax_d.axhline(1.0, color="0.75", lw=0.5, zorder=0)
    ax_d.set_xlabel(r"pulse duration, intensity FWHM (ps)")
    ax_d.set_ylabel(r"$\Delta\Phi/\Delta\Phi_{\mathrm{CW}}$")
    ax_d.set_xlim(tau[0], tau[-1])
    panel_label(ax_d, "(d)")

    fig.tight_layout(h_pad=0.6, w_pad=0.8)
    save(fig, "fig4_single")


if __name__ == "__main__":
    main()
