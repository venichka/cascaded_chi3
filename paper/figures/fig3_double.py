"""Fig. 3 (draft): doubly resonant LiNbO3 metasurface.

(a) Placeholder for the unit cell and mode profiles of the chosen design.
(b) R = Q2 eta F / 6 for LiNbO3 (F = 0.102) versus Q2 and eta at
    Domega = 2 gamma2; the star marks Q2 = 4000, eta = 0.3.
(c) Stored pump energy versus input power for increasing and decreasing
    sweeps of the two-mode equations (Eqs. 1-2 of the main text), with the
    second-harmonic mode (solid) and without it (dashed).  Units beta =
    gamma2 = 1, N1 = N2 = 1/2; gamma1 = 1, gamma1^e = 0.5, pump-frame
    Domega = -2, alpha3 chosen so that R = 20, and in each case a pump
    detuning 1.3 times the minimum for bistability.  gamma1 = gamma2 keeps
    gamma1 below |Domega|, so the second harmonic stores less energy than the
    pump at the folds and the native terms involving it stay small.
(d) Ratio of four-wave-mixing efficiencies with and without the
    second-harmonic mode in the effective-Kerr description,
    |alpha3 + alpha_casc|^2 / alpha3^2 (solid), and (1+R)^2 (dashed),
    versus Domega/gamma2 < 0.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from paperstyle import (apply, WIDTH, BLUE, RED, GREY, panel_label, placeholder,
                        save)
from params import TCMTParams
from tcmt import integrate

F_LN = 0.1018
DW = -2.0                      # pump-frame Domega / gamma2
G1, G1E = 1.0, 0.5
ALPHA_C = -1.0 / (DW + 0.5j)   # cascaded coefficient, beta = gamma2 = 1
R_DESIGN = 20.0
ALPHA3 = abs(ALPHA_C.real) / R_DESIGN


def delta1_min(alpha):
    """Smallest |delta1| for bistability (Supplement 1, Section 4.C)."""
    rho = alpha.imag / abs(alpha.real)
    return (G1 / 2) / np.tan(np.radians(30.0) - np.arctan(rho))


def folds(alpha, d1):
    """Input powers at the two folds of the steady-state cubic."""
    a1, a2 = alpha.real, alpha.imag
    g = G1 / 2
    roots = np.roots([3 * abs(alpha) ** 2, 4 * (a1 * d1 + a2 * g),
                      d1 ** 2 + g ** 2])
    n = np.sort(roots[np.isreal(roots)].real)
    P = n * ((d1 + a1 * n) ** 2 + (g + a2 * n) ** 2) / G1E
    return n, P


def sweep(p, s2_grid, t_step=60.0):
    """Quasi-static sweep: each step starts from the previous state."""
    out, y = [], (0j, 0j)
    for s2 in s2_grid:
        sol = integrate(p, np.sqrt(s2), (0.0, t_step), y0=y)
        y = (sol.a[0, -1], sol.a[1, -1])
        out.append(abs(y[0]) ** 2)
    return np.array(out)


def main():
    apply()
    fig, axs = plt.subplots(2, 2, figsize=(WIDTH, 4.0))
    (ax_a, ax_b), (ax_c, ax_d) = axs

    placeholder(ax_a, "unit cell and mode profiles\nof the chosen design")
    panel_label(ax_a, "(a)")

    # (b) R(Q2, eta)
    Q2 = np.logspace(2, 5, 300)
    eta = np.logspace(-2, 0, 300)
    QQ, EE = np.meshgrid(Q2, eta)
    RR = QQ * EE * F_LN / 6
    pc = ax_b.pcolormesh(QQ, EE, RR, norm=LogNorm(1e-2, 2e3), cmap="Blues",
                         shading="auto", rasterized=True)
    cs = ax_b.contour(QQ, EE, RR, levels=[1, 10, 100, 1000], colors="k",
                      linewidths=0.5)
    for R_lab, q_lab in ((1, 1.6e2), (10, 1.3e3), (100, 1.3e4)):
        e_lab = 6 * R_lab / (q_lab * F_LN)
        ax_b.text(q_lab, e_lab, f"{R_lab:g}", fontsize=6.5, ha="center",
                  va="center", bbox=dict(fc="white", ec="none", pad=0.6))
    ax_b.plot([4000], [0.3], "*", color=RED, ms=7, mec="k", mew=0.4)
    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_xlabel(r"$Q_2$")
    ax_b.set_ylabel(r"$\eta_{\mathrm{geom}}$")
    cb = fig.colorbar(pc, ax=ax_b, pad=0.02)
    cb.set_label(r"$R$")
    cb.ax.tick_params(labelsize=7)
    panel_label(ax_b, "(b)")

    # (c) hysteresis with and without the SH mode
    cases = []
    for beta, alpha, ls, col in ((1.0, ALPHA3 + ALPHA_C, "-", BLUE),
                                 (0.0, ALPHA3 + 0j, "--", GREY)):
        d1 = -1.3 * delta1_min(alpha)          # red-shifting response
        p = TCMTParams(Domega=DW - 2 * d1, gamma1=G1, gamma1_rad=G1E,
                       alpha3=ALPHA3, delta1=d1, beta=beta)
        cases.append((p, alpha, ls, col, d1))
    for p, alpha, ls, col, d1 in cases:
        _, P = folds(alpha, d1)
        P_lo, P_hi = P.min(), P.max()
        s2 = np.logspace(np.log10(0.3 * P_lo), np.log10(1.12 * P_hi), 90)
        up = sweep(p, s2)
        down = sweep(p, s2[::-1])[::-1]
        ax_c.loglog(s2, up, ls, color=col)
        ax_c.loglog(s2, down, ls, color=col)
        n_f, _ = folds(alpha, d1)
        shift = abs(ALPHA_C) * n_f if p.beta else 0 * n_f
        print(f"beta={p.beta}: delta1 = {d1:.2f}, cold mismatch = {DW - 2*d1:.2f}, "
              f"folds at P = {P_lo:.3g}, {P_hi:.3g}; U2/U1 at folds = "
              f"{shift / abs(DW + 0.5j)}")
    ax_c.set_xticks([3, 10, 30])
    ax_c.set_xticklabels(["3", "10", "30"])
    ax_c.set_xlabel(r"input power $|s|^2$ ($\gamma_2^3/|\beta|^2$)")
    ax_c.set_ylabel(r"stored energy $|a_1|^2$ ($\gamma_2^2/|\beta|^2$)")
    panel_label(ax_c, "(c)")

    # (d) FWM efficiency ratio
    dw = np.linspace(-8, -0.5, 300)
    ac = -1.0 / (dw + 0.5j)
    ratio = np.abs(ALPHA3 + ac) ** 2 / ALPHA3 ** 2
    R = np.abs(ac.real) / ALPHA3
    ax_d.semilogy(dw, ratio, "-", color=BLUE)
    ax_d.semilogy(dw, (1 + R) ** 2, "--", color=GREY)
    ax_d.set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax_d.set_ylabel(r"$\eta_{\mathrm{FWM}}/\eta_{\mathrm{FWM}}^{\mathrm{nat}}$")
    ax_d.set_xlim(-8, -0.5)
    panel_label(ax_d, "(d)")

    fig.tight_layout(h_pad=0.6, w_pad=0.8)
    save(fig, "fig3_double")


if __name__ == "__main__":
    main()
