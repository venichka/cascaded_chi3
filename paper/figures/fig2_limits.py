"""Fig. 2 (draft): strength-bandwidth trade-off and host comparison.

(a) |chi_casc| at Domega = 2 gamma2 (left axis) and the bandwidth
    gamma2/2pi (right axis) versus Q2, for LiNbO3 with eta = 0.3.
(b) R = Q2 eta F / 6 (Eq. 6 of the main text) versus Q2 at
    Domega = 2 gamma2 and eta = 0.3 for four hosts; triangles mark R = 1.

Constants: Table S1 (paper/checks/constants/compute_output.txt).
"""

import numpy as np
import matplotlib.pyplot as plt

from paperstyle import apply, WIDTH, BLUE, RED, GREEN, ORANGE, panel_label, save

C0 = 2.99792458e8
LAM = 1.55e-6
W1 = 2 * np.pi * C0 / LAM
F2 = 2 * C0 / LAM                     # SH frequency (Hz)
ETA = 0.3

# LiNbO3: |chi2| = 37.5 pm/V, n2 = n_e(775 nm) = 2.178
CHI2_LN, N2_LN = 37.5e-12, 2.178
C_LN = 2 * W1 * CHI2_LN**2 * ETA / (3 * N2_LN**2)   # m^2 V^-2 s^-1

HOSTS = [  # name, F, line style, colour
    ("LiNbO$_3$", 0.1018, "-", BLUE),
    ("GaP", 4.535e-3, "--", RED),
    ("Al$_{0.2}$Ga$_{0.8}$As", 3.731e-3, "-.", GREEN),
    ("AlN", 2.468e-3, ":", ORANGE),
]


def main():
    apply()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(WIDTH, 2.15))

    # (a) strength and bandwidth versus Q2
    Q2 = np.logspace(2, 5, 200)
    g2 = 2 * W1 / Q2
    chi = C_LN / np.abs(2 * g2 + 0.5j * g2)
    bw = F2 / Q2 / 1e12                # THz
    Qm = np.logspace(2, 5, 7)
    ax1.loglog(Q2, chi, "-", color=BLUE)
    ax1.loglog(Qm, C_LN / np.abs((2 + 0.5j) * 2 * W1 / Qm), "o",
               mfc="white", mec=BLUE, mew=0.8)
    ax1.set_xlabel(r"$Q_2$")
    ax1.set_ylabel(r"$|\chi^{(3)}_{\mathrm{casc}}|$ (m$^2$/V$^2$)", color=BLUE)
    ax1.tick_params(axis="y", colors=BLUE, which="both")
    axb = ax1.twinx()
    axb.loglog(Q2, bw, "-", color=RED)
    axb.loglog(Qm, F2 / Qm / 1e12, "s", mfc="white", mec=RED, mew=0.8)
    axb.set_ylabel(r"$\gamma_2/2\pi$ (THz)", color=RED)
    axb.tick_params(axis="y", colors=RED, which="both")
    ax1.set_xlim(8e1, 1.25e5)
    panel_label(ax1, "(a)")

    # (b) R versus Q2 for the hosts
    Q2 = np.logspace(1, 5, 200)
    for name, F, ls, col in HOSTS:
        R = Q2 * ETA * F / 6
        ax2.loglog(Q2, R, ls, color=col)
        q1 = 6 / (ETA * F)
        ax2.plot([q1], [1.0], "v", color=col, ms=4)
        ax2.plot([], [], ls, color=col, label=name)
    ax2.axhline(1.0, color="0.75", lw=0.5, zorder=0)
    ax2.set_xlim(1e1, 1e5)
    ax2.set_ylim(1e-2, 1e3)
    ax2.set_xlabel(r"$Q_2$")
    ax2.set_ylabel(r"$R$")
    ax2.legend(loc="upper left", handlelength=2.2)
    panel_label(ax2, "(b)")

    fig.tight_layout(w_pad=0.8)
    save(fig, "fig2_limits")


if __name__ == "__main__":
    main()
