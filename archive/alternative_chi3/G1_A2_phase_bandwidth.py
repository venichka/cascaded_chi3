"""
Group 1, Figure A2: Phase-bandwidth locus
==========================================

Shows the three families of operating points on a log-log phi_NL vs Delta_f plot:

  - Native, sweep Q_1 at fixed material:        slope -2 (Q_1^2 / (1/Q_1))
  - Cascade, fix Q_1, sweep Q_2:                slope -1 (Q_2 / (1/Q_2)) [sum rule]
  - Cascade, fix Q_2, sweep Q_1:                slope -2 (Q_1^2 / (1/Q_1))

Plus crossover points for several materials.

This is the central design-decision figure of section 11.
"""

import numpy as np
import matplotlib.pyplot as plt
from core import MATERIALS, cascade_ratio_R, crossover_Q2, adjacent_save


def figure_phase_bandwidth_locus(savepath=None):
    omega1 = 1.22e15  # FM angular freq (LiNbO3 telecom, rad/s)
    eta_geom = 0.3

    # Sweep ranges
    Q1_arr = np.logspace(2, 5, 50)         # FM Q range
    Q2_arr = np.logspace(2, 5, 50)         # SH Q range
    Q1_fixed = 1000.0                       # for cascade-sweep-Q2 family
    Q2_fixed = 4000.0                       # for cascade-sweep-Q1 family

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel (a): LiNbO3 case in detail
    mat = MATERIALS['LiNbO3']
    ax = axes[0]

    # Family A: native only, sweep Q_1
    # phi_NL ~ Q_1^2 * chi3_native;  Delta_f = omega/(2 pi Q_1)
    Df_A = omega1 / (2 * np.pi * Q1_arr)
    phi_A = Q1_arr**2 * mat.chi3_native       # arbitrary units, but proportionality is correct
    ax.loglog(Df_A, phi_A, '-', color='C0', lw=2,
              label=r'A: native, sweep $Q_1$ (slope $-2$)')

    # Family B_II: cascade, fix Q_1, sweep Q_2
    # phi_NL ~ Q_1^2 * |chi^(2)|^2 eta / n_2^2 * Q_2;  Delta_f = omega/(2 pi Q_2)
    Df_B2 = omega1 / (2 * np.pi * Q2_arr)
    phi_B2 = Q1_fixed**2 * (mat.chi2**2 * eta_geom / mat.n2sq) * Q2_arr
    ax.loglog(Df_B2, phi_B2, '-', color='C2', lw=2,
              label=fr'B: cascade, $Q_1={Q1_fixed:.0f}$, sweep $Q_2$ (slope $-1$, sum rule)')

    # Family B_I: cascade, fix Q_2, sweep Q_1
    # phi_NL ~ Q_1^2 * |chi^(2)|^2 eta / n_2^2 * Q_2_fixed;  Delta_f = omega/(2 pi Q_1)
    Df_B1 = omega1 / (2 * np.pi * Q1_arr)
    phi_B1 = Q1_arr**2 * (mat.chi2**2 * eta_geom / mat.n2sq) * Q2_fixed
    ax.loglog(Df_B1, phi_B1, '-', color='C3', lw=2,
              label=fr'B: cascade, $Q_2={Q2_fixed:.0f}$, sweep $Q_1$ (slope $-2$)')

    # Annotate slope regions
    ax.set_xlabel(r'Usable bandwidth $\Delta f$ [Hz]')
    ax.set_ylabel(r'$\phi_{\rm NL}$ at fixed $|s_+|^2$  (arb. units)')
    ax.set_title(r'(a) Phase-bandwidth locus for LiNbO$_3$')
    ax.legend(loc='upper right', fontsize=9, framealpha=0.95)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(1e9, 1e13)

    # Panel (b): A vs A+B for several materials, plot R vs Q_2
    ax = axes[1]
    Q2_arr2 = np.logspace(2, 5, 100)
    for key, mat in MATERIALS.items():
        if mat.chi2 == 0:
            continue  # Si, Si3N4 have no cascade
        R_arr = np.array([cascade_ratio_R(mat, Q2) for Q2 in Q2_arr2])
        ax.loglog(Q2_arr2, R_arr, lw=2, label=mat.name)
        # Mark crossover
        Qc = crossover_Q2(mat)
        if 1e2 < Qc < 1e5:
            ax.axvline(Qc, color=ax.lines[-1].get_color(), ls=':', alpha=0.5)
            ax.annotate(fr'$Q_2^{{\rm cross}}\!\sim\!{Qc:.0e}$',
                        xy=(Qc, 1.0), xytext=(Qc * 1.3, 0.2),
                        color=ax.lines[-1].get_color(), fontsize=8.5)
    ax.axhline(1.0, color='k', ls='--', lw=1.2, alpha=0.7,
               label=r'$R=1$: cascade equals native')
    ax.set_xlabel(r'SH cavity $Q_2$')
    ax.set_ylabel(r'Cascade ratio $R=|\chi^{(3)}_{\rm casc}|/|\chi^{(3)}_{\rm native}|$')
    ax.set_title(r'(b) Material survey: $R$ vs $Q_2$ at $\Delta\omega/\gamma_2=2$, $\eta_{\rm geom}=0.3$')
    ax.legend(loc='upper left', fontsize=8.5, framealpha=0.95)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(1e2, 1e5)
    ax.set_ylim(1e-3, 1e3)

    fig.suptitle(r'Architecture decision: phase-bandwidth tradeoff and $R$ vs $Q_2$',
                 fontsize=11, y=1.01)
    fig.tight_layout()
    if savepath: adjacent_save(fig, savepath)
    return fig


if __name__ == '__main__':
    figure_phase_bandwidth_locus(savepath='G1_A2_phase_bandwidth.png')
