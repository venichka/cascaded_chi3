"""
Group 1, Figure A3: Q-scaling of phi_NL for native vs cascade architectures
===========================================================================

Holding everything else fixed (same FM cavity Q_1, same |chi^(2)|, same pump),
sweep the SH cavity Q_2. The native architecture is independent of Q_2
(phi_NL flat), while the cascade architecture has phi_NL linear in Q_2.

This is the cleanest single visualization of the architecture difference:
the cascade has a free extra knob ($Q_2$) that the native does not.

Two panels:
  (a) phi_NL vs Q_2 for both architectures (and combined A+B), log-log
  (b) Ratio (A+B)/A vs Q_2 -- grows as 1 + a*Q_2 in proportion to material R
"""

import numpy as np
import matplotlib.pyplot as plt
from core import steady_state_cw, alpha3_casc_full, adjacent_save


def measure_phi_NL(beta, gamma1, gamma2, dw_FM, dw_SH, alpha3_native, s_amp):
    a1, a2, _ = steady_state_cw(beta, gamma1, gamma2, dw_FM, dw_SH,
                                alpha3_native, s_amp)
    if abs(a1) < 1e-30:
        return 0.0
    a1_linear = np.sqrt(gamma1) * s_amp / (gamma1 / 2.0 - 1j * dw_FM)
    return np.angle(a1 / a1_linear)


def figure_Q_scaling(savepath=None):
    """Hold gamma_1 = 1, sweep gamma_2 (i.e. Q_2 = omega_2/gamma_2 inversely)."""
    gamma1 = 1.0
    dw_FM = 0.0
    s_pump = 0.05    # small signal -- linear regime
    # Hold beta * gamma_2 such that we represent fixed |chi^(2)| (which sets beta)
    # while varying Q_2 = omega_2/gamma_2 -> 1/gamma_2.
    # In our normalization: beta is a material+overlap parameter, independent of gamma_2
    # gamma_2 is what we sweep. dw_SH = const * gamma_2 (i.e. fixed in units of gamma_2_ref)
    # OR: hold dw_SH fixed in absolute units? Let's hold dw_SH/gamma_2 = const (typical operating).
    beta = 1.0
    gamma2_arr = np.logspace(-1.5, 0.5, 12)   # gamma_2 from 0.03 to 3 (Q_2 from 32 down to 0.3)
    Q2_arr = 1.0 / gamma2_arr                  # in arbitrary units

    # Fix native alpha_3 so that at gamma_2 = 1 (mid-range), R = 5
    gamma2_ref = 1.0
    dw_SH_over_gamma2 = 4.0   # operating point: 4 gamma_2 off resonance
    alpha_3_casc_ref = abs(alpha3_casc_full(beta, gamma2_ref, dw_SH_over_gamma2 * gamma2_ref))
    R_ref = 5.0
    alpha3_native = alpha_3_casc_ref / R_ref

    phi_A = np.zeros_like(gamma2_arr)
    phi_AB = np.zeros_like(gamma2_arr)
    phi_B = np.zeros_like(gamma2_arr)

    for i, g2 in enumerate(gamma2_arr):
        dw_SH = dw_SH_over_gamma2 * g2   # fixed in units of gamma_2 (typical adiabatic)

        # A only: alpha_3_native, no cascade
        phi_A[i] = measure_phi_NL(beta=0.0, gamma1=gamma1, gamma2=g2, dw_FM=dw_FM,
                                  dw_SH=dw_SH, alpha3_native=alpha3_native, s_amp=s_pump)
        # B only
        phi_B[i] = measure_phi_NL(beta=beta, gamma1=gamma1, gamma2=g2, dw_FM=dw_FM,
                                  dw_SH=-dw_SH,  # use negative for constructive sign (same as alpha3_native > 0)
                                  alpha3_native=0.0, s_amp=s_pump)
        # A+B (constructive)
        phi_AB[i] = measure_phi_NL(beta=beta, gamma1=gamma1, gamma2=g2, dw_FM=dw_FM,
                                   dw_SH=-dw_SH,
                                   alpha3_native=alpha3_native, s_amp=s_pump)

    # Compute the actual R at each gamma_2
    R_arr = np.array([abs(alpha3_casc_full(beta, g2, dw_SH_over_gamma2*g2)) / alpha3_native
                       for g2 in gamma2_arr])
    # alpha_3_casc ~ |beta|^2 / |Delta_omega|; with Delta_omega = 4 gamma_2, this scales as 1/gamma_2 = Q_2
    # so R should scale linearly with Q_2

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.5))

    # Panel (a): phi vs Q_2 for three architectures
    ax = axes[0]
    ax.loglog(Q2_arr, np.abs(phi_A), 'o-', color='C0', ms=6, mfc='white', mew=1.5,
              label='A (native only)')
    ax.loglog(Q2_arr, np.abs(phi_B), 's-', color='C3', ms=6, mfc='white', mew=1.5,
              label='B (cascade only)')
    ax.loglog(Q2_arr, np.abs(phi_AB), '^-', color='C2', ms=6, mfc='white', mew=1.5,
              label='A+B (combined)')
    # Reference: phi_native is Q_2-independent
    ax.axhline(np.abs(phi_A[0]), color='C0', ls=':', alpha=0.5,
               label=r'$\phi^A$ flat in $Q_2$')
    # Cascade: phi_B linear in Q_2
    Q2_fine = np.logspace(-0.7, 1.7, 100)
    pred_B = np.abs(phi_B[6]) * (Q2_fine / Q2_arr[6])   # linear in Q_2
    ax.loglog(Q2_fine, pred_B, ':', color='C3', alpha=0.5, lw=1.5,
              label=r'$\phi^B \propto Q_2$')
    ax.set_xlabel(r'SH cavity quality factor $Q_2$ (arb. units)')
    ax.set_ylabel(r'$|\phi_{\rm NL}|$ (rad)')
    ax.set_title(r'(a) $\phi_{\rm NL}$ vs $Q_2$:  native is flat, cascade is linear')
    ax.legend(loc='lower right', fontsize=9, framealpha=0.95)
    ax.grid(True, which='both', alpha=0.3)

    # Panel (b): ratio (A+B)/A vs Q_2 = 1 + R_actual
    ax = axes[1]
    ratio = np.abs(phi_AB) / np.abs(phi_A)
    ax.semilogx(Q2_arr, ratio, 'D', color='C4', ms=8, mfc='white', mew=2,
                label=r'$\phi^{\rm A+B}/\phi^{\rm A}$ (TCMT)')
    ax.semilogx(Q2_arr, 1 + R_arr, '-', color='k', lw=1.5,
                label=r'$1 + R$ (analytic; $R\propto Q_2$)')
    ax.axhline(1.0, color='gray', ls=':', alpha=0.5)
    ax.set_xlabel(r'SH cavity quality factor $Q_2$ (arb. units)')
    ax.set_ylabel(r'$\phi^{\rm A+B}_{\rm NL}\,/\,\phi^{\rm A}_{\rm NL}$')
    ax.set_title(r'(b) Cascade enhancement grows linearly with $Q_2$')
    ax.legend(loc='upper left', framealpha=0.95)
    ax.grid(True, alpha=0.3)

    fig.suptitle(r'Q-scaling: $\phi^A$ independent of $Q_2$, $\phi^B \propto Q_2$  '
                 fr'(fixed FM cavity $Q_1$, $\Delta\omega/\gamma_2={dw_SH_over_gamma2:.0f}$)',
                 fontsize=11, y=1.02)
    fig.tight_layout()
    if savepath:
        adjacent_save(fig, savepath)
    return fig


if __name__ == '__main__':
    figure_Q_scaling(savepath='G1_A3_Q_scaling.png')
