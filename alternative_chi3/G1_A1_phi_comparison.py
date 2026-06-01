"""
Group 1, Figure A1: SPM phase shift versus input power
======================================================

Compares architectures and shows the sign-tunability of the cascade:

  - A    (native-Kerr only):                beta = 0, alpha_3 > 0
  - B    (cascade only, Delta_SH < 0):      cascade self-focusing, same sign as A
  - A+B_constructive (Delta_SH < 0):        cascade reinforces native -> (1+R)
  - A+B_destructive  (Delta_SH > 0):        cascade opposes native -> (1-R)
"""

import numpy as np
import matplotlib.pyplot as plt
from core import steady_state_cw, alpha3_casc_full, OUT, adjacent_save


def measure_phi_NL(beta, gamma1, gamma2, dw_FM, dw_SH, alpha3_native, s_amp):
    a1, a2, _ = steady_state_cw(beta, gamma1, gamma2, dw_FM, dw_SH,
                                alpha3_native, s_amp)
    if abs(a1) < 1e-30:
        return 0.0, 0.0
    a1_linear = np.sqrt(gamma1) * s_amp / (gamma1 / 2.0 - 1j * dw_FM)
    phi = np.angle(a1 / a1_linear)
    return phi, abs(a1)**2


def figure_phi_NL_comparison(savepath=None):
    gamma1 = 1.0
    gamma2 = 1.0
    dw_FM = 0.0
    dw_SH_pos = 4.0
    dw_SH_neg = -4.0
    beta = 1.0

    alpha_3_casc_value = abs(alpha3_casc_full(beta, gamma2, dw_SH_pos))
    R_target = 10.0
    alpha3_native = alpha_3_casc_value / R_target

    s_arr = np.logspace(-3, 0.3, 35)
    phi = {'A': np.zeros_like(s_arr),
           'B_neg_dwSH': np.zeros_like(s_arr),
           'AB_neg_dwSH': np.zeros_like(s_arr),
           'AB_pos_dwSH': np.zeros_like(s_arr)}

    for i, s in enumerate(s_arr):
        phi['A'][i], _ = measure_phi_NL(
            beta=0.0, gamma1=gamma1, gamma2=gamma2, dw_FM=dw_FM,
            dw_SH=dw_SH_pos, alpha3_native=alpha3_native, s_amp=s)
        phi['B_neg_dwSH'][i], _ = measure_phi_NL(
            beta=beta, gamma1=gamma1, gamma2=gamma2, dw_FM=dw_FM,
            dw_SH=dw_SH_neg, alpha3_native=0.0, s_amp=s)
        phi['AB_neg_dwSH'][i], _ = measure_phi_NL(
            beta=beta, gamma1=gamma1, gamma2=gamma2, dw_FM=dw_FM,
            dw_SH=dw_SH_neg, alpha3_native=alpha3_native, s_amp=s)
        phi['AB_pos_dwSH'][i], _ = measure_phi_NL(
            beta=beta, gamma1=gamma1, gamma2=gamma2, dw_FM=dw_FM,
            dw_SH=dw_SH_pos, alpha3_native=alpha3_native, s_amp=s)

    # Small-signal predictions
    phi_A_lin = 8 * alpha3_native / gamma1**2 * s_arr**2
    phi_B_neg_lin = 8 * alpha_3_casc_value / gamma1**2 * s_arr**2
    phi_AB_constructive_lin = phi_A_lin + phi_B_neg_lin
    phi_AB_destructive_lin = phi_A_lin - phi_B_neg_lin

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.5))

    ax = axes[0]
    ax.plot(s_arr**2, phi['A'], 'o', color='C0', ms=5, mfc='white', mew=1.2,
            label=r'A (native only)')
    ax.plot(s_arr**2, phi['B_neg_dwSH'], 's', color='C3', ms=5, mfc='white', mew=1.2,
            label=r'B ($\Delta\omega<0$: cascade)')
    ax.plot(s_arr**2, phi['AB_neg_dwSH'], '^', color='C2', ms=5, mfc='white', mew=1.2,
            label=r'A+B, $\Delta\omega<0$  (constructive)')
    ax.plot(s_arr**2, phi['AB_pos_dwSH'], 'v', color='C1', ms=5, mfc='white', mew=1.2,
            label=r'A+B, $\Delta\omega>0$  (destructive)')
    ax.plot(s_arr**2, phi_A_lin, '-', color='C0', alpha=0.4, lw=1)
    ax.plot(s_arr**2, phi_B_neg_lin, '-', color='C3', alpha=0.4, lw=1)
    ax.plot(s_arr**2, phi_AB_constructive_lin, '-', color='C2', alpha=0.4, lw=1)
    ax.plot(s_arr**2, phi_AB_destructive_lin, '-', color='C1', alpha=0.4, lw=1)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xlabel(r'Input power $|s_+|^2$')
    ax.set_ylabel(r'Nonlinear phase $\phi_{\rm NL}$ (rad)')
    ax.set_title(r'(a) Sign-tunable cascade reinforces or cancels native Kerr')
    ax.legend(loc='upper left', fontsize=8.5, framealpha=0.95)
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    ax.set_xlim(1e-6, max(s_arr**2)*1.2)

    ax = axes[1]
    valid = np.abs(phi['A']) > 1e-12
    ratio_constructive = phi['AB_neg_dwSH'][valid] / phi['A'][valid]
    ratio_destructive = phi['AB_pos_dwSH'][valid] / phi['A'][valid]
    ax.semilogx(s_arr[valid]**2, ratio_constructive, '^-', color='C2',
                ms=5, mfc='white', mew=1.2,
                label=r'$\phi^{\rm A+B}/\phi^{\rm A}$ ($\Delta\omega<0$)')
    ax.semilogx(s_arr[valid]**2, ratio_destructive, 'v-', color='C1',
                ms=5, mfc='white', mew=1.2,
                label=r'$\phi^{\rm A+B}/\phi^{\rm A}$ ($\Delta\omega>0$)')
    ax.axhline(1 + R_target, color='C2', ls='--', lw=1.2, alpha=0.7,
               label=fr'$1+R={1+R_target:.0f}$ (small-signal)')
    ax.axhline(1 - R_target, color='C1', ls='--', lw=1.2, alpha=0.7,
               label=fr'$1-R={1-R_target:.0f}$ (small-signal)')
    ax.axhline(0, color='gray', lw=0.5, ls=':')
    ax.set_xlabel(r'Input power $|s_+|^2$')
    ax.set_ylabel(r'$\phi^{\rm A+B}_{\rm NL}\,/\,\phi^{\rm A}_{\rm NL}$')
    ax.set_title(r'(b) Cascade leverage: small-signal limit $1\pm R$, saturation visible')
    ax.legend(loc='center left', framealpha=0.95)
    ax.grid(True, alpha=0.3)

    fig.suptitle(r'Architecture A vs A+B: cascade reinforces ($\Delta\omega<0$) '
                 r'or cancels ($\Delta\omega>0$) the native $\chi^{(3)}$  '
                 fr'(target $|R|={R_target}$)',
                 fontsize=11, y=1.02)
    fig.tight_layout()
    if savepath: adjacent_save(fig, savepath)
    return fig


if __name__ == '__main__':
    figure_phi_NL_comparison(savepath='G1_A1_phi_comparison.png')
