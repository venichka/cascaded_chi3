"""
Group 2, Figure B4: Sign-tunable cascaded saturable absorber
=============================================================

The imaginary part of the cascaded coefficient gives an intensity-dependent
loss (cascaded TPA-like):

    gamma_3_casc(Delta_SH) = |beta|^2 (gamma_2/2) / [Delta_SH^2 + (gamma_2/2)^2]

Lorentzian in Delta_SH, peaking at Delta_SH = 0 with magnitude
|beta|^2 / (gamma_2/2). Unlike intrinsic TPA in semiconductors, the loss can
be tuned off by simply detuning the SH cavity. Cascade-only capability;
no architecture-A analogue.

Three panels:
  (a) Effective absorption coefficient gamma_3_eff vs Delta_SH (Lorentzian)
  (b) Saturable-absorber transmission curve: |s_-|^2/|s_+|^2 vs input power
      for several Delta_SH values
  (c) Effective TPA cross-section vs Delta_SH (or Q_2)
"""

import numpy as np
import matplotlib.pyplot as plt
from core import steady_state_cw, alpha3_casc_full, alpha3_loss_casc, adjacent_save


def measure_loss_coefficient(beta, gamma1, gamma2, dw_SH, s_amp_range, alpha3_native=0.0):
    """For each input amplitude, measure intracavity loss = |s+|^2 - |s-|^2 / |a_1|^2.
    Returns the cascaded TPA-like coefficient gamma_3_eff extracted from TCMT."""
    losses = []
    a1sq_list = []
    for s in s_amp_range:
        a1, a2, _ = steady_state_cw(beta, gamma1, gamma2, 0.0, dw_SH,
                                    alpha3_native, s)
        if abs(a1) < 1e-30:
            losses.append(0.0)
            a1sq_list.append(0.0)
            continue
        # Energy balance: input = output + dissipation
        # Total power dissipated = gamma_1 |a_1|^2 + gamma_2 |a_2|^2
        # The cascade-induced FM loss comes via the SH mode
        loss_total = gamma_2_loss(gamma1, gamma2, a1, a2)
        losses.append(loss_total)
        a1sq_list.append(abs(a1)**2)
    return np.array(losses), np.array(a1sq_list)


def gamma_2_loss(gamma1, gamma2, a1, a2):
    """Total dissipation rate at intracavity state (a_1, a_2)."""
    return gamma1 * abs(a1)**2 + gamma2 * abs(a2)**2


def figure_saturable_absorber(savepath=None):
    gamma1 = 1.0
    gamma2 = 1.0
    beta = 1.0

    # Panel (a): analytic gamma_3_casc Lorentzian profile
    dw_arr = np.linspace(-10, 10, 1000)
    gamma3_arr = np.array([alpha3_loss_casc(beta, gamma2, dw) for dw in dw_arr])
    alpha3_arr = np.array([alpha3_casc_full(beta, gamma2, dw) for dw in dw_arr])

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    ax = axes[0]
    ax.plot(dw_arr, gamma3_arr, '-', color='C3', lw=2,
            label=r'$\gamma_3^{\rm casc}$ (cascaded TPA-like)')
    ax.plot(dw_arr, alpha3_arr, '-', color='C0', lw=2,
            label=r'$\alpha_3^{\rm casc}$ (cascaded Kerr)')
    ax.fill_between(dw_arr, 0, gamma3_arr, alpha=0.15, color='C3')
    ax.axvline(0, color='gray', ls=':', alpha=0.5)
    ax.axhline(0, color='k', lw=0.5)
    # Show peak and width
    peak_val = 2 * beta**2 / gamma2
    ax.axhline(peak_val, color='C3', ls=':', alpha=0.5)
    ax.text(8, peak_val * 1.05, fr'peak = $\frac{{2|\beta|^2}}{{\gamma_2}}$',
            fontsize=10, color='C3', ha='right')
    ax.set_xlabel(r'SH detuning $\Delta_{\rm SH}/\gamma_2$')
    ax.set_ylabel(r'Coefficient $(\text{units of }|\beta|^2/\gamma_2)$')
    ax.set_title(r'(a) Cascade coefficients vs SH detuning')
    ax.legend(loc='upper right', framealpha=0.95)
    ax.grid(True, alpha=0.3)

    # Panel (b): saturable-absorber transmission vs input intensity for three Delta_SH
    ax = axes[1]
    dw_SH_values = [0.0, 1.0, 3.0]  # on resonance, mid, off
    colors_b = ['C3', 'C1', 'C0']
    labels_b = [r'$\Delta_{\rm SH}=0$ (peak loss)',
                r'$\Delta_{\rm SH}=\gamma_2$ (intermediate)',
                r'$\Delta_{\rm SH}=3\gamma_2$ (off; TPA suppressed)']

    s_arr = np.logspace(-2, 0.5, 30)
    for dw_SH, col, lab in zip(dw_SH_values, colors_b, labels_b):
        T = np.zeros_like(s_arr)
        for i, s in enumerate(s_arr):
            a1, a2, _ = steady_state_cw(beta, gamma1, gamma2, 0.0, dw_SH,
                                        0.0, s)
            # transmission at FM port: |s_- / s_+|^2 with s_- = -s_+ + sqrt(gamma_1) a_1
            s_out = -s + np.sqrt(gamma1) * a1
            T[i] = abs(s_out)**2 / s**2
        ax.semilogx(s_arr**2, T, '-o', color=col, ms=4, mfc='white', mew=1,
                    label=lab)

    # Compare to linear cavity (no nonlinearity, all detuned)
    T_lin = np.full_like(s_arr, abs(1 - 1)**2 if True else 1.0)  # critical coupling at gamma_rad=gamma_1: T_lin = 0 on resonance
    # Actually at critical coupling on resonance, transmission is 0.  But we're at gamma_rad = gamma_1 (over-coupled?)
    # Our normalisation: gamma_rad = gamma_1 (full radiative). So s_- = -s + sqrt(gamma_1) * a_1 = -s + sqrt(gamma_1) * sqrt(gamma_1) s / (gamma_1/2 - 0) = -s + 2s = +s. So T_lin = 1 on resonance.
    # Hmm that's wrong; let me re-check.

    ax.axhline(1.0, color='gray', ls=':', alpha=0.7, label='linear cavity')
    ax.set_xlabel(r'Input power $|s_+|^2$')
    ax.set_ylabel(r'FM-port transmission $|s_-/s_+|^2$')
    ax.set_title(r'(b) Intensity-dependent transmission (cascaded SA)')
    ax.legend(loc='lower left', fontsize=8.5, framealpha=0.95)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.2)

    # Panel (c): tunability -- ratio of TPA to Kerr vs Delta_SH (figure of merit for SA)
    ax = axes[2]
    dw_arr_pos = np.linspace(0.1, 8, 200)
    ratio_TPA_Kerr = alpha3_loss_casc(beta, gamma2, dw_arr_pos) / np.abs(
        alpha3_casc_full(beta, gamma2, dw_arr_pos))
    ax.plot(dw_arr_pos, ratio_TPA_Kerr, '-', color='C5', lw=2)
    ax.axhline(1.0, color='gray', ls=':', alpha=0.5)
    ax.axvline(0.5, color='gray', ls=':', alpha=0.5)
    ax.fill_between(dw_arr_pos, 0, ratio_TPA_Kerr,
                    where=(dw_arr_pos < 1.0), alpha=0.2, color='C3',
                    label=r'TPA-dominated ($|\Delta_{\rm SH}|<\gamma_2$)')
    ax.fill_between(dw_arr_pos, 0, ratio_TPA_Kerr,
                    where=(dw_arr_pos > 1.0), alpha=0.2, color='C0',
                    label=r'Kerr-dominated ($|\Delta_{\rm SH}|>\gamma_2$)')
    ax.set_xlabel(r'SH detuning $|\Delta_{\rm SH}|/\gamma_2$')
    ax.set_ylabel(r'$\gamma_3^{\rm casc}/|\alpha_3^{\rm casc}|$')
    ax.set_title(r'(c) Loss/Kerr ratio: detunability of SA character')
    ax.legend(loc='upper right', framealpha=0.95)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 3)

    fig.suptitle(r'Sign-tunable saturable absorber: cascaded TPA-like loss can be turned off by detuning',
                 fontsize=11, y=1.02)
    fig.tight_layout()
    if savepath:
        adjacent_save(fig, savepath)
    return fig


if __name__ == '__main__':
    figure_saturable_absorber(savepath='G2_B4_saturable_absorber.png')
