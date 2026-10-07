"""
Group 2, Figure B5: All-optical AM modulator via Delta_omega(t)
================================================================

Modulating the SH-cavity detuning (electro-optically, thermo-optically, or
piezo-electrically) at frequency Omega_mod periodically modulates the
cascaded Kerr coefficient alpha_3_casc(t), which in turn modulates the FM
transmission at fixed pump intensity.

Linear-response regime: delta_Omega << gamma_2, |Delta_SH|.

Three panels:
  (a) Time traces (input pump constant, Delta_SH modulated, T modulated)
  (b) Modulation amplitude vs Omega_mod  -> single-pole roll-off at gamma_2/2
  (c) Modulation depth vs working detuning Delta_SH_0 -> 1/Delta^2 scaling
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from core import adjacent_save


def tcmt_rhs_timedep_dwSH(t, y, beta, gamma1, gamma2, dw_FM, dw_SH_func,
                           alpha3_native, s_in_func):
    a1 = y[0] + 1j * y[1]
    a2 = y[2] + 1j * y[3]
    s_in = s_in_func(t)
    dw_SH = dw_SH_func(t)
    da1 = ((1j * dw_FM - gamma1 / 2.0) * a1
           + 1j * np.conj(beta) * np.conj(a1) * a2
           + 1j * alpha3_native * abs(a1) ** 2 * a1
           + np.sqrt(gamma1) * s_in)
    da2 = ((1j * dw_SH - gamma2 / 2.0) * a2 + 1j * beta * a1 ** 2)
    return [da1.real, da1.imag, da2.real, da2.imag]


def simulate_AM(beta, gamma1, gamma2, dw_FM, dw_SH_0, delta_Omega,
                Omega_mod, s_pump, T_total, n_pts=4000):
    dw_SH_func = lambda t: dw_SH_0 + delta_Omega * np.sin(Omega_mod * t)
    s_func = lambda t: s_pump + 0j
    t_eval = np.linspace(0, T_total, n_pts)
    sol = solve_ivp(
        tcmt_rhs_timedep_dwSH, [0, T_total],
        [0., 0., 0., 0.],
        args=(beta, gamma1, gamma2, dw_FM, dw_SH_func, 0.0, s_func),
        rtol=1e-10, atol=1e-13, method='DOP853',
        t_eval=t_eval, max_step=0.2 / max(Omega_mod, gamma2, 1.0),
    )
    a1 = sol.y[0] + 1j * sol.y[1]
    s_out = -s_pump + np.sqrt(gamma1) * a1
    T = np.abs(s_out) ** 2 / s_pump ** 2
    return sol.t, T


def extract_modulation_amplitude(t, T, Omega_mod, settle_fraction=0.5):
    i0 = int(len(t) * settle_fraction)
    t_s = t[i0:]
    T_s = T[i0:] - np.mean(T[i0:])
    cos_p = np.cos(Omega_mod * t_s)
    sin_p = np.sin(Omega_mod * t_s)
    a_c = 2 * np.mean(T_s * cos_p)
    a_s = 2 * np.mean(T_s * sin_p)
    return np.sqrt(a_c ** 2 + a_s ** 2)


def figure_AM_modulator(savepath=None):
    gamma1 = 1.0
    gamma2 = 1.0
    dw_FM = 0.0
    dw_SH_0 = 4.0
    delta_Omega = 0.05    # small modulation, linear response
    beta = 1.0
    s_pump = 0.3

    # Panel (a): time traces
    Omega_mod_demo = 0.5 * gamma2
    T_total = 25 * 2 * np.pi / Omega_mod_demo
    t_arr, T_demo = simulate_AM(beta, gamma1, gamma2, dw_FM, dw_SH_0,
                                  delta_Omega, Omega_mod_demo, s_pump, T_total)

    # Panel (b): transfer function vs Omega_mod
    Omega_mod_arr = np.logspace(-1.5, 1.5, 22) * gamma2
    mod_amplitudes = []
    for Om in Omega_mod_arr:
        T_total = max(30 * 2 * np.pi / Om, 60.0 / gamma2)
        t_a, T_a = simulate_AM(beta, gamma1, gamma2, dw_FM, dw_SH_0,
                                delta_Omega, Om, s_pump, T_total)
        amp = extract_modulation_amplitude(t_a, T_a, Om)
        mod_amplitudes.append(amp)
    mod_amplitudes = np.array(mod_amplitudes)

    # Panel (c): modulation amplitude vs working detuning at fixed (low) Omega_mod
    dw_SH_arr = np.linspace(2, 10, 9)
    Omega_fixed = 0.2 * gamma2
    mod_amp_vs_dw = []
    for dw0 in dw_SH_arr:
        T_total = 60.0 / gamma2
        t_a, T_a = simulate_AM(beta, gamma1, gamma2, dw_FM, dw0,
                                delta_Omega, Omega_fixed, s_pump, T_total)
        amp = extract_modulation_amplitude(t_a, T_a, Omega_fixed)
        mod_amp_vs_dw.append(amp)
    mod_amp_vs_dw = np.array(mod_amp_vs_dw)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    ax = axes[0]
    dw_SH_t = dw_SH_0 + delta_Omega * np.sin(Omega_mod_demo * t_arr)
    i0 = len(t_arr) // 2
    ax2 = ax.twinx()
    ax.plot(t_arr[i0:] * gamma2, T_demo[i0:], '-', color='C2', lw=1.5,
            label=r'FM transmission $|s_-/s_+|^2$')
    ax2.plot(t_arr[i0:] * gamma2, dw_SH_t[i0:], '--', color='C3', lw=1.2, alpha=0.7,
             label=r'$\Delta_{\rm SH}(t)$')
    ax.set_xlabel(r'$t \cdot \gamma_2$')
    ax.set_ylabel(r'Transmission', color='C2')
    ax2.set_ylabel(r'$\Delta_{\rm SH}(t)/\gamma_2$', color='C3')
    ax.tick_params(axis='y', labelcolor='C2')
    ax2.tick_params(axis='y', labelcolor='C3')
    ax.set_title(fr'(a) $\Omega_{{\rm mod}}=0.5\,\gamma_2$, $\delta\Omega = {delta_Omega}$ (linear response)')
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.loglog(Omega_mod_arr / gamma2, mod_amplitudes, 'o-', color='C0',
              ms=6, mfc='white', mew=1.5, label='TCMT measurement')
    Omega_fine = np.logspace(-1.5, 1.5, 200) * gamma2
    rolloff = 1.0 / np.sqrt((gamma2/2)**2 + Omega_fine**2)
    rolloff_norm = mod_amplitudes[0] * rolloff / (1.0 / (gamma2/2))
    ax.loglog(Omega_fine / gamma2, rolloff_norm, '--', color='k', lw=1.2,
              label=r'Single-pole low-pass: $\frac{\gamma_2/2}{\sqrt{(\gamma_2/2)^2+\Omega^2}}$')
    ax.axvline(0.5, color='gray', ls=':', alpha=0.7, label=r'$\Omega = \gamma_2/2$')
    ax.axvline(dw_SH_0, color='C3', ls=':', alpha=0.7,
               label=fr'$\Omega = \Delta_{{\rm SH,0}} = {dw_SH_0:.0f}\gamma_2$ (AC peak)')
    ax.set_xlabel(r'Modulation frequency $\Omega_{\rm mod}/\gamma_2$')
    ax.set_ylabel(r'Transmission modulation amplitude')
    ax.set_title(r'(b) Modulator transfer function: low-pass + AC peak at $\Omega=\Delta_{\rm SH}$')
    ax.legend(loc='lower left', fontsize=9, framealpha=0.95)
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[2]
    ax.loglog(dw_SH_arr / gamma2, mod_amp_vs_dw, 'o-', color='C3', ms=6,
              mfc='white', mew=1.5, label='TCMT measurement')
    dw_fine = np.logspace(np.log10(2), np.log10(10), 100)
    pred = 1.0 / (dw_fine * gamma2) ** 2
    norm = mod_amp_vs_dw[2] / pred[np.argmin(np.abs(dw_fine - dw_SH_arr[2]/gamma2))]
    ax.loglog(dw_fine, pred * norm, '--', color='k', lw=1.2,
              label=r'$\propto 1/\Delta_{\rm SH}^2$')
    ax.set_xlabel(r'Working detuning $\Delta_{\rm SH,0}/\gamma_2$')
    ax.set_ylabel(r'Transmission modulation amplitude')
    ax.set_title(fr'(c) Modulation depth vs detuning')
    ax.legend(loc='upper right', framealpha=0.95)
    ax.grid(True, which='both', alpha=0.3)

    fig.suptitle(r'All-optical AM modulator: modulating $\Delta\omega(t)$ '
                 r'modulates FM transmission  (cascade-only capability)',
                 fontsize=11, y=1.02)
    fig.tight_layout()
    if savepath:
        adjacent_save(fig, savepath)
    return fig


if __name__ == '__main__':
    figure_AM_modulator(savepath='G2_B5_AM_modulator.png')
