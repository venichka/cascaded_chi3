"""
Cascaded chi^(2):chi^(2) -> chi^(3)_eff in a doubly-resonant metasurface
=======================================================================

TCMT simulation of the equations from the derivation note (Section 4):

    dot a_2 = (i Delta_omega - gamma_2/2) a_2 + i beta a_1^2
    dot a_1 = (i delta_p - gamma_1/2) a_1 + i beta* a_1* a_2 + i alpha_3_native |a_1|^2 a_1
              + sqrt(gamma_1_rad) s_+

with delta_p the FM pump-cavity detuning (in the rotating frame at pump freq),
Delta_omega = 2*omega_pump - omega_2 the SH detuning in the rotating frame.

UNITS: gamma_2 = 1 sets the time scale.  All frequencies are in units of gamma_2;
all amplitudes are dimensionless (the absolute scale is set by beta).

DEMONSTRATIONS (one figure each):
  Figure 1 - Validation + sign tunability: TCMT vs analytic alpha_3_casc(Delta_omega)
  Figure 2 - Sum rule: peak alpha_3_eff x bandwidth is Q_2-independent
  Figure 3 - Bistability: hysteresis in |a_1|^2 vs |s_+|^2 (CW, off-resonance pump)
  Figure 4 - Saturation: small-signal a_2/a_1^2 ratio breaks down at high power
  Figure 5 - Pulse-domain bandwidth limit: cascade follows only if pulse is long vs 1/gamma_2
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os

# ---------------------- Plot styling ----------------------
plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 11.5,
    'axes.labelsize': 11,
    'legend.fontsize': 9.5,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi': 110,
    'savefig.dpi': 160,
    'savefig.bbox': 'tight',
    'mathtext.fontset': 'cm',
})

# ============================================================
#                     CORE TCMT INTEGRATOR
# ============================================================

def tcmt_rhs(t, y, beta, gamma1, gamma2, dw_FM, dw_SH, alpha3_native, s_in_func):
    """
    Right-hand side of the TCMT equations in the rotating frame at pump freq.

    State vector y = [Re(a_1), Im(a_1), Re(a_2), Im(a_2)].

    Parameters
    ----------
    beta             : chi^(2) coupling constant (complex; here typically real positive)
    gamma1, gamma2   : intrinsic decay rates (FM, SH)
    dw_FM            : FM detuning omega_pump - omega_1 (positive = blue-pumped)
    dw_SH            : SH detuning 2*omega_pump - omega_2 (= Delta_omega when dw_FM=0)
    alpha3_native    : bare Kerr coefficient (real)
    s_in_func        : callable s_in_func(t) returning complex input amplitude
    """
    a1 = y[0] + 1j*y[1]
    a2 = y[2] + 1j*y[3]
    s_in = s_in_func(t)

    da1 = ((1j*dw_FM - gamma1/2.0) * a1
           + 1j * np.conj(beta) * np.conj(a1) * a2
           + 1j * alpha3_native * abs(a1)**2 * a1
           + np.sqrt(gamma1) * s_in)             # (assume gamma_1_rad = gamma_1)
    da2 = ((1j*dw_SH - gamma2/2.0) * a2
           + 1j * beta * a1**2)

    return [da1.real, da1.imag, da2.real, da2.imag]


def steady_state_cw(beta, gamma1, gamma2, dw_FM, dw_SH, alpha3_native, s_amp,
                    t_settle=None, y0=None):
    """Run TCMT under CW input until steady state. Returns (a1, a2, final_state)."""
    if t_settle is None:
        t_settle = max(60.0/gamma2, 60.0/gamma1)
    if y0 is None:
        y0 = [0.0, 0.0, 0.0, 0.0]
    s_func = lambda t: s_amp + 0j
    sol = solve_ivp(
        tcmt_rhs, [0, t_settle], y0,
        args=(beta, gamma1, gamma2, dw_FM, dw_SH, alpha3_native, s_func),
        rtol=1e-10, atol=1e-13, method='DOP853',
    )
    a1 = sol.y[0, -1] + 1j*sol.y[1, -1]
    a2 = sol.y[2, -1] + 1j*sol.y[3, -1]
    return a1, a2, sol.y[:, -1]


# ============================================================
#                  ANALYTIC PREDICTIONS
# ============================================================

def alpha3_casc_full(beta, gamma2, dw):
    """Real part of cascaded back-action coefficient (Kerr-like).
    Im(back-action / |a1|^2 a1) = - |beta|^2 * Delta_omega / (Delta_omega^2 + gamma_2^2/4)
    """
    return -abs(beta)**2 * dw / (dw**2 + (gamma2/2)**2)


def alpha3_loss_casc(beta, gamma2, dw):
    """Cascaded TPA-like loss coefficient (positive)."""
    return abs(beta)**2 * (gamma2/2.0) / (dw**2 + (gamma2/2)**2)


def alpha3_casc_adiabatic(beta, dw):
    """Adiabatic limit Kerr: -|beta|^2 / Delta_omega."""
    return -abs(beta)**2 / dw


# ============================================================
#                     OUTPUT DIRECTORY
# ============================================================
OUT = '/home/claude/sim/figures'
os.makedirs(OUT, exist_ok=True)


# ============================================================
#   FIGURE 1: Validation + sign tunability (Delta_omega sweep)
# ============================================================

def figure1_validation_sign_tunability(beta=1.0, gamma1=5.0, gamma2=1.0,
                                       s_amp=0.05, savepath=None):
    """
    Sweep SH detuning Delta_omega.  At each value, drive CW with small input
    and extract the effective Kerr coefficient and TPA-like loss from the
    cascaded back-action.  Compare with the closed-form expressions.

    Demonstrates:
      (a) TCMT exactly reproduces the analytic cascaded Kerr -beta^2 Dw/(Dw^2+g2^2/4)
      (b) The Kerr coefficient changes sign as Dw -> 0 (sign tunability)
      (c) The cascaded TPA peaks on resonance (Lorentzian)
      (d) The 'naive' adiabatic formula -beta^2/Dw deviates near Dw=0
    """
    dw_arr = np.concatenate([np.linspace(-10, -0.3, 25),
                             np.linspace(0.3, 10, 25)])
    alpha_kerr_meas = np.zeros_like(dw_arr)
    alpha_loss_meas = np.zeros_like(dw_arr)

    for i, dw in enumerate(dw_arr):
        a1, a2, _ = steady_state_cw(beta, gamma1, gamma2, 0.0, dw, 0.0, s_amp)
        coef = (1j * np.conj(beta) * np.conj(a1) * a2) / (abs(a1)**2 * a1)
        alpha_kerr_meas[i] = coef.imag
        alpha_loss_meas[i] = -coef.real

    dw_fine_neg = np.linspace(-10, -0.05, 400)
    dw_fine_pos = np.linspace(0.05, 10, 400)
    dw_fine_full = np.linspace(-10, 10, 800)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))

    ax = axes[0]
    ax.plot(dw_fine_full/gamma2,
            alpha3_casc_full(beta, gamma2, dw_fine_full)/beta**2,
            'k-', lw=1.5,
            label=r'Full: $-\Delta\omega/(\Delta\omega^2+\gamma_2^2/4)$')
    for dw_seg in [dw_fine_neg, dw_fine_pos]:
        ax.plot(dw_seg/gamma2,
                alpha3_casc_adiabatic(beta, dw_seg)/beta**2,
                '--', color='gray', lw=1,
                label=(r'Adiabatic: $-1/\Delta\omega$' if dw_seg is dw_fine_neg else None))
    ax.plot(dw_arr/gamma2, alpha_kerr_meas/beta**2, 'o', color='C0',
            ms=5.5, mfc='white', mew=1.3, label='TCMT (small-signal)')
    ax.set_xlabel(r'Detuning  $\Delta\omega / \gamma_2$')
    ax.set_ylabel(r'$\alpha_3^{\rm Kerr}/|\beta|^2$')
    ax.set_title(r'(a) Cascaded Kerr — sign-tunable, vanishes on resonance')
    ax.set_ylim(-1.3, 1.3)
    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(0, color='k', lw=0.5)
    ax.legend(loc='upper right', framealpha=0.95)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(dw_fine_full/gamma2,
            alpha3_loss_casc(beta, gamma2, dw_fine_full)/beta**2,
            'k-', lw=1.5,
            label=r'Analytic: $(\gamma_2/2)/(\Delta\omega^2+\gamma_2^2/4)$')
    ax.plot(dw_arr/gamma2, alpha_loss_meas/beta**2, 's', color='C3',
            ms=5.5, mfc='white', mew=1.3, label='TCMT (small-signal)')
    ax.set_xlabel(r'Detuning  $\Delta\omega / \gamma_2$')
    ax.set_ylabel(r'$\gamma_{2{\rm PA}}^{\rm eff}/|\beta|^2$')
    ax.set_title(r'(b) Cascaded TPA — Lorentzian peak at $\Delta\omega = 0$')
    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(0, color='k', lw=0.5)
    ax.legend(loc='upper right', framealpha=0.95)
    ax.grid(alpha=0.3)

    fig.suptitle(r'Validation: TCMT (numerical) vs analytic cascaded coefficients   '
                 fr'($\beta={beta},\,\gamma_1/\gamma_2={gamma1/gamma2:g},\,|s_+|={s_amp}$)',
                 fontsize=11, y=1.02)
    fig.tight_layout()
    if savepath: fig.savefig(savepath); print(f'  saved {savepath}')
    return fig


# ============================================================
#   FIGURE 2: Bandwidth-enhancement sum rule (vary gamma_2)
# ============================================================

def figure2_sum_rule(beta=1.0, dw_ref=2.0, savepath=None):
    """
    Linearise around CW operating point.  For a small AC modulation of |a_1|^2
    at frequency Omega, the effective complex Kerr coefficient is

        alpha_3_eff(Omega) = i |beta|^2 / [gamma_2/2 + i (Omega - Delta_omega)]

    so |alpha_3_eff(Omega)| = |beta|^2 / sqrt[(gamma_2/2)^2 + (Omega-Delta_omega)^2]
    is a (root-)Lorentzian centred at Omega = Delta_omega with bandwidth ~ gamma_2.

    Peak |alpha_3_eff| = 2|beta|^2/gamma_2,  FWHM = sqrt(3)*gamma_2,
    so peak * FWHM = 2 sqrt(3) |beta|^2 (independent of gamma_2).

    Plot: |alpha_3_eff(Omega)| for several gamma_2.  As gamma_2 shrinks
    (Q_2 grows): peak grows, BW shrinks, product is invariant.
    """
    gamma2_list = [0.5, 1.0, 2.0, 4.0]
    colors = plt.cm.viridis(np.linspace(0.10, 0.85, len(gamma2_list)))

    Omega = np.linspace(-6, 10, 4000)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))

    ax = axes[0]
    peaks, fwhms, products = [], [], []
    for gamma2, color in zip(gamma2_list, colors):
        denom = (gamma2/2)**2 + (Omega - dw_ref)**2
        absH = beta**2 / np.sqrt(denom)
        ax.plot(Omega, absH/beta**2, color=color, lw=1.7,
                label=fr'$\gamma_2 = {gamma2:g}$  ($Q_2 \propto {1/gamma2:g}$)')
        peaks.append(2*beta**2 / gamma2)
        fwhms.append(np.sqrt(3) * gamma2)
        products.append(2*np.sqrt(3) * beta**2)
    ax.axvline(dw_ref, color='gray', ls=':', alpha=0.6,
               label=fr'$\Omega = \Delta\omega = {dw_ref}$')
    ax.set_xlabel(r'Modulation frequency  $\Omega$  (units of $\gamma_2^{\rm ref}=1$)')
    ax.set_ylabel(r'$|\alpha_3^{\rm eff}(\Omega)|/|\beta|^2$')
    ax.set_title(r'(a) AC response: smaller $\gamma_2$ → taller, narrower peak')
    ax.legend(loc='upper right', framealpha=0.95)
    ax.grid(alpha=0.3)

    ax = axes[1]
    bar_x = np.arange(len(gamma2_list))
    width = 0.35
    ax.bar(bar_x - width/2, np.array(peaks)/beta**2, width,
           color=colors, edgecolor='k', label='Peak  $|\\alpha_3^{\\rm eff}|/|\\beta|^2$')
    ax.bar(bar_x + width/2, np.array(fwhms),         width,
           color=colors, edgecolor='k', alpha=0.45,
           label=r'FWHM $\Delta\Omega$')
    ax2 = ax.twinx()
    ax2.plot(bar_x, np.array(products)/beta**2, 'k-D', ms=8, lw=1.5,
             label=r'Peak $\times$ FWHM $/|\beta|^2$')
    ax2.set_ylim(0, 5)
    ax2.set_ylabel(r'Peak $\times$ FWHM  /  $|\beta|^2$', color='k')
    ax.set_xticks(bar_x)
    ax.set_xticklabels([fr'$\gamma_2={g:g}$' for g in gamma2_list])
    ax.set_ylabel('Peak / FWHM individually')
    ax.set_title(r'(b) Sum rule: peak grows as $1/\gamma_2$, BW shrinks as $\gamma_2$')
    ax.legend(loc='upper left', framealpha=0.95)
    ax2.legend(loc='upper right', framealpha=0.95)
    ax.grid(alpha=0.3, axis='y')
    ax2.axhline(2*np.sqrt(3), color='r', ls='--', lw=1.2, alpha=0.7)
    ax2.text(len(gamma2_list)-0.5, 2*np.sqrt(3)+0.15,
             r'$2\sqrt{3}\,|\beta|^2$  (universal)',
             color='r', ha='right', fontsize=9.5)

    fig.suptitle(r'Bandwidth–enhancement sum rule '
                 fr'(fixed $\Delta\omega^{{\rm ref}}={dw_ref:g},\,\beta=1$)',
                 fontsize=11, y=1.02)
    fig.tight_layout()
    if savepath: fig.savefig(savepath); print(f'  saved {savepath}')
    return fig


# ============================================================
#   FIGURE 3: Optical bistability (CW pump, off-resonance)
# ============================================================

def figure3_bistability(beta=1.5, gamma1=1.0, gamma2=1.0,
                        dw_SH_intrinsic=4.0, dw_FM=1.5,
                        smax=2.5, n_pts=70, savepath=None):
    """
    Sweep input amplitude up then down at fixed pump frequency
    (positive FM detuning dw_FM > 0, with cascaded Kerr negative for dw_SH > 0:
    the cavity is pulled DOWN by intensity; pumping ABOVE resonance is the
    bistability geometry).

    The hysteresis loop is the experimental signature of the cascaded
    chi^(3): self-phase modulation drives the cavity through bistability
    even though the pump is monochromatic.
    """
    # Effective SH detuning in the rotating-at-pump frame:
    dw_SH_eff = dw_SH_intrinsic + 2*dw_FM

    s_arr = np.linspace(0.01, smax, n_pts)
    a1sq_up = np.zeros(n_pts);   a2sq_up = np.zeros(n_pts)
    a1sq_dn = np.zeros(n_pts);   a2sq_dn = np.zeros(n_pts)

    # Up sweep — start from zero, follow the lower branch
    state = [0., 0., 0., 0.]
    for i, s in enumerate(s_arr):
        a1, a2, state = steady_state_cw(beta, gamma1, gamma2, dw_FM, dw_SH_eff,
                                        0.0, s, t_settle=80.0/min(gamma1, gamma2),
                                        y0=state)
        a1sq_up[i] = abs(a1)**2; a2sq_up[i] = abs(a2)**2

    # Down sweep — start from the upper-branch endpoint
    state_dn = list(state)
    for j, s in enumerate(s_arr[::-1]):
        a1, a2, state_dn = steady_state_cw(beta, gamma1, gamma2, dw_FM, dw_SH_eff,
                                           0.0, s, t_settle=80.0/min(gamma1, gamma2),
                                           y0=state_dn)
        a1sq_dn[n_pts-1-j] = abs(a1)**2; a2sq_dn[n_pts-1-j] = abs(a2)**2

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))

    ax = axes[0]
    ax.plot(s_arr**2, a1sq_up, '-o', color='C0', ms=4, mfc='white', mew=1.0,
            label='Sweep up')
    ax.plot(s_arr**2, a1sq_dn, '-s', color='C3', ms=4, mfc='white', mew=1.0,
            label='Sweep down')
    # Mark the jumps
    iup_jump = int(np.argmax(np.diff(a1sq_up)))
    idn_jump = int(np.argmax(-np.diff(a1sq_dn)))
    ax.axvline(s_arr[iup_jump]**2, color='C0', ls=':', alpha=0.5)
    ax.axvline(s_arr[idn_jump]**2, color='C3', ls=':', alpha=0.5)
    ax.set_xlabel(r'Input power  $|s_+|^2$')
    ax.set_ylabel(r'Cavity intensity  $|a_1|^2$')
    ax.set_title(r'(a) Hysteresis loop — bistable response at FM')
    ax.legend(loc='upper left', framealpha=0.95)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(s_arr**2, a2sq_up, '-o', color='C0', ms=4, mfc='white', mew=1.0,
            label='Sweep up')
    ax.plot(s_arr**2, a2sq_dn, '-s', color='C3', ms=4, mfc='white', mew=1.0,
            label='Sweep down')
    ax.set_xlabel(r'Input power  $|s_+|^2$')
    ax.set_ylabel(r'SH intensity  $|a_2|^2$')
    ax.set_title(r'(b) SH amplitude jumps in lockstep with FM')
    ax.legend(loc='upper left', framealpha=0.95)
    ax.grid(alpha=0.3)

    fig.suptitle(r'Optical bistability driven by cascaded $\chi^{(3)}$  '
                 fr'($\beta={beta},\,\Delta\omega={dw_SH_intrinsic:g},\,\Delta_p={dw_FM:g}$)',
                 fontsize=11, y=1.02)
    fig.tight_layout()
    if savepath: fig.savefig(savepath); print(f'  saved {savepath}')
    return fig


# ============================================================
#   FIGURE 4: Saturation of small-signal description
# ============================================================

def figure4_saturation(beta=2.0, gamma1=1.0, gamma2=1.0, dw_SH=4.0,
                       s_max=2.5, n_pts=70, savepath=None):
    """
    Saturation: as the input grows, the cascaded back-action becomes large
    enough that |a_1|^2 departs from the linear-cavity prediction
    (4|s_+|^2 / gamma_1 at on-resonance pumping), and the SH starts to carry
    a non-negligible fraction of the cavity energy.

    The CW steady-state ratio |a_2|^2 / |a_1|^4 is exactly |beta|^2/|z|^2
    independent of intensity (this is *not* an approximation), so the right
    place to look for saturation is in the FM build-up itself: the cascaded
    Kerr (real) self-detunes the cavity off resonance, and the cascaded TPA
    (imaginary) adds an intensity-dependent loss.

    Pump on FM resonance (Delta_p = 0) so any departure from linear is a
    pure cascaded-chi3 effect, not a simple detuning artifact.
    """
    s_arr = np.linspace(0.01, s_max, n_pts)
    a1sq = np.zeros(n_pts); a2sq = np.zeros(n_pts)

    state = [0., 0., 0., 0.]
    for i, s in enumerate(s_arr):
        a1, a2, state = steady_state_cw(beta, gamma1, gamma2, 0.0, dw_SH, 0.0, s,
                                        t_settle=120.0/min(gamma1, gamma2),
                                        y0=state)
        a1sq[i] = abs(a1)**2; a2sq[i] = abs(a2)**2

    # Linear cavity (no cascade) at on-resonance pumping: |a_1|^2 = 4|s_+|^2 / gamma_1
    a1sq_linear = 4.0 / gamma1 * s_arr**2

    # Energy partition into SH
    sh_fraction = a2sq / (a1sq + a2sq + 1e-20)

    # Saturation parameter |beta a_1|/|Delta_omega|
    sat_param = beta * np.sqrt(a1sq) / abs(dw_SH)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))

    # Panel (a): FM intracavity intensity vs input, with linear-cavity reference
    ax = axes[0]
    ax.plot(s_arr**2, a1sq_linear, '--', color='k', lw=1.3,
            label=r'Linear cavity: $4|s_+|^2/\gamma_1$')
    ax.plot(s_arr**2, a1sq, '-o', color='C0', ms=4.5, mfc='white', mew=1.0,
            label='TCMT (with cascade)')
    ax.set_xlabel(r'Input power  $|s_+|^2$')
    ax.set_ylabel(r'$|a_1|^2$ (FM intracavity intensity)')
    ax.set_title('(a) Cascaded back-action saturates the FM build-up')
    ax.legend(loc='upper left', framealpha=0.95)
    ax.grid(alpha=0.3)

    # Panel (b): SH energy fraction + saturation parameter on twin axis
    ax = axes[1]
    ax.plot(s_arr**2, sh_fraction, '-s', color='C3', ms=4.5, mfc='white', mew=1.0,
            label=r'SH energy fraction  $|a_2|^2/(|a_1|^2+|a_2|^2)$')
    ax.set_xlabel(r'Input power  $|s_+|^2$')
    ax.set_ylabel('SH energy fraction', color='C3')
    ax.set_ylim(0, max(0.5, max(sh_fraction)*1.15))
    ax.tick_params(axis='y', labelcolor='C3')
    ax.grid(alpha=0.3)

    ax2 = ax.twinx()
    ax2.plot(s_arr**2, sat_param, '-^', color='C2', ms=4, mfc='white', mew=1.0,
             label=r'Saturation param $|\beta a_1|/|\Delta\omega|$')
    ax2.axhline(1.0, color='k', ls=':', lw=1.0, alpha=0.6)
    ax2.set_ylabel(r'$|\beta a_1|/|\Delta\omega|$', color='C2')
    ax2.tick_params(axis='y', labelcolor='C2')
    ax2.set_ylim(0, max(1.2, max(sat_param)*1.15))

    # Combined legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2,
              loc='upper left', framealpha=0.95, fontsize=9)
    ax.set_title('(b) SH grows; saturation parameter approaches unity')

    fig.suptitle(r'Saturation: departure from linear-cavity / small-signal cascaded $\chi^{(3)}$ '
                 fr'($\beta={beta},\,\gamma_1={gamma1:g},\,\Delta\omega={dw_SH:g},\,\Delta_p=0$)',
                 fontsize=11, y=1.02)
    fig.tight_layout()
    if savepath: fig.savefig(savepath); print(f'  saved {savepath}')
    return fig


# ============================================================
#   FIGURE 5: Pulse-domain bandwidth limit
# ============================================================

def figure5_pulse_bandwidth(beta=1.0, gamma1=5.0, gamma2=1.0, dw_SH=4.0,
                            s0=0.3, savepath=None):
    """
    Drive with Gaussian pulses of three different widths.  Compare actual
    a_2(t) against the adiabatic prediction a_2_ad = -i beta a_1^2 / (i*Dw - g_2/2).

    For sigma >> 1/gamma_2 (long pulse, narrow band): adiabatic works, SH tracks |a_1|^2.
    For sigma ~ 1/gamma_2: SH lags and overshoots.
    For sigma << 1/gamma_2 (short pulse, broad band): SH cannot track at all,
    response is heavily filtered.
    """
    sigmas = [0.3, 1.5, 8.0]
    titles = [r'Short  ($\sigma\gamma_2=0.3$): cascade cannot follow',
              r'Marginal ($\sigma\gamma_2=1.5$): partial tracking, ringing',
              r'Long  ($\sigma\gamma_2=8$): adiabatic, $a_2$ tracks $a_1^2$']

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 6.0), sharex='col')

    for col, (sigma, title) in enumerate(zip(sigmas, titles)):
        # Pulse: peak at t = 5 sigma, total window 12 sigma + tail for SH ringdown
        t_center = 5.0*sigma
        t_max = max(12.0*sigma, 12.0/gamma2 + t_center)
        s_func = lambda t, sig=sigma, tc=t_center: s0 * np.exp(-((t - tc)/sig)**2)

        sol = solve_ivp(
            tcmt_rhs, [0, t_max], [0., 0., 0., 0.],
            args=(beta, gamma1, gamma2, 0.0, dw_SH, 0.0, s_func),
            rtol=1e-10, atol=1e-13, method='DOP853', dense_output=True,
            max_step=min(sigma, 1.0/gamma2)/15,
        )
        t_plot = np.linspace(0, t_max, 1500)
        ys = sol.sol(t_plot)
        a1 = ys[0] + 1j*ys[1]
        a2 = ys[2] + 1j*ys[3]
        s_in_arr = np.array([s_func(t) for t in t_plot])

        # Adiabatic SH prediction
        a2_ad = -1j*beta*a1**2 / (1j*dw_SH - gamma2/2)

        # Top: input pulse shape and FM response
        ax = axes[0, col]
        ax.plot(t_plot * gamma2, np.abs(s_in_arr)**2 / s0**2, 'k-', lw=1.5,
                label=r'$|s_+|^2$ (norm.)')
        ax.plot(t_plot * gamma2, np.abs(a1)**2 / max(np.abs(a1)**2 + 1e-20),
                color='C0', lw=1.5, label=r'$|a_1|^2$ (norm.)')
        ax.set_title(title, fontsize=10.5)
        ax.set_ylabel('FM intensities (norm.)')
        if col == 2:
            ax.legend(loc='upper right', fontsize=9)
        ax.grid(alpha=0.3)

        # Bottom: SH actual vs adiabatic
        ax = axes[1, col]
        ax.plot(t_plot * gamma2, np.abs(a2)**2, color='C3', lw=1.5,
                label=r'TCMT $|a_2|^2$')
        ax.plot(t_plot * gamma2, np.abs(a2_ad)**2, '--', color='k', lw=1.0,
                label=r'Adiabatic $|a_2^{\rm ad}|^2$')
        ax.set_xlabel(r'$t \times \gamma_2$')
        ax.set_ylabel(r'$|a_2|^2$')
        if col == 2:
            ax.legend(loc='upper right', fontsize=9)
        ax.grid(alpha=0.3)

    fig.suptitle(r'Pulse-domain bandwidth limit: SH (and hence cascaded $\chi^{(3)}$) '
                 r'follows only when pulse spectral width $\lesssim \gamma_2$',
                 fontsize=11, y=1.00)
    fig.tight_layout()
    if savepath: fig.savefig(savepath); print(f'  saved {savepath}')
    return fig


# ============================================================
#                       MAIN
# ============================================================
if __name__ == '__main__':
    print('Generating figures (this may take a minute)...')

    print('  [1/5] Validation + sign tunability ...')
    figure1_validation_sign_tunability(savepath=os.path.join(OUT, 'fig1_validation.png'))

    print('  [2/5] Sum rule ...')
    figure2_sum_rule(savepath=os.path.join(OUT, 'fig2_sum_rule.png'))

    print('  [3/5] Bistability ...')
    figure3_bistability(savepath=os.path.join(OUT, 'fig3_bistability.png'))

    print('  [4/5] Saturation ...')
    figure4_saturation(savepath=os.path.join(OUT, 'fig4_saturation.png'))

    print('  [5/5] Pulse-domain bandwidth ...')
    figure5_pulse_bandwidth(savepath=os.path.join(OUT, 'fig5_pulse_bandwidth.png'))

    print('All figures saved to', OUT)
