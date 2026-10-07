"""
Group 2, Figure B2: Bistability threshold comparison
=====================================================

Demonstrates that the cascaded contribution alpha_3_casc lowers the
bistability threshold proportionally to (1+R)^-1.

Strategy: hold alpha_3_native fixed (smaller than the value required to
trigger bistability alone), then add cascade with increasing |beta|.
The combined alpha_3_eff = alpha_3_native + alpha_3_casc crosses the
bistability cusp at finite R, and from there the threshold drops as
1/(1+R).

Three panels:
  (a) hysteresis loops for several R values
  (b) snap-up threshold vs R, with 1/(1+R) prediction
  (c) intracavity SH energy at snap-up
"""

import numpy as np
import matplotlib.pyplot as plt
from core import steady_state_cw, alpha3_casc_full, adjacent_save


def sweep_input(beta, gamma1, gamma2, dw_FM, dw_SH, alpha3_native,
                s_arr, direction='up'):
    n = len(s_arr)
    a1sq = np.zeros(n)
    a2sq = np.zeros(n)
    state = [0., 0., 0., 0.]
    seq = range(n) if direction == 'up' else range(n - 1, -1, -1)
    for k in seq:
        a1, a2, state = steady_state_cw(
            beta, gamma1, gamma2, dw_FM, dw_SH, alpha3_native, s_arr[k],
            t_settle=60.0 / min(gamma1, gamma2), y0=state)
        a1sq[k] = abs(a1)**2
        a2sq[k] = abs(a2)**2
    return a1sq, a2sq


def find_snap_threshold(s_arr, a1sq_up):
    """Largest jump in up-sweep."""
    di = np.diff(a1sq_up)
    idx = int(np.argmax(di))
    return s_arr[idx + 1]**2, idx


def figure_bistability_comparison(savepath=None):
    gamma1 = 1.0
    gamma2 = 1.0
    # alpha_3 > 0 (typical native dielectric Kerr): need delta_p < 0 (red-detuned pump)
    dw_FM = -2.5
    # For constructive cascade with alpha_3_casc > 0: need Delta_SH < 0
    # Make |Delta_SH| LARGE to keep |beta a_1|/|Delta_SH| << 1 (no parametric oscillation)
    dw_SH = -10.0   # deeper adiabatic
    # alpha_3 large enough that A *alone* is bistable
    alpha3_native = 2.0
    # Cascade beta: pick values keeping |alpha_3_casc| moderate
    beta_values = [0.0, 2.0, 3.0, 4.0]
    R_values = [abs(alpha3_casc_full(b, gamma2, dw_SH)) / alpha3_native
                for b in beta_values]

    s_max = 5.0
    n_pts = 120
    s_arr = np.linspace(0.05, s_max, n_pts)

    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(beta_values)))

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # Panel (a): hysteresis loops
    ax = axes[0]
    thresholds = []
    for i, (beta, R, c) in enumerate(zip(beta_values, R_values, colors)):
        up, _ = sweep_input(beta, gamma1, gamma2, dw_FM, dw_SH,
                            alpha3_native, s_arr, 'up')
        dn, _ = sweep_input(beta, gamma1, gamma2, dw_FM, dw_SH,
                            alpha3_native, s_arr, 'down')
        ax.plot(s_arr**2, up, '-', color=c, lw=1.5,
                label=fr'$\beta={beta:.1f}$, $R={R:.2f}$')
        ax.plot(s_arr**2, dn, '--', color=c, lw=1.5, alpha=0.6)
        sthr, _ = find_snap_threshold(s_arr, up)
        thresholds.append(sthr)
    ax.set_xlabel(r'Input power $|s_+|^2$')
    ax.set_ylabel(r'Intracavity FM intensity $|a_1|^2$')
    ax.set_title(r'(a) Hysteresis loops: increasing $\beta$ (i.e. $R$)')
    ax.legend(loc='upper left', fontsize=9, framealpha=0.95)
    ax.grid(True, alpha=0.3)

    # Panel (b): threshold vs R, prediction 1/(1+R) relative to R=0
    ax = axes[1]
    sthr_native = thresholds[0]
    ratios_measured = np.array(thresholds) / sthr_native
    R_fine = np.linspace(0, max(R_values) * 1.2, 100)
    ratios_predicted = 1.0 / (1.0 + R_fine)
    ax.plot(R_values, ratios_measured, 'o', color='C3', ms=8, mfc='white', mew=2,
            label='TCMT (measured)')
    ax.plot(R_fine, ratios_predicted, '-', color='k', lw=1.5,
            label=r'Predicted: $1/(1+R)$')
    ax.set_xlabel(r'Cascade ratio $R$')
    ax.set_ylabel(r'$|s_+|^{2,\,{\rm A+B}}_{\rm bist}\,/\,|s_+|^{2,\,{\rm A}}_{\rm bist}$')
    ax.set_title(r'(b) Threshold reduction follows $1/(1+R)$')
    ax.legend(loc='upper right', framealpha=0.95)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)

    # Panel (c): absolute threshold vs R, with predicted curve
    ax = axes[2]
    threshold_predicted = sthr_native / (1.0 + R_fine)
    ax.plot(R_values, thresholds, 'o', color='C3', ms=8, mfc='white', mew=2,
            label='TCMT (measured)')
    ax.plot(R_fine, threshold_predicted, '-', color='k', lw=1.5,
            label=r'Predicted: $|s_+|^{2,\,{\rm A}}_{\rm bist}/(1+R)$')
    ax.axhline(sthr_native, color='C0', ls=':', alpha=0.6,
               label=r'Native-only threshold')
    ax.set_xlabel(r'Cascade ratio $R$')
    ax.set_ylabel(r'Snap-up threshold $|s_+|^2_{\rm bist}$')
    ax.set_title(r'(c) Absolute switching energy budget')
    ax.legend(loc='upper right', framealpha=0.95)
    ax.grid(True, alpha=0.3)

    fig.suptitle(r'Bistability: cascade lowers the threshold/switching-energy as $1/(1+R)$  '
                 fr'(constructive sign, $\delta_p={dw_FM},\,\Delta_{{\rm SH}}={dw_SH}$)',
                 fontsize=11, y=1.02)
    fig.tight_layout()
    if savepath:
        adjacent_save(fig, savepath)
    return fig


if __name__ == '__main__':
    figure_bistability_comparison(savepath='G2_B2_bistability_comparison.png')
