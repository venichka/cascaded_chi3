"""Time-domain TCMT pulsed simulation for the SH-only metasurface.

The CW analysis (sample_effects.py) evaluates the steady-state cascaded chi3 at the
pulse PEAK intensity -- valid only when the pump envelope is slow compared with the
SH lifetime 1/gamma2. For a ps pump that assumption fails: tau_p ~ 1-3 ps is
comparable to 1/gamma2 ~ 0.9-1.6 ps for these samples, so the SH cannot fully build
up and the true cascaded response is reduced. This module integrates the SH-only
TCMT in the time domain to quantify that reduction.

Model (note eq:a2-1R), with the FM responding instantaneously (broadband, no FM
resonance) so the drive follows the pulse intensity envelope:

    da2/dt = (i*Domega - gamma2/2) a2 + i*beta(t),   beta(t) ∝ |E_loc(t)|^2 ∝ I(t)

Dimensionless (tau = gamma2 t, delta = Domega/gamma2, drive b normalised to peak 1):

    da2/dtau = (i*delta - 1/2) a2 + i*b(tau),   b(tau) = exp(-4 ln2 (tau/W)^2),
    W = gamma2 * tau_p   (pulse intensity FWHM in units of 1/gamma2)

The cascaded back-polarisation is P_casc(t) ∝ a2(t) E_loc*(t), so the instantaneous
effective susceptibility is chi3_eff(t) ∝ a2(t)/E_loc(t)^2 ∝ a2(t)/b(t). Normalising
to the CW limit (a2_ss/b = i/(1/2 - i*delta)) defines the transient factor

    F(tau) = [a2(tau)/b(tau)] * (1/2 - i*delta)/i      (-> 1 in CW steady state),

and chi3_eff(tau) = F(tau) * chi3_CW(delta). The measurable, intensity-weighted
nonlinear phase relative to the quasi-static (CW-at-peak) value is

    r_phi(delta, W) = Int Re[F(tau) chi3hat] I(tau) dtau / ( Re[chi3hat] Int I dtau ),
    chi3hat(delta) = -1/(delta + i/2)   (dimensionless CW lineshape; real consts cancel).

Assumes a transform-limited (real-envelope) pulse. Run to generate figures + table.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from style import apply_style
apply_style()

from sample_params import build_sample
from sample_effects import SAMPLE_COLORS, _save, DOP

FOURLN2 = 4.0 * np.log(2.0)
_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz   # numpy 2.x renamed trapz


def chi3hat(delta):
    """Dimensionless CW cascaded lineshape (real prefactors drop out of ratios)."""
    return -1.0 / (delta + 0.5j)


def solve_a2(delta, W, n=8000, span=8.0):
    """Integrate the SH amplitude for a Gaussian-intensity pulse of FWHM W (=gamma2*tau_p).

    Returns (tau, b, a2) on a uniform grid spanning +/- span*W around the pulse peak.
    """
    tmax = max(span * W, 12.0)
    tau = np.linspace(-tmax, tmax, n)
    b = np.exp(-FOURLN2 * (tau / W) ** 2)

    def drive(t):
        return np.exp(-FOURLN2 * (t / W) ** 2)

    def rhs(t, y):
        x, yy = y
        bt = drive(t)
        return [-delta * yy - 0.5 * x, delta * x - 0.5 * yy + bt]

    sol = solve_ivp(rhs, (tau[0], tau[-1]), [0.0, 0.0], t_eval=tau,
                    method="DOP853", rtol=1e-9, atol=1e-12)
    a2 = sol.y[0] + 1j * sol.y[1]
    return tau, b, a2


def a2_ss(delta):
    """CW steady-state a2 at unit drive."""
    return 1j / (0.5 - 1j * delta)


def reduction_factors(delta, W):
    """Return (r_phi, r_peak, sh_buildup, eff_kerr) for detuning delta, width W=gamma2*tau_p.

    Note: the *energy* (I-weighted) nonlinear phase is conserved for any pulse length
    (a sum rule -- the SH delivers all its action, just spread over ~1/gamma2). The
    physically reduced quantity, and what a closed-aperture Z-scan with an energy
    detector measures, is the I^2-weighted phase (phase ∝ I, detection weight ∝ I):

        eff_kerr = Int Re[c a2(t)] b(t) dt / Int b(t)^2 dt ,  c = (1/2 - i*delta)/i * chi3hat,
        r_phi    = eff_kerr / Re[chi3hat]                  (-> 1 in CW; < 1 for short pulses).

    r_peak     : instantaneous Kerr fraction Re[chi3_eff(0)]/Re[chi3_CW] at the pulse peak.
    sh_buildup : peak SH energy |a2|^2 / CW |a2_ss|^2 (can exceed 1 from resonant ringing).
    eff_kerr   : effective (I^2-weighted) Re chi3 in CW-lineshape units (for the sign-flip plot).
    """
    tau, b, a2 = solve_a2(delta, W)
    ch = chi3hat(delta)
    c = (0.5 - 1j * delta) / 1j * ch          # so that c * a2_ss = chi3hat
    b2 = _trapz(b * b, tau)
    eff_kerr = _trapz(np.real(c * a2) * b, tau) / b2
    rch = np.real(ch)
    r_phi = eff_kerr / rch if abs(rch) > 1e-9 else np.nan
    i0 = int(np.argmin(np.abs(tau)))           # pulse centre
    F0 = (a2[i0] / b[i0]) * (0.5 - 1j * delta) / 1j
    r_peak = np.real(F0 * ch) / rch if abs(rch) > 1e-9 else np.nan
    sh = np.max(np.abs(a2) ** 2) / np.abs(a2_ss(delta)) ** 2
    return r_phi, r_peak, sh, eff_kerr


# ----------------------------------------------------------------------
def fig_timedomain(samples, delta=DOP):
    """SH lag (left) and the delayed/reduced phase accumulation (right) for S4."""
    p = samples[4]
    inv_g_ps = 1e12 / p["gamma2"]                       # 1/gamma2 in ps
    ch = chi3hat(delta)
    c = (0.5 - 1j * delta) / 1j * ch
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

    # left: SH energy build-up vs the pulse, tau_p = 1, 2, 3 ps
    for tp, ls in zip((1.0, 2.0, 3.0), ("-", "--", ":")):
        W = p["gamma2"] * tp * 1e-12
        tau, b, a2 = solve_a2(delta, W)
        t_ps = tau / p["gamma2"] * 1e12
        sh = np.abs(a2) ** 2 / np.abs(a2_ss(delta)) ** 2
        axes[0].plot(t_ps, sh, ls, color=SAMPLE_COLORS[4], label=fr"$\tau_p={tp:.0f}$ ps ($W={W:.2f}$)")
    W3 = p["gamma2"] * 3e-12
    tau, b, _ = solve_a2(delta, W3)
    axes[0].plot(tau / p["gamma2"] * 1e12, b, color="0.6", lw=1.0, label="drive $I(t)$, 3 ps")
    axes[0].axhline(1.0, color="0.4", lw=0.8, ls="--")
    axes[0].set_xlim(-12, 14)
    axes[0].text(-11.5, 1.04, "CW level", fontsize=8, color="0.4")
    axes[0].set_xlabel("time [ps]"); axes[0].set_ylabel(r"SH energy $|a_2|^2/|a_2^{\rm ss}|^2$")
    axes[0].set_title(fr"SH build-up lags the pulse (S4, $1/\gamma_2={inv_g_ps:.2f}$ ps)")
    axes[0].legend(fontsize=8.5, loc="upper right")

    # right: Z-scan phase-accumulation density (no 1/b division -> no edge blow-up)
    for tp, col, lab in zip((1.0, 3.0), (SAMPLE_COLORS[4], "#2ca02c"), ("1", "3")):
        W = p["gamma2"] * tp * 1e-12
        tau, b, a2 = solve_a2(delta, W)
        t_ps = tau / p["gamma2"] * 1e12
        rho_p = np.abs(np.real(c * a2) * b)             # pulsed phase density
        rho_qs = np.abs(np.real(ch) * b ** 2)           # quasi-static (∝ I^2)
        nrm = rho_qs.max()
        axes[1].plot(t_ps, rho_qs / nrm, color=col, lw=1.2, ls="--", alpha=0.7)
        axes[1].plot(t_ps, rho_p / nrm, color=col, lw=2.0, label=fr"$\tau_p={lab}$ ps")
    axes[1].set_xlim(-8, 12)
    axes[1].set_xlabel("time [ps]")
    axes[1].set_ylabel(r"Z-scan phase density (norm.)")
    axes[1].set_title("Phase: pulsed (solid) lags + below quasi-static (dashed)")
    axes[1].legend(fontsize=8.5, loc="upper right")
    _save(fig, "sample_pulsed_timedomain")


def fig_reduction(samples, delta=DOP):
    """Pulse phase reduction r_phi vs tau_p for all samples; mark 1-3 ps band."""
    tp = np.linspace(0.3, 6.0, 120) * 1e-12
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    rows = {}
    for i, p in samples.items():
        W = p["gamma2"] * tp
        r = np.array([reduction_factors(delta, w)[0] for w in W])
        ax.plot(tp * 1e12, r, color=SAMPLE_COLORS[i],
                label=fr"S{i} ($1/\gamma_2={1e12/p['gamma2']:.2f}$ ps)")
        rows[i] = np.interp([1, 2, 3], tp * 1e12, r)
    ax.axvspan(1.0, 3.0, color="#1f77b4", alpha=0.10)
    ax.text(2.0, 0.06, "ps-laser\n(1--3 ps)", ha="center", fontsize=8.5, color="#1f77b4")
    ax.axhline(1.0, color="0.4", lw=1.0, ls="--", label="CW / quasi-static")
    ax.set_xlabel(r"pulse duration $\tau_p$ (intensity FWHM) [ps]")
    ax.set_ylabel(r"phase reduction $r_\varphi=\varphi_{\rm pulsed}/\varphi_{\rm quasi\text{-}static}$")
    ax.set_title(fr"Transient suppression of the cascaded Kerr phase ($\Delta\omega={delta:.0f}\gamma_2$)")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=8.5, loc="lower right")
    _save(fig, "sample_pulsed_reduction")
    return rows


def fig_signflip(samples):
    """Pulse-integrated effective Re chi3 vs detuning for tau_p=1,2,3 ps (S4) vs CW."""
    p = samples[4]
    delta = np.linspace(-6, 6, 241)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    # CW Re chi3 (normalised lineshape)
    cw = np.real(chi3hat(delta))
    norm = np.max(np.abs(cw))
    ax.plot(delta, cw / norm, color="0.4", lw=2, ls="--", label="CW (quasi-static)")
    for tp, col in zip((1.0, 2.0, 3.0), ("#d62728", "#ff7f0e", "#2ca02c")):
        W = p["gamma2"] * tp * 1e-12
        eff = []
        for d in delta:
            eff.append(reduction_factors(d, W)[3])   # I^2-weighted effective Re chi3
        ax.plot(delta, np.array(eff) / norm, color=col, label=fr"pulsed $\tau_p={tp:.0f}$ ps")
    ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax.set_ylabel(r"$\mathrm{Re}\,\chi^{(3)}_{\rm eff}$ (normalised to CW peak)")
    ax.set_title(r"Sign flip survives pulsing (S4) -- reduced and broadened")
    ax.legend(fontsize=8.5, loc="upper right")
    _save(fig, "sample_pulsed_signflip")


def report(samples, delta=DOP):
    print("=== TIME-DOMAIN PULSED TCMT  (delta = %.0f gamma2) ===" % delta)
    print("SH lifetime 1/gamma2 [ps]:  " + ", ".join("S%d %.2f" % (i, 1e12 / samples[i]["gamma2"]) for i in samples))
    print("  (r_phi = I^2-weighted / Z-scan phase reduction; r_peak = peak-instant Kerr fraction)")
    print("\n%-2s %8s | %7s %7s %7s | %7s %7s %7s" %
          ("S", "1/g2[ps]", "rphi@1", "rphi@2", "rphi@3", "rpk@1", "rpk@2", "rpk@3"))
    for i, p in samples.items():
        rphi, rpk = [], []
        for tp in (1.0, 2.0, 3.0):
            W = p["gamma2"] * tp * 1e-12
            a, b_, _, _ = reduction_factors(delta, W)
            rphi.append(a); rpk.append(b_)
        print("%-2d %8.2f | %7.3f %7.3f %7.3f | %7.3f %7.3f %7.3f"
              % (i, 1e12 / p["gamma2"], *rphi, *rpk))


def main():
    samples = {i: build_sample(i) for i in (1, 2, 3, 4)}
    report(samples)
    fig_timedomain(samples)
    fig_reduction(samples)
    fig_signflip(samples)


if __name__ == "__main__":
    main()
