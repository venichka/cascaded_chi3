"""Cross-phase modulation (XPM): cascade vs native, verifying (1+R)
ratio and the factor 2 over SPM.

Tier 1, item 2 of the note's Section 13. A weak signal at omega_1 +
Omega (distinguishable from the pump by frequency offset
|Omega| < gamma_1) acquires a Kerr-induced phase shift driven by
|a_1^pump|^2. The cascade picture mediates this through pump->SH
->DFG-with-signal; the native picture goes directly through chi^(3).

Predicted scalings (Boyd convention, scalar TCMT):
  - alpha_3^XPM = 2 * alpha_3^SPM (factor 2 from index-permutation
    degeneracy with one distinguishable field).
  - alpha_3^{XPM,A+B} / alpha_3^{XPM,A} = 1 + R, same as SPM ratio.

Implementation. Same bichromatic infrastructure as the FWM figure:
drive at pump + weak signal, integrate to quasi-periodic steady state,
extract the +Omega Fourier component of a_1. The phase of that
component relative to the linear-cavity prediction is the XPM phase
shift on the signal.

Two panels:
  (a) XPM phase shift on the signal vs pump intensity, cascade vs native
      vs linear cavity. Slope is unity in log-log (linear in pump
      intensity at small signal).
  (b) Ratio (1+R) between cascade-XPM and native-XPM, plus the
      factor 2 between XPM (this fig) and SPM (recomputed at the same
      parameters).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from params import TCMTParams
from tcmt import integrate, steady_state

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def bichromatic_pump(s_p, s_s, Omega):
    def f(t):
        return s_p + s_s * np.exp(1j * Omega * t)
    return f


def extract_signal_phase(p: TCMTParams, s_p: float, s_s: float,
                         Omega: float, t_settle: float = 400.0,
                         n_periods: int = 30, n_per_period: int = 60):
    """Integrate the bichromatic TCMT; return the +Omega component of
    a_1 (referenced to t=0) and the linear-cavity prediction. The
    difference in phase between them is the XPM-induced phase shift."""
    T_period = 2 * np.pi / Omega
    t_total = t_settle + n_periods * T_period
    N = n_periods * n_per_period
    t_eval = t_settle + np.arange(N) * T_period * n_periods / N
    sp = bichromatic_pump(s_p, s_s, Omega)
    sol = integrate(p, sp, (0.0, t_total), t_eval=t_eval,
                    max_step=min(0.5, T_period / 30))
    a1 = sol.a[0]
    # +Omega bin of a_1 (FFT references the window's start time t_settle;
    # multiply by exp(+i*Omega*t_settle) to get the rotating-frame
    # amplitude in the same phase convention as a_lin below).
    A = np.fft.fft(a1) / N
    # FFT bin at m = n_periods gives C_1 * exp(+i Omega t_settle) where
    # C_1 is the rotating-frame +Omega Fourier coefficient referenced to
    # t=0. Strip the exp factor to recover C_1 itself.
    a_sig = A[n_periods] * np.exp(-1j * Omega * t_eval[0])
    # Linear cavity response at +Omega in the rotating frame:
    # d a_s/dt = -(gamma_1/2) a_s + sqrt(gamma_rad) s_s, with a_s
    # oscillating as exp(i Omega t) -- in steady state
    # a_s = sqrt(gamma_rad) s_s / (gamma_1/2 + i Omega).
    a_lin = np.sqrt(p.gamma1_rad) * s_s / (0.5 * p.gamma1 + 1j * Omega)
    return a_sig, a_lin


def main():
    # alpha_3 = +0.04 with Delta omega = -2 (constructive cascade: same
    # sign as native), so the combined coefficient is alpha_3 +
    # alpha_3^casc and XPM/SPM ratios give 1+R rather than |1-R|.
    # R = 0.471 / 0.04 = 11.8 matches the worked example.
    base = dict(gamma1=10.0, gamma1_rad=5.0, alpha3=0.04, delta1=0.0)
    p_nat = TCMTParams(Domega=-2.0, beta=0.0, **base)
    p_cas = TCMTParams(Domega=-2.0, beta=1.0, **base)

    Omega = 0.5      # signal offset within FM cavity (gamma_1=10)
    s_s = 1e-3       # weak signal

    s_p_grid = np.geomspace(0.01, 0.5, 16)
    phi_XPM_nat = np.empty_like(s_p_grid)
    phi_XPM_cas = np.empty_like(s_p_grid)
    phi_SPM_nat = np.empty_like(s_p_grid)
    phi_SPM_cas = np.empty_like(s_p_grid)

    for i, s_p in enumerate(s_p_grid):
        # XPM: phase shift on signal at +Omega
        a_sig_n, a_lin = extract_signal_phase(p_nat, s_p, s_s, Omega)
        a_sig_c, _ = extract_signal_phase(p_cas, s_p, s_s, Omega)
        phi_XPM_nat[i] = np.angle(a_sig_n / a_lin)
        phi_XPM_cas[i] = np.angle(a_sig_c / a_lin)

        # SPM: phase shift on pump itself at zero-frequency (DC)
        # Re-use the cw steady_state result with single-tone pump.
        a1_n, _ = steady_state(p_nat, s_p)
        a1_c, _ = steady_state(p_cas, s_p)
        a_pump_lin = np.sqrt(p_nat.gamma1_rad) * s_p / (0.5 * p_nat.gamma1)
        phi_SPM_nat[i] = np.angle(a1_n / a_pump_lin)
        phi_SPM_cas[i] = np.angle(a1_c / a_pump_lin)

    # ---- Panel (a): XPM phase shift vs pump intensity ---------
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.7))

    axs[0].loglog(s_p_grid ** 2, np.abs(phi_XPM_nat), "o-", color="C2",
                  ms=5, mfc="none", mew=1.2,
                  label=r"native XPM")
    axs[0].loglog(s_p_grid ** 2, np.abs(phi_XPM_cas), "s-", color="C0",
                  ms=5, mfc="none", mew=1.2,
                  label=r"cascade XPM")
    # Slope +1 guide
    s_ref = np.array([s_p_grid[0] ** 2, s_p_grid[-1] ** 2])
    y0 = abs(phi_XPM_nat[0])
    axs[0].loglog(s_ref, y0 * s_ref / s_ref[0], ":", color="0.55", lw=0.9)
    axs[0].text(3e-4, 3e-7, "slope $+1$\n(small-signal)",
                fontsize=9, color="0.4")
    axs[0].set_xlabel(r"$|s_+^{(p)}|^2$  (pump power)")
    axs[0].set_ylabel(r"XPM phase shift on signal  $|\phi_{\rm XPM}|$  (rad)")
    axs[0].set_title("(a) XPM phase shift: cascade exceeds native by $\\sim(1+R)$")
    axs[0].grid(True, which="both", alpha=0.3)
    axs[0].legend(loc="lower right", fontsize=9.5)

    # ---- Panel (b): Architecture ratio + XPM/SPM ratio --------------
    ratio_arch = np.abs(phi_XPM_cas) / np.abs(phi_XPM_nat)
    ratio_XPM_SPM_nat = np.abs(phi_XPM_nat) / np.abs(phi_SPM_nat)
    ratio_XPM_SPM_cas = np.abs(phi_XPM_cas) / np.abs(phi_SPM_cas)

    R_predicted = 12
    axs[1].semilogx(s_p_grid ** 2, ratio_arch, "o-", color="C0", ms=5,
                    mfc="none", mew=1.2,
                    label=r"$|\phi_{\rm cas}^{\rm XPM}| / |\phi_{\rm nat}^{\rm XPM}|$ (TCMT)")
    axs[1].axhline(1 + R_predicted, color="C3", ls="--", lw=1.6,
                   label=rf"$1 + R = {1 + R_predicted}$")
    axs[1].semilogx(s_p_grid ** 2, ratio_XPM_SPM_nat, "^-", color="C2",
                    ms=5, mfc="none", mew=1.2,
                    label=r"$|\phi_{\rm nat}^{\rm XPM}|/|\phi_{\rm nat}^{\rm SPM}|$ (TCMT)")
    axs[1].semilogx(s_p_grid ** 2, ratio_XPM_SPM_cas, "v-", color="C1",
                    ms=5, mfc="none", mew=1.2,
                    label=r"$|\phi_{\rm cas}^{\rm XPM}|/|\phi_{\rm cas}^{\rm SPM}|$ (TCMT)")
    axs[1].axhline(2.0, color="0.4", ls=":", lw=1.0)
    axs[1].text(s_p_grid[0]**2 * 1.5, 2.1, r"factor $2$ (Boyd, XPM/SPM)",
                fontsize=9, color="0.4")
    axs[1].set_xlabel(r"$|s_+^{(p)}|^2$  (pump power)")
    axs[1].set_ylabel(r"ratio")
    axs[1].set_title(r"(b) $(1+R)$ cascade/native and factor-2 XPM/SPM")
    axs[1].set_ylim(0, max(35, ratio_arch.max() * 1.1))
    axs[1].grid(True, which="both", alpha=0.3)
    axs[1].legend(loc="upper left", fontsize=9, ncol=1)

    fig.suptitle(
        r"XPM: cascade $\times(1+R)$ over native; XPM $\times 2$ over SPM "
        r"-- $\Omega = 0.5\gamma_2$, $\widetilde\alpha_3 = +0.04$, "
        r"$\Delta\omega = -2\gamma_2$ (constructive; $R\simeq 12$)",
        fontsize=11.5, y=1.00,
    )
    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "xpm"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    print(f"  at |s+|^2 = {s_p_grid[0]**2:.3e}:")
    print(f"    XPM ratio cascade/native = {ratio_arch[0]:.2f}   (predicted 1+R = {1+R_predicted})")
    print(f"    XPM/SPM ratio (native)   = {ratio_XPM_SPM_nat[0]:.2f}  (predicted 2.0)")
    print(f"    XPM/SPM ratio (cascade)  = {ratio_XPM_SPM_cas[0]:.2f}  (predicted 2.0)")
    plt.close(fig)


if __name__ == "__main__":
    main()
