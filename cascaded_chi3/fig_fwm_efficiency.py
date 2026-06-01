"""FWM efficiency: cascade vs native, verifying the (1+R)^2 ratio.

Tier 1, item 4 of the note's Section 13. The cascade has a quadratic
leverage in R for FWM efficiency because the conversion amplitude is
proportional to the effective Kerr coefficient, and the efficiency is
that amplitude squared. This is the strongest single-observable test
of the cascade mechanism.

Setup. Drive with a bichromatic input at the FM port:
    s_+(t) = s_p + s_s exp(i Omega t)
in the rotating frame at omega_p = omega_1. Pump is at DC; signal at
Omega > 0 (offset by Omega from the FM resonance, with |Omega| < gamma_1
so it sits inside the cavity linewidth). FWM generates an idler at
-Omega:
    pump + pump - signal -> idler   (FWM phase-matching in the FM cavity)
mediated by either the bulk Kerr alpha_3 (native) or the cascade
beta -> SH -> DFG-back (cascade). For both architectures the FM
output s_-(t) inherits a component at -Omega; its modulus-squared,
normalised to the signal input power, is the conversion efficiency

    eta_FWM = |s_minus_(-Omega)|^2 / |s_+_(Omega)|^2.

The architecture ratio at small signal is (1+R)^2 from the note's
eq. for eta_FWM.

Implementation. We integrate the TCMT to a quasi-periodic steady state,
then FFT the last few modulation periods to extract the idler-frequency
component of s_-(t).
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
from tcmt import integrate

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def bichromatic_pump(s_p, s_s, Omega):
    def f(t):
        return s_p + s_s * np.exp(1j * Omega * t)
    return f


def extract_idler_amplitude(p: TCMTParams, s_p: float, s_s: float,
                            Omega: float, t_settle: float = 400.0,
                            n_periods: int = 20, n_per_period: int = 60):
    """Integrate the TCMT under bichromatic drive to quasi-periodic
    steady state, then FFT the last n_periods to extract the
    idler-frequency component of s_minus(t).

    Uses np.fft.fft on exactly n_periods * n_per_period samples
    (NO duplicate endpoint) so that the +Omega and -Omega Fourier
    bins are exactly orthogonal -- this is critical because the +Omega
    signal is ~10^4 larger than the -Omega idler and ANY spectral
    leakage drowns out the idler.
    """
    T_period = 2 * np.pi / Omega
    t_total = t_settle + n_periods * T_period
    N = n_periods * n_per_period
    # N evenly-spaced samples spanning EXACTLY n_periods periods, no
    # duplicate endpoint. dt = (n_periods * T) / N = T / n_per_period.
    t_eval = t_settle + np.arange(N) * T_period * n_periods / N

    sp = bichromatic_pump(s_p, s_s, Omega)
    sol = integrate(p, sp, (0.0, t_total), t_eval=t_eval,
                    max_step=min(0.5, T_period / 30))
    a1 = sol.a[0]
    t = sol.t
    s_in = np.array([sp(ti) for ti in t])
    s_out = -s_in + np.sqrt(p.gamma1_rad) * a1  # input-output relation

    # Discrete Fourier transform. With N samples over n_periods
    # periods, the integer-frequency bins are at k = ..., -n_periods,
    # 0, +n_periods, ... corresponding to Fourier components at
    # frequencies -Omega, 0, +Omega, ...
    S = np.fft.fft(s_out) / N
    idler_amp = S[-n_periods]                # bin at -Omega
    signal_amp = S[n_periods]                # bin at +Omega
    return idler_amp, signal_amp


def main():
    # FM cavity matched across architectures. Use Delta omega = -2 so
    # the cascade has the SAME sign as native alpha_3 > 0 (constructive
    # configuration): combined alpha_3^eff = alpha_3 + alpha_3^casc, and
    # FWM ~ |alpha_3^eff|^2 gives the clean (1+R)^2 ratio rather than the
    # cancellation case |1-R|^2.
    # alpha_3 = 0.04 = LiNbO3 native chi^(3) in dim. TCMT units
    # at the conservative design point |Delta omega| = 2 gamma_2,
    # Q_2 = 4000 -- yields full-Lorentzian R = 0.471/0.04 = 11.8 ~ 12.
    base = dict(gamma1=10.0, gamma1_rad=5.0, alpha3=0.04, delta1=0.0)

    p_nat = TCMTParams(Domega=-2.0, beta=0.0, **base)
    p_cas = TCMTParams(Domega=-2.0, beta=1.0, **base)

    Omega = 0.5  # signal offset within FM linewidth gamma_1 = 10
    s_s = 1e-3   # weak signal

    # Sweep pump power; small enough to stay in low-saturation
    s_p_grid = np.geomspace(0.01, 1.0, 25)

    eta_nat = np.zeros_like(s_p_grid)
    eta_cas = np.zeros_like(s_p_grid)

    for i, s_p in enumerate(s_p_grid):
        idler_n, _ = extract_idler_amplitude(p_nat, s_p, s_s, Omega)
        idler_c, _ = extract_idler_amplitude(p_cas, s_p, s_s, Omega)
        eta_nat[i] = abs(idler_n) ** 2 / s_s ** 2
        eta_cas[i] = abs(idler_c) ** 2 / s_s ** 2

    # Predicted ratio: (1+R)^2 with R = 12 -> 13^2 = 169
    R_predicted = 12
    ratio_predicted = (1 + R_predicted) ** 2

    fig, axs = plt.subplots(1, 2, figsize=(12.0, 4.7))

    axs[0].loglog(s_p_grid ** 2, eta_nat, "o-", color="C2", ms=5,
                  mfc="none", mew=1.2,
                  label=r"native ($\beta = 0$)")
    axs[0].loglog(s_p_grid ** 2, eta_cas, "s-", color="C0", ms=5,
                  mfc="none", mew=1.2,
                  label=r"cascade ($\beta = 1$)")
    # Slope guides -- FWM efficiency ~ |a_1|^4 ~ |s_+|^4 at small signal
    s2_ref = s_p_grid[[0, -1]] ** 2
    eta_ref = eta_nat[0] * (s2_ref / s2_ref[0]) ** 2
    axs[0].loglog(s2_ref, eta_ref, ":", color="0.45", lw=0.9)
    axs[0].text(s2_ref[0] * 2, eta_ref[0] * 3.5, "slope $+2$\n(small-signal)",
                fontsize=9, color="0.4")
    axs[0].set_xlabel(r"$|s_+|^2$  (pump power)")
    axs[0].set_ylabel(r"$\eta_{\rm FWM} = |s_-^{(i)}|^2 / |s_+^{(s)}|^2$")
    axs[0].set_title("(a) FWM conversion efficiency vs pump power")
    axs[0].grid(True, which="both", alpha=0.3)
    axs[0].legend(fontsize=9.5, loc="lower right")

    # Ratio plot
    ratio = eta_cas / eta_nat
    axs[1].semilogx(s_p_grid ** 2, ratio, "o-", color="C0", ms=5,
                    mfc="none", mew=1.2,
                    label=r"$\eta_{\rm cas}/\eta_{\rm nat}$ (TCMT)")
    axs[1].axhline(ratio_predicted, color="C3", lw=1.8, ls="--",
                   label=rf"$(1+R)^2 = ({R_predicted+1})^2 = {ratio_predicted}$")
    axs[1].set_xlabel(r"$|s_+|^2$  (pump power)")
    axs[1].set_ylabel(r"ratio  $\eta_{\rm cas}/\eta_{\rm nat}$")
    axs[1].set_title(r"(b) Architectural ratio matches $(1+R)^2$ at small signal")
    axs[1].grid(True, which="both", alpha=0.3)
    axs[1].legend(fontsize=9.5, loc="upper right")
    axs[1].set_ylim(0, ratio_predicted * 1.4)

    fig.suptitle(
        r"FWM efficiency: cascade outperforms native by $(1+R)^{2}$"
        r" -- $\Omega/\gamma_2 = 0.5$, $\widetilde{\alpha}_3 = +0.04$, "
        r"$\Delta\omega = -2\gamma_2$ (constructive; $R\simeq 12$)",
        fontsize=11.5, y=1.00,
    )
    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "fwm_efficiency"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    # Spot-check the ratio at the smallest power
    print(f"  eta_nat/{s_p_grid[0]**2:.3e} : {eta_nat[0]:.3e}")
    print(f"  eta_cas/{s_p_grid[0]**2:.3e} : {eta_cas[0]:.3e}")
    print(f"  ratio (small signal)         : {ratio[0]:.1f}")
    print(f"  predicted (1+R)^2            : {ratio_predicted}")
    plt.close(fig)


if __name__ == "__main__":
    main()
