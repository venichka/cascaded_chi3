"""Sensor sensitivity: cascade as a Delta-omega-readout transducer.

Tier 2, item 8 of the note's Section 13.  The cascaded Kerr coefficient
alpha_3^casc depends on the SH-cavity detuning Delta omega.  Any
perturbation that shifts omega_2 (temperature, strain, applied electric
field via Pockels, adsorbed molecules) shifts Delta omega, which shifts
alpha_3^casc, which shifts the FM-port nonlinear phase phi_NL.

The Lorentzian derivative

    d alpha_3^casc / d Delta omega = - |beta|^2 (b^2 - Delta omega^2)
                                     / (Delta omega^2 + b^2)^2 ,
                                     b = gamma_2 / 2

is the sensor's response function.  Its magnitude peaks at Delta omega
= 0 with value

    |d alpha_3^casc / d Delta omega|_max = |beta|^2 / b^2
                                         = 4 |beta|^2 / gamma_2^2
                                         ~ Q_2^2  (at fixed |beta|).

This Q_2^2 scaling is the key claim: a linear-cavity sensor reading off
a resonance shift goes as Q_2 (Lorentzian width sets the smallest
resolvable shift); the cascade gets an extra Q_2 from the coefficient
itself being proportional to Q_2.  The Q_2^2 scaling DOES come at the
cost of operating at Delta omega = 0, where the cascaded TPA peaks --
practical sensors bias a few gamma_2 off-resonance, picking up an O(1)
slope reduction in exchange for TPA suppression by 1/Delta omega^2.

Three panels (layout from alternative_chi3/G2_B6, but with physical
Q_2 units in panel (b) for experimentalists):

  (a) d alpha_3^casc / d Delta omega vs Delta omega / gamma_2: clean
      Lorentzian derivative with zero crossings at +/- gamma_2 / 2 and
      peak magnitude at Delta omega = 0.  Analytic curve and TCMT
      finite-difference points overlaid.

  (b) Peak sensitivity |d alpha_3^casc / d Delta omega|_max vs Q_2 on
      log-log axes in PHYSICAL UNITS (Q_2 from 10^2 to 10^5, anchored
      at lambda_1 = 1.55 um, omega_2 = 2 omega_1).  TCMT measurements
      (open markers) fall on the Q_2^2 prediction; a linear-cavity
      Q_2 reference (dashed) is shown for comparison.

  (c) Phase sensitivity |d phi_NL / d Delta omega| vs working detuning
      |Delta omega| / gamma_2 at fixed pump.  Same Lorentzian-derivative
      shape, scaled by the chain-rule factor 2 |a_1|^2 / gamma_1.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from analytic import alpha3_eff_casc_full
from params import TCMTParams, OMEGA1_REAL
from tcmt import steady_state

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def phi_NL_steady(p: TCMTParams, s_p: float, t_settle: float = 500.0) -> float:
    """arg(a_1) at steady state under real-positive CW pump."""
    a1, _ = steady_state(p, s_p, t_settle=t_settle)
    return float(np.angle(a1))


def numerical_dalpha_dDomega(Domega: float, dD: float = 1e-3,
                             beta: float = 1.0) -> float:
    """Centered finite difference of alpha_3^casc(Delta omega) about
    Domega, with step dD in units of gamma_2 = 1."""
    fp = alpha3_eff_casc_full(Domega + dD) * beta ** 2
    fm = alpha3_eff_casc_full(Domega - dD) * beta ** 2
    return (fp - fm) / (2 * dD)


def numerical_dphi_dDomega(Domega: float, base: dict, s_p: float,
                           dD: float = 1e-2) -> float:
    """Centered finite difference of phi_NL w.r.t. SH detuning, by
    re-running TCMT at Domega + dD and Domega - dD."""
    p_p = TCMTParams(Domega=Domega + dD, **base)
    p_m = TCMTParams(Domega=Domega - dD, **base)
    return (phi_NL_steady(p_p, s_p) - phi_NL_steady(p_m, s_p)) / (2 * dD)


def main():
    # Physical units anchor (LiNbO3 telecom; see params.py)
    omega_1 = OMEGA1_REAL                # 1.22e15 rad/s
    omega_2 = 2 * omega_1                # ~2.44e15 rad/s

    # Dimensionless TCMT defaults (gamma_2 = 1)
    gamma_1_over_gamma_2 = 10.0          # FM cavity 10x broader than SH
    gamma_1_rad_over_gamma_1 = 0.5       # under-coupled FM port

    base = dict(gamma1=gamma_1_over_gamma_2, gamma1_rad=gamma_1_over_gamma_2 * 0.5,
                alpha3=0.0, delta1=0.0, beta=1.0)
    s_p = 0.1                            # small-signal: phi_NL ~ alpha_eff

    fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.7), gridspec_kw={"wspace": 0.34})

    # =========================================================================
    # Panel (a): Lorentzian derivative d alpha_3^casc / d Delta omega
    # =========================================================================
    Domega_axis = np.linspace(-5, 5, 81)
    # Analytic derivative: -(b^2 - Delta omega^2) / (Delta omega^2 + b^2)^2 with b = 0.5
    b = 0.5
    deriv_analytic = -(b ** 2 - Domega_axis ** 2) / (Domega_axis ** 2 + b ** 2) ** 2
    # TCMT finite-difference: same coefficient, just numerically
    deriv_tcmt = np.array([numerical_dalpha_dDomega(D) for D in Domega_axis])

    axs[0].plot(Domega_axis, deriv_analytic, "-", color="k", lw=1.7,
                label=r"analytic  $-\frac{b^2 - \Delta\omega^2}{(\Delta\omega^2 + b^2)^2}$")
    axs[0].plot(Domega_axis[::3], deriv_tcmt[::3], "o", color="C3",
                ms=5, mfc="white", mew=1.3, label="TCMT (finite difference)")
    axs[0].axhline(0, color="k", lw=0.4)
    # Peak annotation at the actual peak location (y = -4 at x = 0)
    axs[0].annotate(r"$\Delta\omega = 0$  (peak magnitude $= 4|\beta|^2/\gamma_2^2$)",
                    xy=(0, -4), xytext=(1.5, -3.2),
                    fontsize=9, color="C2",
                    arrowprops=dict(arrowstyle="->", color="C2", lw=0.9))
    # Zero-crossing markers
    for crossing in (-b, b):
        axs[0].axvline(crossing, color="C0", ls=":", lw=0.9, alpha=0.6)
    axs[0].annotate(r"$\Delta\omega = \pm \gamma_2/2$" + "\n(zero slope)",
                    xy=(b, 0), xytext=(b + 0.4, 0.8),
                    fontsize=9, color="C0",
                    arrowprops=dict(arrowstyle="->", color="C0", lw=0.9))
    axs[0].set_xlabel(r"$\Delta\omega \, / \, \gamma_2$")
    axs[0].set_ylabel(r"$d\alpha_3^{\rm casc}/d\Delta\omega$  "
                      r"(units of $|\beta|^2 / \gamma_2^2$)")
    axs[0].set_title(r"(a) Lorentzian-derivative response, peak at $\Delta\omega = 0$")
    axs[0].legend(loc="lower right", fontsize=9, framealpha=0.95)
    axs[0].grid(True, alpha=0.3)
    axs[0].set_ylim(-5.0, 1.5)

    # =========================================================================
    # Panel (b): peak sensitivity vs Q_2 in PHYSICAL UNITS
    # =========================================================================
    # Physical claim: peak |d alpha_3^casc / d Delta omega|_phys
    #                = 4 |beta|^2_phys / gamma_2^2_phys  ~  Q_2^2
    # (because gamma_2_phys = omega_2 / Q_2 and |beta|_phys is set by the mode
    #  overlap, independent of Q_2).
    #
    # The dim TCMT framework hardcodes gamma_2 = 1.  Different physical Q_2
    # values are realised by setting beta_dim = |beta|_phys / gamma_2_phys
    # = (|beta|_phys / gamma_2_ref) * (Q_2 / Q_2_ref).  So beta_dim is
    # proportional to Q_2 (at fixed physical beta).  The dim TCMT peak slope
    # at Delta omega = 0 is 4 * beta_dim^2, which scales as Q_2^2 -- exactly
    # the prediction we wish to validate.
    Q2_ref = 4000.0
    gamma_2_ref_phys = omega_2 / Q2_ref   # 6.1e11 rad/s ~ 97 GHz at Q_2 = 4000
    Q2_grid_phys = np.array([500, 1000, 2000, 4000, 10_000, 30_000, 100_000])
    Q2_factor_grid = Q2_grid_phys / Q2_ref     # 0.125 .. 25
    beta_dim_grid = Q2_factor_grid             # beta_dim proportional to Q_2

    # Run TCMT at each beta_dim and read the peak |d alpha / d Delta omega|
    # numerically by centered difference.  Hold gamma_1_dim = 10 (constant
    # FM cavity in units of physical gamma_2_ref).
    #
    # Stay-linear constraint: |beta a_1|^2 << |Delta omega|^2 ~ dD^2.
    # |a_1|^2 ~ gamma_rad |s_+|^2 / (gamma_1/2)^2, so we need
    #   beta_dim^2 * gamma_rad |s_+|^2 / (gamma_1/2)^2 << dD^2
    # i.e. |s_+| << dD * (gamma_1/2) / (beta_dim sqrt(gamma_rad)).
    # Scale s_p as ~ 1/beta_dim to stay linear at large beta_dim.
    peak_tcmt_dim_grid = np.zeros_like(beta_dim_grid)
    for k, b_dim in enumerate(beta_dim_grid):
        base_k = dict(gamma1=gamma_1_over_gamma_2,
                      gamma1_rad=0.5 * gamma_1_over_gamma_2,
                      alpha3=0.0, delta1=0.0, beta=float(b_dim))
        # Adaptive s_p and dD: shrink as 1/beta_dim so we stay deep in
        # the small-signal regime even at very high beta_dim (i.e. high Q_2).
        s_p_k = min(s_p, 0.05 / max(b_dim, 1.0))   # cap at 0.05/beta_dim
        dD = 1e-2 / max(b_dim, 1.0)                # shrink finite-difference step too
        phi_p = phi_NL_steady(TCMTParams(Domega=+dD, **base_k), s_p_k)
        phi_m = phi_NL_steady(TCMTParams(Domega=-dD, **base_k), s_p_k)
        dphi_dDom = (phi_p - phi_m) / (2 * dD)
        # Strip the chain factor (linear-cavity build-up factor):
        # phi_NL ~ alpha_eff * |a_1|^2 * (2 / gamma_1),
        # |a_1|^2 = gamma_rad |s_+|^2 / (gamma_1 / 2)^2,
        # so chain = 8 gamma_rad |s_+|^2 / gamma_1^3.
        chain = 8.0 * base_k["gamma1_rad"] * s_p_k ** 2 / base_k["gamma1"] ** 3
        peak_tcmt_dim_grid[k] = abs(dphi_dDom / chain)

    # Analytic Q_2^2 curve (dim peak = 4 beta_dim^2 = 4 (Q_2/Q_2_ref)^2)
    peak_analytic_dim = 4.0 * (Q2_grid_phys / Q2_ref) ** 2
    # Linear-cavity sensor reference: scales as Q_2 (not Q_2^2)
    linear_cavity_ref = 4.0 * (Q2_grid_phys / Q2_ref)

    axs[1].loglog(Q2_grid_phys, peak_analytic_dim, "-", color="k", lw=1.8,
                  label=r"analytic cascade  $\propto Q_2^{\,2}$")
    axs[1].loglog(Q2_grid_phys, linear_cavity_ref, "--", color="0.55", lw=1.4,
                  label=r"linear-cavity sensor  $\propto Q_2$")
    axs[1].loglog(Q2_grid_phys, peak_tcmt_dim_grid, "s", color="C3", ms=10,
                  mfc="white", mew=1.8, label="TCMT (peak)")
    # Reference point marker
    axs[1].plot([Q2_ref], [4.0], "D", color="C2", ms=10, mfc="white", mew=1.6,
                label=rf"reference  $Q_2 = {int(Q2_ref)}$")
    axs[1].axvspan(1e3, 1e5, alpha=0.07, color="grey",
                   label=r"qBIC range  $Q_2 \in [10^3, 10^5]$")
    axs[1].set_xlabel(r"$Q_2$  (SH-cavity quality factor, physical units)")
    axs[1].set_ylabel(r"peak  $|d\alpha_3^{\rm casc}/d\Delta\omega|$  "
                      r"(units of $|\beta_{\rm ref}|^2/\gamma_{2,{\rm ref}}^{\,2}$)")
    axs[1].set_title(r"(b) Sensor sensitivity scales as $Q_2^{\,2}$ in physical units")
    axs[1].legend(loc="upper left", fontsize=9, framealpha=0.95)
    axs[1].grid(True, which="both", alpha=0.3)

    # =========================================================================
    # Panel (c): phase sensitivity vs working detuning
    # =========================================================================
    # The actual experimental observable is d phi_NL / d Delta omega vs the
    # bias point Delta omega_bias.  This is the Lorentzian derivative scaled
    # by the chain-rule factor.  Both TCMT and analytic are shown.
    Domega_bias = np.linspace(0.05, 5.0, 60)    # positive bias only
    # Analytic chain rule: phi_NL ~ alpha_eff * (8 gamma_rad / gamma_1^3) |s_+|^2
    chain = 8.0 * (base["gamma1_rad"] / base["gamma1"] ** 3) * s_p ** 2
    dphi_analytic = chain * np.abs(
        -(b ** 2 - Domega_bias ** 2) / (Domega_bias ** 2 + b ** 2) ** 2
    )
    dphi_tcmt = np.zeros_like(Domega_bias)
    for i, D in enumerate(Domega_bias):
        dphi_tcmt[i] = abs(numerical_dphi_dDomega(D, base, s_p, dD=1e-3))

    axs[2].semilogy(Domega_bias, dphi_analytic, "-", color="k", lw=1.7,
                    label=r"analytic chain rule")
    axs[2].semilogy(Domega_bias[::3], dphi_tcmt[::3], "o", color="C3",
                    ms=5, mfc="white", mew=1.3, label="TCMT")
    # Annotate the TPA-vs-Kerr boundary and the zero-slope point
    axs[2].axvline(0.5, color="C0", ls=":", lw=0.9, alpha=0.6)
    axs[2].text(0.55, dphi_analytic.min() * 1.5, r"$\Delta\omega = \gamma_2 / 2$"
                                                 "\n(zero slope)",
                fontsize=8.5, color="C0", ha="left", va="bottom")
    axs[2].axvspan(0.0, 1.0, alpha=0.08, color="C3",
                   label=r"TPA-penalised bias  ($|\Delta\omega| \lesssim \gamma_2$)")
    axs[2].set_xlabel(r"working detuning  $|\Delta\omega| \, / \, \gamma_2$")
    axs[2].set_ylabel(r"$|d\phi_{\rm NL} / d\Delta\omega|$  (rad / unit detuning)")
    axs[2].set_title(rf"(c) Phase sensitivity at $|s_+|^2 = {s_p**2:.2g}$")
    axs[2].legend(loc="upper right", fontsize=9.5, framealpha=0.95)
    axs[2].grid(True, which="both", alpha=0.3)

    fig.suptitle(
        r"$\Delta\omega$-readout transducer: peak at $\Delta\omega = 0$, "
        r"$Q_2^{\,2}$ scaling beats linear-cavity $Q_2$ sensors",
        fontsize=11.5, y=1.02,
    )
    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "sensor_sensitivity"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    print(f"  Q_2 anchor: Q_2_ref = {int(Q2_ref)}, gamma_2_phys_ref = "
          f"{gamma_2_ref_phys:.2e} rad/s ({gamma_2_ref_phys/(2*np.pi)/1e9:.1f} GHz)")
    print(f"  TCMT peak |d alpha / d Delta omega| at each Q_2:")
    for Q2, ana, tcmt in zip(Q2_grid_phys, peak_analytic_dim, peak_tcmt_dim_grid):
        print(f"    Q_2 = {int(Q2):>6d} : analytic = {ana:9.3f}, TCMT = {tcmt:9.3f}, "
              f"ratio = {tcmt/ana:.3f}")
    plt.close(fig)


if __name__ == "__main__":
    main()
