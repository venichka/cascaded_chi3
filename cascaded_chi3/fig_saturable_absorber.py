"""Cascade as a sign-tunable saturable absorber.

Tier 2, item 5 of the note's Section 13. The cascade simultaneously
generates a Kerr-like (dispersive) and a TPA-like (absorptive)
contribution on the FM, controlled by a single SH-cavity detuning:

    alpha_3^casc(Delta omega) = - |beta|^2 * Delta omega
                                  / (Delta omega^2 + (gamma_2 / 2)^2)
                                                                (dispersive, Kerr)
    gamma_3^casc(Delta omega) =   |beta|^2 * gamma_2 / 2
                                  / (Delta omega^2 + (gamma_2 / 2)^2)
                                                                (absorptive, TPA-like)

The TPA Lorentzian PEAKS at Delta omega = 0 with magnitude
|beta|^2 / (gamma_2 / 2) = 2 |beta|^2 / gamma_2 (no Kerr at all there);
the Kerr coefficient is dispersive and peaks at Delta omega = +/- gamma_2 / 2
with magnitude |beta|^2 / gamma_2. The two coincide at |Delta omega| =
gamma_2 / 2.  Sweeping Delta omega across these scales tunes the
nonlinear response continuously from pure-SA (Delta omega = 0) to
pure-Kerr (|Delta omega| >> gamma_2), without any change in material.

Comparison with SESAMs (semiconductor saturable absorber mirrors): the
cascaded SA has recovery time ~ 1/gamma_2 (ps at Q_2 = 4000) and is
all-optical with no excited-state physics or thermal drift; SESAMs have
a fixed, material-determined recovery time.  The trade-off is that
cascaded SA requires a chi^(2)-active substrate and double-resonance
fabrication.

Three panels (layout from alternative_chi3/G2_B4):

  (a) Cascade Kerr and TPA coefficients vs Delta omega / gamma_2.
      Both on the same axes, with the TPA peak at zero and the Kerr
      S-shape with zero crossings at the same point.

  (b) Intensity-dependent FM-port transmission |s_- / s_+|^2 vs input
      |s_+|^2 for three Delta omega values: 0 (peak SA), gamma_2
      (intermediate), 3 gamma_2 (off-resonance, TPA suppressed). A
      linear-cavity reference (no NL) is overlaid.

  (c) Loss / Kerr ratio gamma_3^casc / |alpha_3^casc| vs |Delta omega|
      / gamma_2. The ratio crosses 1 at |Delta omega| = gamma_2 / 2,
      marking the boundary between SA-dominated (TPA > Kerr) and
      Kerr-dominated operation. The TPA-dominated and Kerr-dominated
      regions are shaded for visual orientation.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from analytic import alpha3_eff_casc_full, gamma_2pa_eff
from params import TCMTParams
from tcmt import integrate, steady_state

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def transmission_cw(p: TCMTParams, s_p: float, t_settle: float = 800.0):
    """Steady-state FM-port transmission |s_- / s_+|^2 under CW pump
    (assumes a stable fixed point)."""
    a1, _ = steady_state(p, s_p, t_settle=t_settle)
    s_minus = -s_p + np.sqrt(p.gamma1_rad) * a1
    return float(abs(s_minus) ** 2 / s_p ** 2)


def transmission_attractor(p: TCMTParams, s_p: float,
                           t_settle: float = 250.0, t_avg: float = 250.0,
                           n_samples: int = 1200):
    """Integrate full TCMT long enough to reach the attractor (fixed point
    OR parametric limit cycle), then time-average |s_-/s_+|^2 over the
    last t_avg interval.  Works regardless of whether the trivial fixed
    point is stable -- limit cycles give a well-defined time average."""
    t_eval = np.linspace(t_settle, t_settle + t_avg, n_samples)
    sol = integrate(p, s_p, (0.0, t_settle + t_avg), t_eval=t_eval,
                    max_step=0.5)
    a1 = sol.a[0]
    s_minus = -s_p + np.sqrt(p.gamma1_rad) * a1
    T_t = np.abs(s_minus) ** 2 / s_p ** 2
    return float(np.mean(T_t)), float(np.std(T_t))


def main():
    # FM cavity: over-coupled (gamma_1_rad = gamma_1, no intrinsic loss).
    # This makes the LINEAR-cavity transmission T_lin = 1 on resonance,
    # matching the conventional SA convention used in alternative_chi3/G2_B4:
    # adding cascade TPA pulls T DOWN from 1 toward critical coupling,
    # which is the conventional "intensity-dependent loss" picture.
    gamma_1 = 1.0
    gamma_1_rad = 1.0       # over-coupled: no intrinsic FM loss
    beta = 1.0

    fig, axs = plt.subplots(1, 3, figsize=(15.5, 5.1), gridspec_kw={"wspace": 0.34})

    # =========================================================================
    # Panel (a): Cascade coefficients vs SH detuning
    # =========================================================================
    Domega_axis = np.linspace(-10, 10, 1000)
    # Dimensionless |β|=1, γ_2=1 ⇒ b = 0.5
    alpha_kerr = alpha3_eff_casc_full(Domega_axis)    # dispersive, S-shaped
    gamma_tpa  = gamma_2pa_eff(Domega_axis)            # absorptive, Lorentzian

    axs[0].plot(Domega_axis, alpha_kerr, "-", color="C0", lw=2.0,
                label=r"$\alpha_3^{\rm casc}$ (cascaded Kerr, dispersive)")
    axs[0].plot(Domega_axis, gamma_tpa, "-", color="C3", lw=2.0,
                label=r"$\gamma_3^{\rm casc}$ (cascaded TPA, absorptive)")
    axs[0].fill_between(Domega_axis, 0, gamma_tpa, alpha=0.10, color="C3")
    axs[0].axhline(0, color="k", lw=0.5)
    axs[0].axvline(0, color="0.5", ls="-.", lw=0.9)
    # Peak marker for the TPA Lorentzian (inside panel, not above)
    peak_tpa = 2.0  # = 2 |β|^2 / γ_2 with β=γ_2=1
    axs[0].axhline(peak_tpa, color="C3", ls=":", lw=0.9, alpha=0.6)
    axs[0].annotate(rf"TPA peak $= 2|\beta|^2/\gamma_2$",
                    xy=(0, peak_tpa), xytext=(2.5, peak_tpa - 0.15),
                    fontsize=9, color="C3", ha="left", va="top",
                    arrowprops=dict(arrowstyle="->", color="C3", lw=0.8))
    axs[0].set_xlabel(r"$\Delta\omega / \gamma_2$")
    axs[0].set_ylabel(r"coefficient  (units of $|\beta|^2 / \gamma_2$)")
    axs[0].set_title(r"(a) Cascade coefficients vs SH detuning")
    axs[0].legend(loc="lower right", fontsize=9.5, framealpha=0.95)
    axs[0].grid(True, alpha=0.3)
    axs[0].set_ylim(-1.2, 2.5)

    # =========================================================================
    # Panel (b): intensity-dependent transmission, three Delta omega values
    # =========================================================================
    # Full TCMT (no adiabatic restriction) integrated to attractor; T is the
    # time-average of |s_-/s_+|^2 over the last 250/gamma_2 of the integration.
    # This is well-defined whether the trivial fixed point is stable (regime A)
    # or has lost stability to parametric oscillation (regime B); the latter
    # gives a limit-cycle time average rather than a fixed-point value.
    # The parametric-instability boundary |beta a_1| = |Delta omega| is shown
    # as a shaded region on each curve.
    dw_values = [0.5, 1.0, 3.0]      # avoid Delta omega = 0 exactly (singular)
    colors_b  = ["C3", "C1", "C0"]
    labels_b  = [
        r"$\Delta\omega = 0.5\,\gamma_2$  (near peak TPA)",
        r"$\Delta\omega = \gamma_2$  (intermediate)",
        r"$\Delta\omega = 3\,\gamma_2$  (off-resonance, TPA suppressed)",
    ]

    # Pump-power grid (log-spaced); cover linear -> NL -> deep saturation.
    s_p_grid = np.geomspace(1e-2, 20.0, 22)
    T_curves = {}
    n1_curves = {}
    for dw, color, label in zip(dw_values, colors_b, labels_b):
        p = TCMTParams(Domega=dw, gamma1=gamma_1, gamma1_rad=gamma_1_rad,
                       alpha3=0.0, delta1=0.0, beta=beta)
        T_arr = np.zeros_like(s_p_grid)
        n1_arr = np.zeros_like(s_p_grid)
        for i, sp in enumerate(s_p_grid):
            T_mean, _ = transmission_attractor(p, sp)
            T_arr[i] = T_mean
            # Also need |a_1| at the attractor to mark parametric instability:
            a1, _ = steady_state(p, sp, t_settle=400.0)
            n1_arr[i] = abs(a1) ** 2
        T_curves[dw] = T_arr
        n1_curves[dw] = n1_arr
        # |beta a_1| = |Delta omega| boundary: split curve into regime-A
        # (stable, solid) and regime-B (limit-cycle / unstable, open markers).
        unstable = (beta * np.sqrt(n1_arr)) > dw
        # Plot the safe (regime A) range as solid line + filled markers:
        axs[1].plot(s_p_grid[~unstable] ** 2, T_arr[~unstable], "o-",
                    color=color, ms=5, mfc=color, mew=0,
                    lw=1.5, label=label)
        # Plot the unstable (regime B) range as dashed line + open markers:
        if unstable.any():
            axs[1].plot(s_p_grid[unstable] ** 2, T_arr[unstable], "s--",
                        color=color, ms=5, mfc="white", mew=1.2,
                        lw=1.0, alpha=0.7)

    # Linear-cavity reference (beta=0): T = 1 on resonance for over-coupled cavity.
    p_lin = TCMTParams(Domega=0.0, gamma1=gamma_1, gamma1_rad=gamma_1_rad,
                       alpha3=0.0, delta1=0.0, beta=0.0)
    T_lin = transmission_cw(p_lin, 1.0)
    axs[1].axhline(T_lin, color="0.4", ls="--", lw=1.2,
                   label=rf"linear cavity (no NL): $T = {T_lin:.2f}$")

    axs[1].set_xscale("log")
    axs[1].set_xlabel(r"input power  $|s_+|^2$")
    axs[1].set_ylabel(r"FM-port transmission  $\langle|s_- / s_+|^2\rangle_t$")
    axs[1].set_title(r"(b) Intensity-dependent transmission (full TCMT)")
    axs[1].legend(loc="lower left", fontsize=8.5, framealpha=0.95)
    axs[1].grid(True, which="both", alpha=0.3)
    axs[1].set_ylim(-0.02, 1.1)
    # Marker-convention note in upper-middle of the panel
    axs[1].text(0.50, 0.55,
                "filled: stable fixed point\n"
                r"open: past $|\beta a_1|\!=\!|\Delta\omega|$,"
                "\n   limit-cycle time-average",
                transform=axs[1].transAxes,
                fontsize=8.5, color="0.20", ha="center", va="top",
                bbox=dict(facecolor="white", edgecolor="0.6",
                          boxstyle="round,pad=0.4"))

    # =========================================================================
    # Panel (c): loss-to-Kerr ratio vs |Delta omega|
    # =========================================================================
    dw_pos = np.linspace(0.01, 8.0, 500)
    ratio  = gamma_2pa_eff(dw_pos) / np.abs(alpha3_eff_casc_full(dw_pos))
    # The ratio simplifies to (gamma_2/2) / |Delta omega| -- monotonic.
    crossing = 0.5    # gamma_2 / 2 in units of gamma_2

    axs[2].plot(dw_pos, ratio, "-", color="0.15", lw=2.0)
    axs[2].axhline(1.0, color="0.5", ls="-.", lw=0.9)
    axs[2].axvline(crossing, color="0.5", ls="-.", lw=0.9)
    # Shading: TPA-dominated for |Δω| < γ_2/2 (ratio > 1); Kerr-dominated otherwise.
    axs[2].fill_between(dw_pos, 0, ratio, where=(dw_pos < crossing),
                        color="C3", alpha=0.18,
                        label=r"TPA-dominated  ($|\Delta\omega| < \gamma_2 / 2$)")
    axs[2].fill_between(dw_pos, 0, ratio, where=(dw_pos >= crossing),
                        color="C0", alpha=0.18,
                        label=r"Kerr-dominated  ($|\Delta\omega| > \gamma_2 / 2$)")
    axs[2].annotate(r"crossing at $|\Delta\omega| = \gamma_2 / 2$",
                    xy=(crossing, 1.0),
                    xytext=(crossing + 0.4, 2.2),
                    fontsize=9.5, color="0.25",
                    arrowprops=dict(arrowstyle="->", color="0.3", lw=0.8))
    axs[2].set_xlabel(r"$|\Delta\omega| / \gamma_2$")
    axs[2].set_ylabel(r"$\gamma_3^{\rm casc} \,/\, |\alpha_3^{\rm casc}|$")
    axs[2].set_title(r"(c) SA character: loss / Kerr ratio")
    axs[2].legend(loc="upper right", fontsize=9.5, framealpha=0.95)
    axs[2].grid(True, alpha=0.3)
    axs[2].set_ylim(0, 5)
    axs[2].set_xlim(0, 4)

    fig.suptitle(
        r"Cascaded saturable absorber: $\Delta\omega$ tunes between absorptive (TPA) "
        r"and dispersive (Kerr) cascade response",
        fontsize=12, y=1.00,
    )
    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "saturable_absorber"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    # Sanity checks
    print(f"  alpha_kerr peak at Domega = +-gamma_2/2 = +-0.5: "
          f"|alpha| = {abs(alpha3_eff_casc_full(0.5)):.3f} = |beta|^2/gamma_2  "
          f"(expected 1.0)")
    print(f"  gamma_tpa peak at Domega = 0: "
          f"gamma = {gamma_2pa_eff(0.0):.3f} = 2|beta|^2/gamma_2  "
          f"(expected 2.0)")
    print(f"  ratio at |Domega|/gamma_2 = 0.5: "
          f"{gamma_2pa_eff(0.5)/abs(alpha3_eff_casc_full(0.5)):.3f}  "
          f"(expected 1.0)")
    plt.close(fig)


if __name__ == "__main__":
    main()
