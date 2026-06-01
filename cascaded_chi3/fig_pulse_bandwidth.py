"""Pulsed pump: the SH cannot follow if the pulse is shorter than 1/gamma2.

Drive with a Gaussian pump pulse:
    s+(t) = s0 * exp(-(t - t0)^2 / (2 sigma_t^2))

at the FM carrier (delta1 = 0, Domega fixed). Vary sigma_t. For
sigma_t * gamma2 >> 1 the SH tracks |a1|^2 in time and the cascaded
nonlinearity acts as an instantaneous Kerr. For sigma_t * gamma2 ~ 1
the SH cannot build up; the adiabatic prediction over-estimates |a2|.

Two panels:
  (a) sample time traces |a2(t)|^2 (full TCMT) overlaid with the
      adiabatic prediction |a1(t)|^4 / (1/4 + Domega^2), for one short
      and one long pulse.
  (b) instantaneous SH-tracking ratio sampled at the moment of FM peak:
          R = (|a2|^2 / |a1|^4)|_{t_peak} * (1/4 + Domega^2)
      where t_peak = argmax |a1|^2. This is the cleanest indicator: at
      the instant the FM is brightest, how close is the SH to its
      instantaneous adiabatic value? R -> 1 for long pulses, R -> 0 for
      pulses shorter than 1/gamma2 (SH lags behind). Crossover at
      sigma_t * gamma2 ~ 1, the time-domain face of Section 6.
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


def gaussian_pump(s0, t0, sigma_t):
    def s(t):
        return s0 * np.exp(-((t - t0) ** 2) / (2 * sigma_t ** 2))
    return s


def run_pulse(p, s0, sigma_t, factor=8.0, n_samples=8000):
    t0 = factor * sigma_t + 5.0 / p.gamma1
    t_end = 2 * t0
    t_eval = np.linspace(0, t_end, n_samples)
    sp = gaussian_pump(s0, t0, sigma_t)
    sol = integrate(p, sp, (0, t_end), t_eval=t_eval,
                    max_step=min(sigma_t / 5, 0.5))
    return sol, t0


def main():
    p = TCMTParams(Domega=2.0, gamma1=1.0, gamma1_rad=0.5, alpha3=0.0,
                   delta1=0.0)

    norm_adiab = 1.0 / (0.25 + p.Domega ** 2)  # |a2|^2 / |a1|^4 in adiabatic

    # 2x2 layout: top row = short and long pulse traces; bottom = duration
    # sweep across the full width. Aspect ratio is closer to square so
    # the per-panel detail is legible in print/slides.
    fig = plt.figure(figsize=(11.5, 7.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.1], hspace=0.35,
                          wspace=0.25)
    ax_short = fig.add_subplot(gs[0, 0])
    ax_long = fig.add_subplot(gs[0, 1])
    ax_sweep = fig.add_subplot(gs[1, :])

    # --------- panels (a, b): time traces, real magnitudes ----------
    s0 = 0.5
    for ax, sigma_t, label_pulse in [
        (ax_short, 0.3, r"(a) short pulse, $\sigma_t \gamma_2 = 0.3$"),
        (ax_long, 8.0, r"(b) long pulse, $\sigma_t \gamma_2 = 8$"),
    ]:
        sol, t0 = run_pulse(p, s0, sigma_t)
        a1, a2 = sol.a
        n2 = np.abs(a2) ** 2
        n2_adiab = np.abs(a1) ** 4 * norm_adiab
        sp_t = s0 * np.exp(-((sol.t - t0) ** 2) / (2 * sigma_t ** 2))

        # x-axis in pulse-width units so both panels look comparable
        x = (sol.t - t0) / sigma_t
        ax.plot(x, n2_adiab, ":", color="C0", lw=1.4,
                label=r"adiabatic $|a_1|^4/(1/4+\Delta\omega^2)$")
        ax.plot(x, n2, "-", color="C3", lw=1.6, label=r"full TCMT $|a_2|^2$")

        # pump shape on a twin axis for shape comparison
        axp = ax.twinx()
        axp.plot(x, np.abs(sp_t) ** 2, "-", color="grey", lw=0.9, alpha=0.5)
        axp.set_ylabel(r"$|s_+|^2$", color="grey", fontsize=9)
        axp.tick_params(axis="y", labelcolor="grey", labelsize=8)
        axp.set_ylim(0, 1.05 * (np.abs(sp_t) ** 2).max())

        ax.set_xlabel(r"$(t-t_0)/\sigma_t$")
        ax.set_ylabel(r"$|a_2|^2$")
        ax.set_title(label_pulse)
        ax.set_xlim(-4, 4)
        ax.legend(loc="upper left", fontsize=9.5)
        ax.grid(True, alpha=0.3)

    # --------- panel (c): tracking quality vs pulse duration ------------
    sigma_grid = np.geomspace(0.05, 30.0, 30)
    s0_b = 0.3  # well below saturation
    R = np.empty_like(sigma_grid)

    for i, sigma_t in enumerate(sigma_grid):
        sol, t0 = run_pulse(p, s0_b, sigma_t, factor=8.0, n_samples=8000)
        a1, a2 = sol.a
        n1 = np.abs(a1) ** 2
        idx = int(np.argmax(n1))
        n1_peak = n1[idx]
        n2_at_peak = abs(a2[idx]) ** 2
        # adiabatic prediction at that instant: |a1|^4 / (1/4 + Domega^2)
        n2_adiab_at_peak = n1_peak ** 2 * norm_adiab
        R[i] = n2_at_peak / n2_adiab_at_peak if n2_adiab_at_peak > 0 else 0.0

    ax_sweep.semilogx(sigma_grid, R, "o-", color="C0", ms=4.5,
                      label=r"$|a_2|^2 / |a_{2,{\rm adiab}}|^2$ at FM peak")
    ax_sweep.axhline(1.0, color="grey", lw=1.2, ls="--", label="adiabatic limit")
    ax_sweep.axvline(1.0, color="C3", lw=0.8, ls=":")
    ax_sweep.text(1.05, 0.06, r"$\sigma_t \gamma_2 = 1$",
                  color="C3", fontsize=10)
    ax_sweep.set_xlabel(r"$\sigma_t \, \gamma_2$  (pulse duration in SH lifetimes)")
    ax_sweep.set_ylabel(r"peak SH energy / adiabatic prediction")
    ax_sweep.set_title("(c) Bandwidth limit (time-domain sum rule)")
    ax_sweep.set_ylim(0, 1.4)
    ax_sweep.grid(True, which="both", alpha=0.3)
    ax_sweep.legend(loc="lower right")

    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "pulse_bandwidth"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
