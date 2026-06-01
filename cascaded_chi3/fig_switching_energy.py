"""All-optical switching: control-pulse energy required to flip the
bistable element from lower to upper branch.

Tier 1, item 3 of the note's Section 13. The cascade has a 1/(1+R) =
1/13 lower bistability threshold (for LiNbO_3 at the conservative
operating point Delta omega = 2 gamma_2, R ~ 12) than the native-Kerr
architecture, and the same factor applies to the minimum energy of
the control pulse required to switch between branches.

Setup. Bias the system on the lower branch of the bistable region
with a CW pump just below the snap-up threshold. Apply a brief
Gaussian "set" pulse on top of the CW: if its energy exceeds a
threshold, the system snaps to the upper branch and stays there.
Sweep the set-pulse amplitude and find the minimum.

Both architectures use a CW bias such that the system sits 95% of
the way to the snap-up cusp (in input-intensity terms). The native
architecture (alpha_3 = R x cascade_alpha_3_native; we boost alpha_3
to bring its bistability cusp into reach) has its cusp at much higher
absolute power. We report ABSOLUTE switching energies (E ~ s_p^2 *
duration) -- the cascade is 1/(1+R) of the native, in line with the
analytic prediction.

Two panels:
  (a) Time traces of the switching event for cascade and native at
      their respective minimum control-pulse energies. Same shape;
      different absolute amplitudes.
  (b) Switching-energy threshold vs alpha_3 (sweep R indirectly via
      varying the native Kerr coefficient at fixed cascade design).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from analytic import cw_kerr_only_intracavity
from params import TCMTParams
from tcmt import integrate

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def find_snap_up(p: TCMTParams, s_max: float = 6.0, n: int = 800):
    """Snap-up power of the analytic Kerr+TPA cubic (Kerr-only model)."""
    s2_axis = np.linspace(0.001, s_max, n)
    n_roots = np.array([len(cw_kerr_only_intracavity(p, np.sqrt(s)))
                        for s in s2_axis])
    bist = s2_axis[n_roots >= 3]
    if bist.size == 0:
        return None
    return float(bist.max())


def switch_with_pulse(p: TCMTParams, s_bias: float, s_pulse_peak: float,
                      sigma_t: float, t0: float, t_total: float):
    """Drive with CW bias + Gaussian pulse; return time trace of |a_1|^2."""
    def s_plus(t):
        return s_bias + s_pulse_peak * np.exp(-((t - t0) / sigma_t) ** 2)
    t_eval = np.linspace(0, t_total, 4000)
    sol = integrate(p, s_plus, (0, t_total), t_eval=t_eval, max_step=0.5)
    return sol.t, np.abs(sol.a[0]) ** 2


def find_switching_threshold(p: TCMTParams, s_bias: float, sigma_t: float,
                             pulse_grid: np.ndarray) -> float | None:
    """Find minimum pulse peak amplitude that snaps the system to the
    upper branch (and STAYS there long after the pulse passes)."""
    t0 = 5 * sigma_t + 30.0 / p.gamma1
    t_total = t0 + 250.0 / p.gamma1
    for s_peak in pulse_grid:
        t, n1 = switch_with_pulse(p, s_bias, s_peak, sigma_t, t0, t_total)
        # Check |a_1|^2 late: did it stay on upper branch?
        n1_final = np.mean(n1[int(0.92 * len(n1)):])
        # The lower-branch and upper-branch values: use the analytic cubic
        roots = cw_kerr_only_intracavity(p, s_bias)
        if roots.size >= 3:
            n1_low, n1_mid, n1_up = roots[0], roots[1], roots[-1]
            mid = 0.5 * (n1_mid + n1_up)
            if n1_final > mid:
                return float(s_peak)
        else:
            # not bistable -- treat any large n1 as switched
            if n1_final > 2 * np.mean(n1[:200]):
                return float(s_peak)
    return None


def main():
    # SAME host material in both architectures (alpha_3 = bare LiNbO_3
    # Kerr in dimensionless units). The note's "1/(1+R) switching-energy
    # ratio" applies to this comparison: same alpha_3, same FM cavity,
    # difference is whether the SH cavity is engineered (beta=1) or not
    # (beta=0). Use alpha_3 sign chosen so the effective Kerr is
    # self-defocusing in both schemes -- this is the geometry where
    # bistability emerges with a blue-detuned pump.
    # Choose Domega large enough that the cascaded TPA is negligible
    # compared with the Kerr (gamma_2pa/|alpha_eff| = (gamma_2/2)/|Domega|
    # is small at Domega = 10*gamma_2). This is the regime where the
    # note's "1/(1+R)" prediction is exact.
    alpha3_native = -0.02   # |alpha_3| representing LiNbO_3 bulk Kerr
    Domega = 10.0
    delta1 = 0.8
    gamma1, gamma1_rad = 0.4, 0.2

    p_cas = TCMTParams(Domega=Domega, gamma1=gamma1, gamma1_rad=gamma1_rad,
                       alpha3=alpha3_native, delta1=delta1, beta=1.0)
    snap_cas = find_snap_up(p_cas, s_max=15.0, n=1500)
    s_bias_cas = np.sqrt(0.75 * snap_cas)
    print(f"Cascade snap-up at |s+|^2 = {snap_cas:.3f}, "
          f"bias at |s+|^2 = {0.75*snap_cas:.3f}")

    p_nat = TCMTParams(Domega=Domega, gamma1=gamma1, gamma1_rad=gamma1_rad,
                       alpha3=alpha3_native, delta1=delta1, beta=0.0)
    snap_nat = find_snap_up(p_nat, s_max=100.0, n=3000)
    s_bias_nat = np.sqrt(0.75 * snap_nat)
    print(f"Native (alpha_3={alpha3_native}) snap-up at |s+|^2 = {snap_nat:.3f}, "
          f"bias at |s+|^2 = {0.75*snap_nat:.3f}")
    print(f"Bistability-threshold ratio nat/cas = {snap_nat/snap_cas:.2f} (predicted 1+R)")

    # ---- Panel (a): time traces at minimum switching energy ---------
    sigma_t = 10.0
    # Each architecture is biased close to its own snap-up; the
    # required pulse amplitude scales roughly with sqrt(bias) so
    # we set the grid wide enough for the native case (~sqrt(snap_nat)).
    pulse_max = max(2.0, 1.2 * np.sqrt(snap_nat))
    pulse_grid = np.linspace(0.02, pulse_max, 80)
    s_peak_cas = find_switching_threshold(p_cas, s_bias_cas, sigma_t, pulse_grid)
    s_peak_nat = find_switching_threshold(p_nat, s_bias_nat, sigma_t, pulse_grid)
    print(f"Cascade minimum switching pulse amplitude : {s_peak_cas:.3f}")
    print(f"Native  minimum switching pulse amplitude : {s_peak_nat:.3f}")

    if s_peak_cas is None or s_peak_nat is None:
        print("WARNING: could not find switching threshold within sweep")
        return

    # Switching energy ~ s_peak^2 * sqrt(pi) * sigma_t
    E_cas = s_peak_cas ** 2 * np.sqrt(np.pi) * sigma_t
    E_nat = s_peak_nat ** 2 * np.sqrt(np.pi) * sigma_t

    t0 = 5 * sigma_t + 30 / p_cas.gamma1
    t_total = t0 + 250.0 / p_cas.gamma1
    t_c, n1_c = switch_with_pulse(p_cas, s_bias_cas, s_peak_cas, sigma_t, t0, t_total)
    t_n, n1_n = switch_with_pulse(p_nat, s_bias_nat, s_peak_nat, sigma_t, t0, t_total)

    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.7))

    axs[0].plot(t_c, n1_c, color="C0", lw=1.7, label="cascade")
    axs[0].plot(t_n, n1_n, color="C2", lw=1.7, label="native")
    axs[0].axvline(t0, color="0.5", ls="-.", lw=0.9)
    axs[0].text(t0 + 5, 0.5 * n1_c.max(), "set pulse",
                fontsize=9.5, color="0.4")
    axs[0].set_xlabel(r"$t \cdot \gamma_2$")
    axs[0].set_ylabel(r"$|a_1|^2$  (intracavity FM energy)")
    axs[0].set_title("(a) Switching events at minimum control-pulse energy")
    axs[0].legend(loc="lower right", fontsize=9.5)
    axs[0].grid(True, alpha=0.3)

    # ---- Panel (b): pulse-energy sweep, find threshold for each arch ----
    # Run each architecture at its OWN bias; sweep the control-pulse
    # amplitude and plot the late-time |a_1|^2 to find the threshold.
    pulse_axis = np.linspace(0.02, pulse_max, 60)
    n1_final_cas = np.zeros_like(pulse_axis)
    n1_final_nat = np.zeros_like(pulse_axis)
    for i, sp in enumerate(pulse_axis):
        _, n1c = switch_with_pulse(p_cas, s_bias_cas, sp, sigma_t,
                                    t0, t_total)
        _, n1n = switch_with_pulse(p_nat, s_bias_nat, sp, sigma_t,
                                    t0, t_total)
        n1_final_cas[i] = np.mean(n1c[int(0.92 * len(n1c)):])
        n1_final_nat[i] = np.mean(n1n[int(0.92 * len(n1n)):])

    axs[1].plot(pulse_axis ** 2 * np.sqrt(np.pi) * sigma_t,
                n1_final_cas, "o-", color="C0", ms=4, mfc="none", mew=1.0,
                label="cascade")
    axs[1].plot(pulse_axis ** 2 * np.sqrt(np.pi) * sigma_t,
                n1_final_nat, "s-", color="C2", ms=4, mfc="none", mew=1.0,
                label="native")
    axs[1].axvline(E_cas, color="C0", ls=":", lw=1.0)
    axs[1].axvline(E_nat, color="C2", ls=":", lw=1.0)
    axs[1].set_xlabel(r"control-pulse energy  $\sim s_{\rm peak}^2 \sigma_t$  (arb.)")
    axs[1].set_ylabel(r"late-time $|a_1|^2$ on bias")
    axs[1].set_title("(b) Threshold sweep: cascade switches at much lower energy")
    axs[1].set_xscale("log")
    axs[1].grid(True, which="both", alpha=0.3)
    axs[1].legend(loc="upper left", fontsize=9.5)

    fig.suptitle(
        rf"All-optical switching: cascade switches at {E_cas/E_nat*100:.1f}\% "
        rf"of the native energy (ratio $= 1/(1+R)$)",
        fontsize=11.5, y=1.00,
    )
    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "switching_energy"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    print(f"  E_switch ratio cascade/native = {E_cas/E_nat:.3f}")
    plt.close(fig)


if __name__ == "__main__":
    main()
