"""Optical bistability driven by the cascaded effective Kerr.

For positive Domega the cascaded Kerr is negative (self-defocusing-like):
the effective FM resonance shifts *down* as |a1|^2 grows. To get a
hysteresis loop, place the pump *above* the FM resonance (delta1 > 0,
"blue-detuned"); as input power grows, the shifted resonance climbs
toward the pump, the response runs away, and the system snaps to the
upper branch. Reducing the power lets it snap back at a lower threshold.

Method: stepwise quasi-static sweep. At each input level we integrate
the full TCMT to steady state, initialised from the previous level's
final state. This traces the lower branch on the way up (snapping at
the upper saddle-node) and the upper branch on the way back down
(snapping at the lower saddle-node), without the transient ringing of
a continuous time-domain ramp.
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


def sweep_both_ways(p: TCMTParams, s2_grid: np.ndarray, t_settle: float):
    """Quasi-static up-then-down hysteresis sweep, returning (n1_up, n2_up,
    n1_dn, n2_dn)."""
    n_steps = len(s2_grid)
    n1_up, n2_up = np.zeros(n_steps), np.zeros(n_steps)
    n1_dn, n2_dn = np.zeros(n_steps), np.zeros(n_steps)
    state = (0 + 0j, 0 + 0j)
    for i, sq in enumerate(s2_grid):
        sol = integrate(p, np.sqrt(sq), (0.0, t_settle), y0=state,
                        max_step=2.0)
        a1, a2 = sol.a[0, -1], sol.a[1, -1]
        state = (a1, a2)
        n1_up[i] = abs(a1) ** 2
        n2_up[i] = abs(a2) ** 2
    for j, sq in enumerate(s2_grid[::-1]):
        sol = integrate(p, np.sqrt(sq), (0.0, t_settle), y0=state,
                        max_step=2.0)
        a1, a2 = sol.a[0, -1], sol.a[1, -1]
        state = (a1, a2)
        n1_dn[n_steps - 1 - j] = abs(a1) ** 2
        n2_dn[n_steps - 1 - j] = abs(a2) ** 2
    return n1_up, n2_up, n1_dn, n2_dn


def find_snap_up_threshold(p: TCMTParams, s2_max: float, n_steps: int,
                           t_settle: float) -> float:
    """Quasi-static up-sweep; return |s_+|^2 at the largest jump in |a_1|^2."""
    s2_grid = np.linspace(1e-4, s2_max, n_steps)
    state = (0 + 0j, 0 + 0j)
    n1 = np.zeros(n_steps)
    for i, sq in enumerate(s2_grid):
        sol = integrate(p, np.sqrt(sq), (0.0, t_settle), y0=state, max_step=2.0)
        a1, a2 = sol.a[0, -1], sol.a[1, -1]
        state = (a1, a2)
        n1[i] = abs(a1) ** 2
    # largest jump in |a_1|^2 marks the snap-up
    di = np.diff(n1)
    idx = int(np.argmax(di))
    return float(s2_grid[idx + 1])


def main():
    # Push deep into the cascaded regime so saturation (|beta a1| ~ |Domega|)
    # is well-separated from the bistability cusp.
    p = TCMTParams(Domega=5.0, gamma1=0.4, gamma1_rad=0.2, alpha3=0.0,
                   delta1=0.8)  # blue-detuned by 4*(gamma1/2)

    # Locate the analytic bistable band, then bracket it with margin.
    s2_axis = np.linspace(0.0, 6.0, 800)
    n_roots = np.array([len(cw_kerr_only_intracavity(p, np.sqrt(s)))
                        for s in s2_axis])
    bist_band = s2_axis[n_roots >= 3]
    s2_max = min(6.0, 1.5 * bist_band.max()) if bist_band.size else 6.0

    n_steps = 70
    s2_grid = np.linspace(1e-3, s2_max, n_steps)
    t_settle = 250.0 / min(p.gamma1, 1.0)

    n1_up, n2_up, n1_dn, n2_dn = sweep_both_ways(p, s2_grid, t_settle)

    # ---- analytic Kerr+2PA cubic branches -------------------------------
    s2_an = np.linspace(0.0, s2_max, 800)
    branches = [cw_kerr_only_intracavity(p, np.sqrt(s)) for s in s2_an]

    # ---- Panel (d) data: clean threshold-vs-R sweep at SHARED params -----
    # Fix Delta omega = 5 gamma_2 (cascade gives self-defocusing,
    # alpha_eff_casc = -0.198) and pick alpha_3 = -0.05 (self-defocusing
    # native) so the cascade ADDS constructively to native. With blue-
    # detuned pump (delta_p > 0) and combined alpha_eff < 0, bistability
    # works in the standard self-defocusing geometry for all beta. Sweep
    # beta to vary R = |alpha_eff_casc| / |alpha_3|. Threshold scales as
    # 1 / (1 + R) (cascade lowers snap-up power proportionally).
    Domega_R = 5.0
    alpha3_R = -0.05  # bare Kerr (negative, self-defocusing) -- bistable alone
    delta_R  = 0.8
    # alpha_eff_casc at Delta omega = 5 with beta = 1: -5/(25+0.25) = -0.198
    # combined alpha_eff = -0.05 + (-0.198 beta^2)  (both negative)
    # R(beta) = |alpha_eff_casc| / |alpha_3| = 0.198 * beta^2 / 0.05
    #        = 3.96 * beta^2
    beta_grid = np.array([0.0, 0.3, 0.55, 0.8, 1.0, 1.25])
    R_grid = beta_grid**2 * 0.198 / abs(alpha3_R)
    threshold_R = np.zeros_like(beta_grid)
    for k, b in enumerate(beta_grid):
        p_k = TCMTParams(Domega=Domega_R, gamma1=0.4, gamma1_rad=0.2,
                         alpha3=alpha3_R, delta1=delta_R, beta=float(b))
        # Bracket the cusp: at beta=0 it sits at moderate |s+|^2; the
        # cusp shrinks as 1/(1+R) so reduce the sweep range correspondingly.
        if b == 0:
            s2_max_k = 30.0
        else:
            s2_max_k = 30.0 / (1 + R_grid[k]) * 2.5   # ~2.5x predicted threshold
        threshold_R[k] = find_snap_up_threshold(p_k, s2_max_k, 80, t_settle)

    fig, axs = plt.subplots(2, 2, figsize=(13.0, 9.2))
    axs = axs.ravel()

    low, mid, high = [], [], []
    for s, roots in zip(s2_an, branches):
        rs = sorted(roots)
        if len(rs) == 1:
            low.append((s, rs[0]))
        elif len(rs) >= 3:
            low.append((s, rs[0]))
            mid.append((s, rs[1]))
            high.append((s, rs[2]))
    if low:
        a = np.array(low)
        axs[0].plot(a[:, 0], a[:, 1], "-", color="grey", lw=1.0,
                    label="analytic stable")
    if high:
        a = np.array(high)
        axs[0].plot(a[:, 0], a[:, 1], "-", color="grey", lw=1.0)
    if mid:
        a = np.array(mid)
        axs[0].plot(a[:, 0], a[:, 1], "--", color="grey", lw=1.0,
                    label="analytic unstable")

    axs[0].plot(s2_grid, n1_up, "o-", color="C0", ms=3.5, mfc="white",
                mew=1.0, label="TCMT, sweep up")
    axs[0].plot(s2_grid, n1_dn, "s-", color="C3", ms=3.5, mfc="white",
                mew=1.0, label="TCMT, sweep down")
    axs[0].set_xlabel(r"$|s_+|^2$  (input power, normalised)")
    axs[0].set_ylabel(r"$|a_1|^2$  (intracavity FM energy)")
    axs[0].set_title("(a) FM hysteresis loop")
    axs[0].legend(loc="lower right", fontsize=9)
    axs[0].grid(True, alpha=0.3)

    axs[1].plot(s2_grid, n2_up, "o-", color="C0", ms=3.5, mfc="white",
                mew=1.0, label="sweep up")
    axs[1].plot(s2_grid, n2_dn, "s-", color="C3", ms=3.5, mfc="white",
                mew=1.0, label="sweep down")
    axs[1].set_xlabel(r"$|s_+|^2$  (input power, normalised)")
    axs[1].set_ylabel(r"$|a_2|^2$  (intracavity SH energy)")
    axs[1].set_title("(b) SH amplitude jumps in lockstep with FM")
    axs[1].legend(loc="upper left", fontsize=9)
    axs[1].grid(True, alpha=0.3)

    # ---- Panel (c): threshold vs R, validating the 1/(1+R) prediction ----
    R_fine = np.linspace(0, R_grid.max() * 1.05, 200)
    pred = threshold_R[0] / (1.0 + R_fine)
    axs[2].plot(R_grid, threshold_R, "o", color="C3", ms=9, mfc="white",
                mew=1.6, label=r"TCMT (snap-up threshold)")
    axs[2].plot(R_fine, pred, "-", color="k", lw=1.4,
                label=r"$|s_+|^{2,{\rm A+B}}_{\rm bist} / |s_+|^{2,{\rm A}}_{\rm bist} = 1/(1+R)$")
    axs[2].set_xlabel(r"$R = |\alpha_3^{\rm casc}| / \alpha_3$")
    axs[2].set_ylabel(r"snap-up threshold  $|s_+|^2_{\rm bist}$")
    axs[2].set_title(r"(c) Threshold drops as $1/(1+R)$  (cascade vs native, same $\alpha_3$)")
    axs[2].grid(True, alpha=0.3)
    axs[2].legend(loc="upper right", fontsize=9.5)
    # Annotate the R=0 point ("native only")
    axs[2].annotate(r"$\beta=0$ (native only)",
                    xy=(0, threshold_R[0]),
                    xytext=(R_grid.max() * 0.18, threshold_R[0] * 1.05),
                    fontsize=9, color="C2",
                    arrowprops=dict(arrowstyle="->", color="C2", lw=0.9))

    # ---- Panel (d): threshold ratio (TCMT/prediction), residuals ----
    ratio = threshold_R / threshold_R[0]      # normalised
    pred_ratio = 1.0 / (1.0 + R_grid)
    axs[3].semilogy(R_grid, ratio, "o", color="C3", ms=9, mfc="white", mew=1.6,
                    label=r"TCMT  $|s_+|^2_{\rm bist}(R) / |s_+|^2_{\rm bist}(0)$")
    axs[3].semilogy(R_fine, 1.0 / (1.0 + R_fine), "-", color="k", lw=1.4,
                    label=r"$1/(1+R)$  prediction")
    axs[3].set_xlabel(r"$R$")
    axs[3].set_ylabel(r"normalised snap-up threshold")
    axs[3].set_title(r"(d) Normalised view: orders-of-magnitude drop on log axis")
    axs[3].grid(True, which="both", alpha=0.3)
    axs[3].legend(loc="upper right", fontsize=9.5)
    # Print residuals for sanity check (also captured in stdout log)
    for k, (Rk, rk, pk) in enumerate(zip(R_grid, ratio, pred_ratio)):
        print(f"  beta={beta_grid[k]:.2f}, R={Rk:.2f}, "
              f"meas={rk:.3f}, 1/(1+R)={pk:.3f}, "
              f"residual={(rk-pk)/pk*100:+.1f}%")

    fig.suptitle(
        rf"Cascaded Kerr bistability  ($\Delta\omega/\gamma_2={p.Domega}$, "
        rf"$\delta_1/\gamma_2={p.delta1}$, blue-detuned).  Panels (c,d) "
        r"validate the $1/(1+R)$ threshold scaling across a range of $R$.",
        fontsize=11.5, y=1.00,
    )
    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "bistability"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
