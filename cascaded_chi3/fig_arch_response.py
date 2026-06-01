"""Side-by-side TCMT response of two metasurface architectures.

  A. Native-Kerr resonance (single FM cavity, alpha_3 = alpha_3_native,
     beta = 0).
  B. Cascaded chi^(2):chi^(2) (FM + SH cavities, same alpha_3 retained,
     beta = 1).

Both run through the *same* TCMT under *identical* CW input and the
*same* FM cavity (same gamma_1, gamma_1_rad, delta_p). The only
difference between the two runs is the SH coupling beta. This isolates
the cascaded contribution from the FM cavity build-up which is common
to both schemes.

Two panels:
  (a) Intracavity FM energy |a_1|^2 vs incident pump |s_+|^2.
  (b) Magnitude of the steady-state nonlinear phase shift on a_1.

Story: at low input, the cascaded curve outperforms native by the
host-dependent enhancement factor (~25x for LiNbO_3 at Q_2=4000).
At higher input the cascaded curve saturates (parametric back-action,
|beta a_1| -> |Delta omega|) while native keeps climbing linearly.
The crossover is the regime where the architecture choice flips.
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
from tcmt import attractor, integrate

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def attractor_phi_and_n1(p: TCMTParams, sp: complex,
                         t_settle: float = 200.0, t_avg: float = 200.0):
    """Full-TCMT attractor: returns (<|a_1|^2>, <|a_2|^2>, phi_avg, is_limit_cycle).

    Time-averages over t_avg after letting transients die for t_settle.
    For a stable fixed point this is just the FP value; for the
    parametric-instability regime it is the limit-cycle time average."""
    res = attractor(p, complex(sp), t_settle=t_settle, t_avg=t_avg,
                    n_samples=800, max_step=0.5)
    return (res["n1_mean"], res["n2_mean"],
            float(np.angle(res["a1_last"])),
            res["is_limit_cycle"])


def main():
    # LiNbO_3 worked-example design point (Section 9 of the note):
    #   |Delta omega| = 2 gamma_2 (conservative, fully adiabatic),
    #   gamma_1/gamma_2 = 10, gamma_1_rad/gamma_1 = 0.5.
    # We pick Delta omega = -2 (FM pump red-detuned from omega_2/2) so
    # that the cascade contribution has the SAME sign as native alpha_3
    # > 0 -- the "constructive" configuration where the cascade enhances
    # rather than cancels native, giving the clean (1+R) phase ratio.
    # The bare-Kerr value alpha_3_native = 0.04 is calibrated so that the
    # full-Lorentzian cascade coefficient |alpha_3^casc| = 0.471 (beta=1,
    # gamma_2=1, |Delta omega|=2) yields R = 0.471/0.04 = 11.8, matching
    # the worked-example R = 12 from |chi3_casc| = Q_2/(24 n_2^2) eta
    # |chi^(2)|^2.
    alpha_3_native = 0.04
    base = dict(Domega=-2.0, gamma1=10.0, gamma1_rad=5.0,
                alpha3=alpha_3_native, delta1=0.0)
    p_nat = TCMTParams(**base, beta=0.0)
    p_cas = TCMTParams(**base, beta=1.0)

    # Coarser pump grid (32 points instead of 70) -- attractor-averaging
    # is more expensive than single steady_state calls.
    s2_grid = np.geomspace(1e-3, 200.0, 32)

    n1_nat = np.empty_like(s2_grid)
    n1_cas = np.empty_like(s2_grid)
    phi_nat = np.empty_like(s2_grid)
    phi_cas = np.empty_like(s2_grid)
    is_lc_cas = np.zeros_like(s2_grid, dtype=bool)

    for i, s2 in enumerate(s2_grid):
        sp = np.sqrt(s2)
        n1_nat[i], _, phi_nat[i], _ = attractor_phi_and_n1(p_nat, sp)
        n1_cas[i], _, phi_cas[i], is_lc_cas[i] = attractor_phi_and_n1(p_cas, sp)

    # Linear-cavity reference for intracavity FM:
    g1 = 0.5 * p_nat.gamma1
    n1_linear = p_nat.gamma1_rad * s2_grid / g1**2

    # Saturation marker (cascaded only):
    n1_sat = p_cas.Domega ** 2  # |beta a_1| = |Domega| <=> |a_1|^2 = Domega^2
    unstable_cas = n1_cas > n1_sat

    fig, axs = plt.subplots(1, 2, figsize=(12.0, 4.8),
                            gridspec_kw={"wspace": 0.28})

    # ---- Panel (a): intracavity FM energy --------------------------------
    axs[0].loglog(s2_grid, n1_linear, "-", color="0.6", lw=1.4,
                  label="linear cavity (no NL)")
    axs[0].loglog(s2_grid, n1_nat, "o", color="C2", ms=5, mfc="C2", mew=0,
                  label=r"native ($\beta = 0$, $\widetilde\alpha_3=0.04$)")
    # Cascade: split stable (filled) vs past saturation (open)
    axs[0].loglog(s2_grid[~unstable_cas], n1_cas[~unstable_cas], "s",
                  color="C0", ms=5, mfc="C0", mew=0,
                  label=r"cascaded ($\beta = 1$, same $\widetilde\alpha_3$)")
    if unstable_cas.any():
        axs[0].loglog(s2_grid[unstable_cas], n1_cas[unstable_cas], "s",
                      color="C0", ms=5, mfc="white", mew=1.3, alpha=0.7,
                      label=r"  past $|\beta a_1|=|\Delta\omega|$")
    axs[0].axhline(n1_sat, color="C3", ls=":", lw=0.9,
                   label=r"$|\beta a_1|^2 = |\Delta\omega|^2$ (sat.)")
    axs[0].set_xlabel(r"$|s_+|^2$  (input power)")
    axs[0].set_ylabel(r"$\langle|a_1|^2\rangle_t$  (intracavity FM energy)")
    axs[0].set_title("(a) FM energy: same input, same FM cavity")
    axs[0].grid(True, which="both", alpha=0.3)
    axs[0].legend(loc="lower right", fontsize=8.5)

    # ---- Panel (b): nonlinear phase shift -------------------------------
    axs[1].loglog(s2_grid, np.abs(phi_nat), "o", color="C2", ms=5,
                  mfc="C2", mew=0, label=r"native: $|\phi_{\rm NL}|$")
    axs[1].loglog(s2_grid[~unstable_cas], np.abs(phi_cas[~unstable_cas]),
                  "s", color="C0", ms=5, mfc="C0", mew=0,
                  label=r"cascaded: $|\phi_{\rm NL}|$")
    if unstable_cas.any():
        axs[1].loglog(s2_grid[unstable_cas], np.abs(phi_cas[unstable_cas]),
                      "s", color="C0", ms=5, mfc="white", mew=1.3,
                      alpha=0.7, label=r"  past sat. (limit-cycle avg)")
    # small-signal slope guides (moved out of the data envelope)
    s_ref = np.array([1e-3, 1e0])
    axs[1].loglog(s_ref, 5e-3 * (s_ref / s_ref[0]) ** 1, ":", color="0.55", lw=0.9)
    axs[1].text(2e-3, 2e-2, "slope $+1$\n(small-signal)", fontsize=9, color="0.4")
    axs[1].set_xlabel(r"$|s_+|^2$  (input power)")
    axs[1].set_ylabel(r"$|\phi_{\rm NL}|$  (rad)")
    axs[1].set_title(r"(b) Nonlinear phase shift on $a_1$")
    axs[1].grid(True, which="both", alpha=0.3)
    axs[1].legend(loc="lower right", fontsize=9)

    fig.suptitle(
        r"Native vs cascade-added -- same FM cavity, same input "
        r"($\widetilde{\alpha}_3 = +0.04$, $\Delta\omega = -2\gamma_2$, "
        r"constructive cascade; $R\simeq 12$, ratio $\simeq 1+R\simeq 13$)",
        fontsize=11.5, y=1.02,
    )
    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "arch_response"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    print(f"  at |s+|^2 = 1.0  : phi_nat = {phi_nat[np.argmin(np.abs(s2_grid - 1.0))]:+.4f}, "
          f"phi_cas = {phi_cas[np.argmin(np.abs(s2_grid - 1.0))]:+.4f}")
    print(f"  at |s+|^2 = 10.0 : phi_nat = {phi_nat[np.argmin(np.abs(s2_grid - 10.0))]:+.4f}, "
          f"phi_cas = {phi_cas[np.argmin(np.abs(s2_grid - 10.0))]:+.4f}")
    ratio = np.abs(phi_cas[np.argmin(np.abs(s2_grid - 1.0))]) / \
            np.abs(phi_nat[np.argmin(np.abs(s2_grid - 1.0))])
    print(f"  cascaded/native phase ratio at |s+|^2 = 1: {ratio:.2f} "
          f"(expected ~12 from R = 0.471 / 0.04 at Delta omega = 2 gamma_2)")
    plt.close(fig)


if __name__ == "__main__":
    main()
