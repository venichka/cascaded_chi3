"""Phase-bandwidth locus: native vs cascaded architectures.

Following Section 12 of the main note: per-photon nonlinear phase shift
phi_NL/|s_+|^2 scales as

  Architecture A (native, sweep Q_1):
        phi_NL ~ Q_1^2 * chi^(3)_native,   bandwidth ~ 1/Q_1
        => slope -2 on log(phi/pump) vs log(bandwidth).

  Architecture B_II (cascade, fix Q_1, sweep Q_2):
        phi_NL ~ Q_1^2 * Q_2 * |chi^(2)|^2 eta / (n_2^2 chi^(3)_native) * chi^(3)_native,
        bandwidth ~ 1/Q_2
        => slope -1. (The cascaded sum rule.)

  Architecture B_I (cascade, fix Q_2, sweep Q_1):
        phi_NL ~ Q_1^2 * (fixed cascaded coefficient),
        bandwidth ~ 1/Q_1 (FM bottlenecks now)
        => slope -2, parallel to A, vertically offset by the structural
        ratio R = |chi^(3)_casc|/|chi^(3)_native|.

The cascade is thus a 2D *region* in the (bandwidth, phi/pump) plane,
bounded by B_I on the high-phase side and B_II on the wide-bandwidth
side. The native locus A is a single line. The vertical gap between
A and B_I is the structural enhancement R, which depends only on the
host's chi^(2)^2 / chi^(3)_native ratio.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from params import ETA_GEOM, HOSTS, OMEGA1_REAL

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"


def enhancement_ratio(Q2: float, host: dict, eta: float = ETA_GEOM) -> float:
    """R = |chi^(3)_eff,casc| / |chi^(3)_native| at Delta omega = 2 gamma_2.

    The prefactor 1/24 comes from omega_1/Delta omega = Q_2/4 at this
    conservative operating point (gamma_2 = omega_2/Q_2 = 2 omega_1/Q_2).
    """
    return (Q2 / (24.0 * host["n"] ** 2)) * eta * host["chi2"] ** 2 / host["chi3"]


def phi_per_pump_native(Q1: float, host: dict) -> float:
    """Up to common geometric prefactor (cancels in ratios).

    phi/pump ~ Q_1^2 * chi^(3)_native / n^4 (the n^4 factor comes from
    |f_1|^4 mode normalisation that does cancel for the cascade in the
    final formula but does NOT cancel for the native architecture --
    native is fully gated by n_1 since only one mode is involved).
    """
    return Q1 ** 2 * host["chi3"] / host["n"] ** 4


def phi_per_pump_cascade(Q1: float, Q2: float, host: dict) -> float:
    """phi/pump for cascade. Same Q_1^2 FM build-up factor; cascaded
    coefficient ~ Q_2 * chi^(2)^2 eta / (n_2^2 chi^(3)_native) * chi^(3)_native
    relative to native, so phi_cascade = R(Q_2) * phi_native(Q_1)."""
    return enhancement_ratio(Q2, host) * phi_per_pump_native(Q1, host)


def bandwidth_native(Q1: float, omega1: float = OMEGA1_REAL) -> float:
    return omega1 / (2 * np.pi * Q1)  # gamma_1/(2 pi)


def bandwidth_casc(Q2: float, Q1: float | None = None,
                   omega1: float = OMEGA1_REAL) -> float:
    """Bandwidth = min(gamma_1, gamma_2) / (2 pi). With omega_2 = 2 omega_1:
    gamma_1 = omega_1/Q_1, gamma_2 = 2 omega_1/Q_2."""
    g2 = 2 * omega1 / Q2
    if Q1 is None:
        return g2 / (2 * np.pi)
    g1 = omega1 / Q1
    return min(g1, g2) / (2 * np.pi)


def truncate_at_bottleneck(bw_arr, phi_arr, Q_other, vary="Q2",
                           omega1: float = OMEGA1_REAL):
    """Keep only points where bandwidth is naturally set by the SWEPT Q
    (not bottlenecked by the OTHER fixed Q). Removes the artificial
    vertical cliffs that appear when min(gamma_1, gamma_2) flips which
    gamma is the bottleneck.

    vary='Q2': sweep Q_2; cascade-limited (B_II valid) when gamma_2 <
        gamma_1, i.e. BW < gamma_1/(2 pi). Use Q_other = Q_1 (fixed).
    vary='Q1': sweep Q_1; FM-limited (B_I valid) when gamma_1 < gamma_2,
        i.e. BW < gamma_2/(2 pi). Use Q_other = Q_2 (fixed).
    """
    g_other = omega1 / Q_other if vary == "Q2" else 2 * omega1 / Q_other
    bw_ceiling = g_other / (2 * np.pi)
    mask = bw_arr < bw_ceiling
    return bw_arr[mask], phi_arr[mask]


def main():
    Q1_grid_A = np.geomspace(10, 1e5, 80)            # family A: vary Q_1
    Q2_grid_BII = np.geomspace(100, 1e5, 80)          # family B_II: vary Q_2
    Q1_grid_BI = np.geomspace(10, 1e5, 80)            # family B_I: vary Q_1

    Q1_fixed_BII = 1_000      # fixed Q_1 for family B_II
    Q2_fixed_BI = 4_000       # fixed Q_2 for family B_I

    # Restrict to materials with chi^(2) > 0 (Si, Si3N4 dropped here:
    # cascade not available; native curve trivial).
    cascade_hosts = {k: v for k, v in HOSTS.items() if v["chi2"] > 0}

    fig, axs = plt.subplots(1, 2, figsize=(13.0, 5.0), gridspec_kw={"wspace": 0.25})

    for ax, hosts_to_show, title in [
        (axs[0], {"LiNbO3": cascade_hosts["LiNbO3"]}, r"(a) LiNbO$_3$: cascade $R\simeq 13$ at $Q_2=4{,}000$"),
        (axs[1], {"AlGaAs": cascade_hosts["AlGaAs"]}, "(b) AlGaAs: native-Kerr wins at any practical $Q_2$"),
    ]:
        for name, host in hosts_to_show.items():
            color = host["color"]
            # Family A: vary Q_1
            bw_A = bandwidth_native(Q1_grid_A)
            phi_A = np.array([phi_per_pump_native(Q1, host) for Q1 in Q1_grid_A])
            ax.loglog(bw_A, phi_A, "-", color="C2", lw=2.2,
                      label=r"A: native (vary $Q_1$, slope $-2$)")

            # Family B_II: fix Q_1, vary Q_2. Truncate where gamma_1
            # becomes the bottleneck (BW no longer set by Q_2).
            bw_BII = np.array([bandwidth_casc(Q2, Q1_fixed_BII) for Q2 in Q2_grid_BII])
            phi_BII = np.array([phi_per_pump_cascade(Q1_fixed_BII, Q2, host)
                                for Q2 in Q2_grid_BII])
            bw_BII, phi_BII = truncate_at_bottleneck(bw_BII, phi_BII,
                                                    Q1_fixed_BII, vary="Q2")
            ax.loglog(bw_BII, phi_BII, "--", color="C0", lw=2.2,
                      label=rf"B$_{{\rm II}}$: cascade, fix $Q_1={Q1_fixed_BII}$, "
                            rf"vary $Q_2$ (slope $-1$)")

            # Family B_I: fix Q_2, vary Q_1. Truncate where gamma_2
            # becomes the bottleneck.
            bw_BI = np.array([bandwidth_casc(Q2_fixed_BI, Q1) for Q1 in Q1_grid_BI])
            phi_BI = np.array([phi_per_pump_cascade(Q1, Q2_fixed_BI, host)
                               for Q1 in Q1_grid_BI])
            bw_BI, phi_BI = truncate_at_bottleneck(bw_BI, phi_BI,
                                                  Q2_fixed_BI, vary="Q1")
            ax.loglog(bw_BI, phi_BI, ":", color="C3", lw=2.2,
                      label=rf"B$_{{\rm I}}$: cascade, fix $Q_2={Q2_fixed_BI}$, "
                            rf"vary $Q_1$ (slope $-2$)")

            # Annotate the structural ratio R at the matched-Q corner
            Q_match = Q2_fixed_BI
            R_val = enhancement_ratio(Q_match, host)
            ax.annotate(rf"vertical offset $= R \approx {R_val:.0f}$"
                        if R_val > 1 else rf"$R \approx {R_val:.2f}$ (native wins)",
                        xy=(0.04, 0.95), xycoords="axes fraction",
                        fontsize=9, color="0.3",
                        verticalalignment="top")

        ax.set_xlabel(r"operating bandwidth  $\Delta f$  [Hz]")
        ax.set_ylabel(r"$\phi_{\rm NL}/|s_+|^2$  (arb.)")
        ax.set_title(title)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="lower left", fontsize=8.5)
        ax.set_xlim(1e9, 5e13)
        ax.set_ylim(1e-21, 1e-7)

    fig.suptitle(
        "Phase-bandwidth locus: native (A) is one line, cascade is a 2D region (B$_{\\rm I}$, B$_{\\rm II}$)",
        fontsize=11.5, y=1.00,
    )
    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "phase_shift_per_pump"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"wrote {out}.{{png,pdf}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
