"""Three-way trade-off: bulk vs resonant native vs cascaded.

All three regimes plotted on the same (bandwidth, enhancement) plane,
where enhancement is the *effective* chi^(3) of the structure
normalised to the host's bulk native chi^(3). For each scheme:

  - Bulk phase-mismatched PPLN: SVEA cascade gives
    chi^(3)_casc,bulk = -(2 omega_1 / (3 c n_2 Delta_k)) |chi^(2)|^2 F(Delta_k L),
    bandwidth ~ |Delta_k| v_g / (2 pi) ~ few THz.
    (NOTE: the formula has 1/n_2 only -- the n_1 powers from the FM
    SVEA projection cancel against the n_1 in the cascade kernel,
    paralleling the single-resonance metasurface result.)
  - Resonant native (single FM cavity, no SH): chi^(3)_eff = chi^(3)_native,
    so enhancement = 1; bandwidth = gamma_1 / (2 pi) = omega_1 / (2 pi Q_1).
    The cavity boosts intracavity intensity (by Q_1^2) but does not
    change the local chi^(3) of the host material.
  - Cascaded doubly-resonant metasurface:
    enhancement = (Q_2 / (12 n_2^2)) eta_geom |chi^(2)|^2 / chi^(3)_native,
    bandwidth = gamma_2 / (2 pi) = 2 omega_1 / (2 pi Q_2).

Three takeaways visible at a glance:
  1. Resonant native sits on the enhancement = 1 line for any Q_1 --
     a horizontal slice along the bandwidth axis (no chi^(3) boost).
  2. Cascaded sits on the same hyperbola as bulk -- the sum rule is
     identical for both, just operated at different points.
  3. The metasurface route is only worth the effort when the host's
     chi^(2)^2 / chi^(3)_native ratio is favourable (see
     fig_cascaded_vs_native.py for the host-dependent crossover).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from style import apply_style
apply_style()
import numpy as np

from params import (CHI2_REAL, ETA_GEOM, N2_REAL, OMEGA1_REAL)

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures"

CHI3_NATIVE = 2e-21  # m^2/V^2, LiNbO3 native Kerr (Section 9)


def metasurface_point(Q2):
    enh = (Q2 / (12.0 * N2_REAL ** 2)) * ETA_GEOM * CHI2_REAL ** 2 / CHI3_NATIVE
    gamma2 = 2.0 * OMEGA1_REAL / Q2
    bw_hz = gamma2 / (2 * np.pi)
    return bw_hz, enh


def bulk_pplnpoint(Delta_k_L, Delta_k=3.0e5, v_g=2.0e8, n2=N2_REAL):
    # 'Bulk PPLN at the same |chi^(2)|': note Section 10 / eq:chi3casc-bulk.
    # chi^(3)_casc,bulk = -2 omega_1 / (3 c n_2 Delta_k) * |chi^(2)|^2 * F(Delta_k L).
    # F -> 1 for Delta_k L >> 1 (broadband cascading); F -> Delta_k L / 2
    # for Delta_k L << 1 (SHG conversion competes).
    c = 3.0e8
    # Phase factor F(Delta_k L), interpolated between the two limits:
    F = 1.0 - np.sin(Delta_k_L) / Delta_k_L * np.cos(Delta_k_L)
    chi3_bulk = abs(2.0 * OMEGA1_REAL * CHI2_REAL ** 2 * F
                    / (3.0 * c * n2 * Delta_k))
    enh = chi3_bulk / CHI3_NATIVE
    # Bandwidth: |Delta_k| v_g / (2 pi). For Delta_k = 3e5 /m, v_g ~ c/n
    # ~ 1.3e8 m/s, this is ~6 THz -- the "~10 THz" of the note.
    bw_hz = abs(Delta_k) * v_g / (2 * np.pi)
    return bw_hz, enh


def native_resonant_point(Q1: float):
    """Resonant native chi^(3) surface: enhancement = 1 (no chi^(3) boost),
    bandwidth set by FM linewidth gamma_1 / (2 pi) = omega_1/(2 pi Q_1)."""
    bw_hz = OMEGA1_REAL / (2 * np.pi * Q1)
    enh = 1.0
    return bw_hz, enh


def main():
    Q2_list = [1_000, 4_000, 10_000, 30_000, 100_000]
    metasurface_pts = np.array([metasurface_point(Q2) for Q2 in Q2_list])

    Q1_list = [100, 1_000, 10_000, 100_000]
    native_pts = np.array([native_resonant_point(Q1) for Q1 in Q1_list])

    bulk_DkL_list = [0.5, 1.0, 3.0]
    bulk_pts = np.array([bulk_pplnpoint(d) for d in bulk_DkL_list])

    # Sum-rule hyperbola: pick the constant from the metasurface points
    # (they are self-consistent within the formula). Constant K such that
    # enhancement * bw_hz = K for the metasurface family.
    K = metasurface_pts[0, 0] * metasurface_pts[0, 1]
    bw_grid = np.geomspace(1e9, 3e13, 400)
    enh_curve = K / bw_grid

    fig, ax = plt.subplots(figsize=(8.0, 5.4))

    ax.loglog(bw_grid, enh_curve, "-", color="grey", lw=1.5,
              label=r"$|\chi^{(3)}_{\rm casc}|/|\chi^{(3)}_{\rm native}| \cdot \Delta f = $ const")
    ax.axhline(1.0, color="k", lw=0.6, ls=":")
    ax.text(2e9, 1.15, r"$|\chi^{(3)}_{\rm casc}| = |\chi^{(3)}_{\rm native}|$",
            fontsize=9, color="k")

    ax.loglog(metasurface_pts[:, 0], metasurface_pts[:, 1], "o", color="C0",
              ms=10, mec="k", mew=1.0, label="metasurface (varying $Q_2$)")
    for (bw, enh), Q2 in zip(metasurface_pts, Q2_list):
        ax.annotate(f"$Q_2={Q2}$", xy=(bw, enh),
                    xytext=(bw * 0.55, enh * 1.6), fontsize=8.5, color="C0")

    ax.loglog(bulk_pts[:, 0], bulk_pts[:, 1], "s", color="C3",
              ms=10, mec="k", mew=1.0, label=r"bulk PPLN (phase mismatch)")
    for (bw, enh), DkL in zip(bulk_pts, bulk_DkL_list):
        ax.annotate(rf"$\Delta k\,L={DkL:g}$", xy=(bw, enh),
                    xytext=(bw * 1.05, enh * 1.5), fontsize=8.5, color="C3")

    # Resonant native: enhancement = 1, bandwidth = gamma_1 / (2 pi)
    ax.loglog(native_pts[:, 0], native_pts[:, 1], "D", color="C2",
              ms=9, mec="k", mew=1.0, label="resonant native (varying $Q_1$)")
    for (bw, enh), Q1 in zip(native_pts, Q1_list):
        # Alternate above/below so labels don't overlap data
        y_off = 0.45 if Q1_list.index(Q1) % 2 == 0 else 2.0
        ax.annotate(f"$Q_1={Q1}$", xy=(bw, enh),
                    xytext=(bw, enh * y_off), fontsize=8.5, color="C2",
                    ha="center")

    ax.set_xlabel(r"usable bandwidth  $\Delta f$  [Hz]")
    ax.set_ylabel(r"$|\chi^{(3)}_{\rm eff}| / |\chi^{(3)}_{\rm native}|$")
    ax.set_title("Three-way trade-off: bulk vs resonant native vs cascaded")
    ax.set_xlim(1e8, 3e13)
    ax.set_ylim(2e-2, 1e3)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower left", fontsize=10)

    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "bulk_vs_metasurface"
    fig.savefig(out.with_suffix(".png"), dpi=150)
    fig.savefig(out.with_suffix(".pdf"))

    print(f"wrote {out}.{{png,pdf}}")
    print(f"  sum-rule constant K = {K:.3e}  (enh * BW_Hz)")
    for (bw, enh), Q2 in zip(metasurface_pts, Q2_list):
        print(f"    Q_2={Q2:>6}: enh = {enh:6.1f}, BW = {bw/1e9:7.1f} GHz, prod = {bw*enh:.3e}")
    for (bw, enh), DkL in zip(bulk_pts, bulk_DkL_list):
        print(f"    bulk DkL={DkL:.2f}: enh = {enh:6.2f}, BW = {bw/1e12:.2f} THz, prod = {bw*enh:.3e}")
    for (bw, enh), Q1 in zip(native_pts, Q1_list):
        print(f"    native Q_1={Q1:>6}: enh = 1.00, BW = {bw/1e9:7.1f} GHz")
    plt.close(fig)


if __name__ == "__main__":
    main()
