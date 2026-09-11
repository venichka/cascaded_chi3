"""Combined 5-row comparison: SH-only samples 1-4 + Design 6.

Places Design 6 alongside the four earlier LiNbO3 samples in one set of comparison
figures + a master table, for the summary deck. Reuses the p-keyed physics from
sample_effects / sample_pulsed; builds samples via sample_params.build_sample and
Design 6 via design6_params.build_design6.

IMPORTANT normalisation caveat (verified in the pedantic check):
  * samples 1-4 : eta_geom = (|overlap| / AREA )^2   -- 2-D cross-section overlap
  * Design 6    : eta_geom = (|overlap| / VOLUME)^2   -- 3-D volume overlap
Both are dimensionless SHG coherence factors with C_bridge = 1, but they come from
different integrals, so the ABSOLUTE cross-comparison of R across the 2-D/3-D boundary
rests on the shared O(1) C_bridge. Within each convention the ranking is robust.

    conda run -n ml_base python combined_effects.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from style import apply_style
apply_style()

from sample_params import build_sample, chi3_casc, eta_geom, CHI3_NATIVE
from sample_effects import phi_NL, I_1GW
from design6_params import build_design6
from sample_pulsed import reduction_factors

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures" / "combined"

COLORS = {1: "#1f77b4", 2: "#ff7f0e", 3: "#2ca02c", 4: "#d62728", 6: "#9467bd"}
TAG = {1: "S1", 2: "S2", 3: "S3", 4: "S4", 6: "D6"}
DOP = 2.0


def get_all():
    P = {i: build_sample(i) for i in (1, 2, 3, 4)}
    P[6] = build_design6()
    return P


def _save(fig, name):
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"{name}.{ext}")
    print(f"wrote figures/combined/{name}.{{png,pdf}}")
    plt.close(fig)


def _lbl(i, p):
    return f"{TAG[i]} ({p['lam1_um']*1e3:.0f} nm, $Q_2$={p['Q2']:.0f})"


def overlap_at(p, dop):
    f = p["f0"] + dop * p["g_THz"]
    f = min(max(f, p["freq_arr"].min()), p["freq_arr"].max())
    return float(p["ov_interp"](f))


# ----------------------------------------------------------------------
def master_table(P):
    print("\n=== MASTER TABLE (Domega = +2 gamma2, C_bridge = 1) ===")
    print("%-3s %8s %7s %7s %7s %9s %8s %9s %9s %9s" %
          ("smp", "lam1nm", "Q2", "ovlp", "eta", "|chi3|", "R", "Df[GHz]", "phi@1GW", "sig_w/w"))
    rows = {}
    for i, p in P.items():
        c3 = chi3_casc(p, DOP)
        R = abs(c3) / CHI3_NATIVE
        rows[i] = dict(R=R, ov=overlap_at(p, DOP), eta=eta_geom(p, DOP), c3=abs(c3),
                       phi=abs(phi_NL(p, DOP, I_1GW)), tol=1.0 / p["Q2"])
        print("%-3s %8.0f %7.0f %7.3f %7.4f %9.2e %8.2f %9.1f %9.2e %9.1e" %
              (TAG[i], p["lam1_um"] * 1e3, p["Q2"], rows[i]["ov"], rows[i]["eta"],
               rows[i]["c3"], R, p["g_THz"] * 1e3, rows[i]["phi"], rows[i]["tol"]))
    best = max(rows, key=lambda k: rows[k]["R"])
    print("  ideal-R winner: %s (R=%.2f).  Note S1-4 area-norm vs D6 volume-norm (C_bridge caveat)."
          % (TAG[best], rows[best]["R"]))
    return rows


def fig_R_comparison(P, rows):
    dop = np.linspace(-6, 6, 1201)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    for i, p in P.items():
        lw = 2.4 if i == 6 else 1.4
        axes[0].plot(dop, [abs(chi3_casc(p, d)) / CHI3_NATIVE for d in dop],
                     color=COLORS[i], lw=lw, label=_lbl(i, p))
    axes[0].axhline(1.0, color="0.4", lw=1.0, ls="--", label="native Kerr ($R=1$)")
    for d in (-DOP, DOP):
        axes[0].axvline(d, color="0.6", lw=0.8, ls=":")
    axes[0].axvspan(-1, 1, color="0.9", alpha=0.6)
    axes[0].set_xlabel(r"$\Delta\omega/\gamma_2$")
    axes[0].set_ylabel(r"$R^{(1\mathrm{R})}=|\chi^{(3)}_{\rm casc}|/\chi^{(3)}_{\rm native}$")
    axes[0].set_title("Cascaded enhancement vs detuning")
    axes[0].legend(fontsize=8, loc="upper center", ncol=2)

    idx = list(P)
    xs = np.arange(len(idx))
    Rp = [rows[i]["R"] for i in idx]
    axes[1].bar(xs, Rp, 0.6, color=[COLORS[i] for i in idx])
    axes[1].axhline(1.0, color="0.4", lw=1.0, ls="--")
    axes[1].set_xticks(xs); axes[1].set_xticklabels([TAG[i] for i in idx])
    axes[1].set_ylabel(r"$R^{(1\mathrm{R})}$ at $\Delta\omega=+2\gamma_2$")
    axes[1].set_title("Enhancement at the operating point")
    for x, i in zip(xs, idx):
        axes[1].text(x, Rp[idx.index(i)] + 0.03, f"{Rp[idx.index(i)]:.2f}", ha="center", fontsize=8.5)
    fig.suptitle("Cascaded $\\chi^{(3)}$: 4 earlier samples + Design 6 "
                 "(S1–4 area-normalised overlap, D6 volume-normalised — see caveat)",
                 y=1.02, fontsize=10.5)
    _save(fig, "combined_R_comparison")


def fig_tradeoff(P, rows):
    """The two routes to R: high-Q x low-overlap (S3,S4) vs low-Q x high-overlap (D6)."""
    # per-sample label offsets so the near-coincident S1/S2 don't collide
    OFF = {1: (-42, 14), 2: (12, -24), 3: (12, 12), 4: (12, 12), 6: (12, 12)}
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    for i, p in P.items():
        R = rows[i]["R"]
        ax.scatter(p["Q2"], rows[i]["ov"], s=90 + 320 * R, color=COLORS[i],
                   alpha=0.8, edgecolor="k", lw=0.6, zorder=3)
        ax.annotate(f"{TAG[i]}, R={R:.2f}", (p["Q2"], rows[i]["ov"]),
                    textcoords="offset points", xytext=OFF[i], fontsize=8.5, color=COLORS[i])
    # iso-R guide: R ∝ Q2 * overlap^2  -> overlap ∝ 1/sqrt(Q2)
    Q = np.linspace(700, 6200, 200)
    for Rref, ls in ((1.0, ":"), (2.0, "--")):
        # calibrate the iso-R curve to S4 (R=2.49 at Q=5633, ov=0.367): const = ov^2 * Q / R
        const = rows[4]["ov"] ** 2 * P[4]["Q2"] / rows[4]["R"]
        ax.plot(Q, np.sqrt(Rref * const / Q), color="0.6", lw=1.0, ls=ls,
                label=f"iso-$R$={Rref:.0f} ($R\\propto Q_2\\,{{\\rm overlap}}^2$)")
    ax.set_xlabel(r"SH quality factor $Q_2$")
    ax.set_ylabel(r"SHG overlap at $\Delta\omega=+2\gamma_2$")
    ax.set_title("Two routes to the same $R$: high-$Q$/low-overlap vs low-$Q$/high-overlap")
    ax.legend(fontsize=8.5, loc="upper right")
    _save(fig, "combined_tradeoff")


def fig_real_device(P):
    """Real-device robustness: (left) pulse r_phi vs tau_p, (right) disorder tolerance."""
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.7))
    tp = np.linspace(0.3, 6.0, 90) * 1e-12
    for i, p in P.items():
        lw = 2.4 if i == 6 else 1.4
        r = np.array([reduction_factors(DOP, p["gamma2"] * w)[0] for w in tp])
        axes[0].plot(tp * 1e12, r, color=COLORS[i], lw=lw,
                     label=fr"{TAG[i]} ($1/\gamma_2$={1e12/p['gamma2']:.2f} ps)")
    axes[0].axvspan(1.0, 3.0, color="#1f77b4", alpha=0.10)
    axes[0].text(2.0, 0.15, "ps laser\n(1–3 ps)", ha="center", fontsize=8.5, color="#1f77b4")
    axes[0].axhline(1.0, color="0.4", lw=1.0, ls="--")
    axes[0].set_ylim(0, 1.15)
    axes[0].set_xlabel(r"pulse duration $\tau_p$ [ps]")
    axes[0].set_ylabel(r"$r_\varphi=\varphi_{\rm pulsed}/\varphi_{\rm quasi\text{-}static}$")
    axes[0].set_title("Pulse-bandwidth robustness (broad pulse $\\Rightarrow$ over-excites high-$Q$)")
    axes[0].legend(fontsize=8, loc="lower right")

    idx = list(P)
    xs = np.arange(len(idx))
    tol = [1.0 / P[i]["Q2"] for i in idx]
    axes[1].bar(xs, tol, 0.6, color=[COLORS[i] for i in idx])
    axes[1].set_xticks(xs); axes[1].set_xticklabels([TAG[i] for i in idx])
    axes[1].set_ylabel(r"max disorder $\sigma_\omega/\omega \lesssim 1/Q_2$")
    axes[1].set_title("Fabrication-disorder tolerance (higher = more forgiving)")
    for x, i in zip(xs, idx):
        axes[1].text(x, tol[idx.index(i)] * 1.02, f"{tol[idx.index(i)]*1e4:.1f}$\\times10^{{-4}}$",
                     ha="center", va="bottom", fontsize=8)
    fig.suptitle(r"Real-device robustness: Design 6's low $Q_2$ wins on BOTH pulse bandwidth and disorder",
                 y=1.02, fontsize=10.5)
    _save(fig, "combined_real_device")


def main():
    P = get_all()
    rows = master_table(P)
    fig_R_comparison(P, rows)
    fig_tradeoff(P, rows)
    fig_real_device(P)


if __name__ == "__main__":
    main()
