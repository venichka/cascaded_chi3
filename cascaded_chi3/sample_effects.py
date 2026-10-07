"""SH-only cascaded-chi3 effects for the designer's 4 LiNbO3 samples.

Consumes the verified per-sample parameters of ``sample_params.py`` and the
data-derived eta_geom(Domega) (from the COMSOL overlap arrays), then evaluates
every observable that *survives* in the singly-resonant case
(``notes/Cascaded_chi3_SH_only_note.tex`` sections sec:chi3-1R .. sec:applications-1R):

  * chi3_casc(Domega)        Re (self-focus/defocus sign flip), Im (2-photon-like loss)
  * R^(1R) = |chi3_casc|/chi3_native
  * sum rule |chi3_casc| * Df
  * single-pass per-photon Kerr phase phi_NL(I_in)         + damage-threshold check
  * single-pass induced loss (saturable absorber) dT/T(I_in)
  * broadband cascaded FWM efficiency eta_FWM(I_in)
  * SHG conversion efficiency eta_SHG(I_in)                 (bridge+coupling caveated)

Effects that DIE in 1R (bistability, switching, parametric limiter, AM-modulator
readout) are not computed -- there is no FM cavity (note sec:applications-1R).

Run directly to (re)generate all figures into  figures/samples/  and print the
operating-point summary table.

    conda run -n ml_base python sample_effects.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from style import apply_style
apply_style()

import sample_params as sp
from sample_params import CHI2, CHI3_NATIVE, C0, A_CELL, build_sample, chi3_casc, eta_geom

EPS0 = 8.8541878128e-12

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures" / "samples"          # dedicated figure folder

SAMPLE_COLORS = {1: "#1f77b4", 2: "#ff7f0e", 3: "#2ca02c", 4: "#d62728"}
DOP = +2.0                                        # default operating point Domega = +2 gamma2

# benchmark incident intensities [W/m^2]   (1 GW/cm^2 = 1e13 W/m^2)
I_1GW = 1e13
I_10GW = 1e14
# LiNbO3 damage / usable-intensity regimes (W/m^2):
#   CW @ green is photorefractive-limited far below the bulk threshold;
#   ns/CW bulk damage ~ 1-10 GW/cm^2;
#   ps pulses (tau_p >~ 1/gamma2 ~ ps, still adiabatic) reach a fluence-limited
#   ceiling ~ 0.1-1 J/cm^2 / few-ps ~ 0.1-1 TW/cm^2 before damage.
DAMAGE_BAND = (1e13, 1e14)          # ns / CW bulk  (1-10 GW/cm^2)
DAMAGE_PS = (1e15, 1e16)            # ps-pulse ceiling (~0.1-1 TW/cm^2)
# bulk PPLN sum-rule reference (note sec:numerical, after 1/n2 correction)
BULK_PPLN_SUMRULE = 2e-8                          # m^2 Hz / V^2

# SHG baseline assumption: radiative fraction of gamma2 (order-of-magnitude)
F_RAD = 0.5


# ----------------------------------------------------------------------
# physics (single-pass, SH-only)
# ----------------------------------------------------------------------
def E_loc_sq(p, I_in, T_geom=1.0):
    """|E_loc|^2 [V^2/m^2] inside the dielectric for incident intensity I_in [W/m^2].

    Boyd amplitude convention (E~ = E e^{-iwt} + c.c., I = 2 n eps0 c |E|^2):
    |E_0|^2 = I_in / (2 eps0 c) in vacuum, |E_loc| = T_geom |E_0| / n1, so
        |E_loc|^2 = T_geom^2 I_in / (2 eps0 c n1^2).
    (Until 2026-10-06 this used the 1/2-convention value 2 T^2 I/(eps0 c n1^2),
    4x larger; that error cancelled against the 4x-too-small chi3_casc in
    phi_NL, induced_loss and eta_FWM, which are therefore unchanged.)
    """
    return T_geom ** 2 * I_in / (2.0 * EPS0 * C0 * p["n1"] ** 2)


def phi_NL(p, dop, I_in, T_geom=1.0):
    """Single-pass per-photon Kerr phase [rad] (note eq:phi-1R)."""
    chi3 = chi3_casc(p, dop)
    return (3.0 * p["omega1"] / (2.0 * C0 * p["n1"])) * chi3.real * E_loc_sq(p, I_in, T_geom) * p["L_meta"]


def induced_loss(p, dop, I_in, T_geom=1.0):
    """Single-pass intensity-dependent loss |dT/T| (note sec:applications-1R item 6).

    alpha_NL L = (3 omega1 / (n1 c)) Im chi3 |E_loc|^2 L  (twice the phase coefficient).
    """
    chi3 = chi3_casc(p, dop)
    return (3.0 * p["omega1"] / (C0 * p["n1"])) * chi3.imag * E_loc_sq(p, I_in, T_geom) * p["L_meta"]


def eta_FWM(p, dop, I_in, T_geom=1.0):
    """Degenerate-pump cascaded FWM efficiency |A_i/A_s|^2 (note eq:FWM-1R-SVEA)."""
    chi3 = chi3_casc(p, dop)
    coeff = 3.0 * p["omega1"] / (2.0 * C0 * p["n1"])
    return (coeff * abs(chi3) * E_loc_sq(p, I_in, T_geom) * p["L_meta"]) ** 2


def eta_SHG(p, dop, I_in, T_geom=1.0, f_rad=F_RAD):
    """SHG conversion efficiency P_SH/P_in (order-of-magnitude; bridge+coupling caveated).

    |M~|^2 = eta_geom V_cell |chi2|^2 / (eps0 n2^2)  (N2 = int eps0 eps |f2|^2 = 1);
    beta^ = eps0 omega1 |E_loc|^2 M~ / N2;   |a2|^2 = |beta^|^2 / (Domega^2 + (gamma2/2)^2);
    stored energy U2 = 2 N2 |a2|^2 (Boyd);   P_SH = f_rad gamma2 U2.
    (Corrected 2026-10-06: the old chain beta^ = (eps0 w1/4)|E|^2 M~ with U2 = |a2|^2
    underestimated P_SH by 2.)
    """
    eta = eta_geom(p, dop)
    V_cell = A_CELL * p["L_meta"]
    M2 = eta * V_cell * CHI2 ** 2 / (EPS0 * p["n2"] ** 2)        # |M~|^2
    Els = E_loc_sq(p, I_in, T_geom)
    beta_abs = EPS0 * p["omega1"] * Els * np.sqrt(M2)            # |beta^|, N2 = 1
    domega = dop * p["gamma2"]
    a2_sq = beta_abs ** 2 / (domega ** 2 + (p["gamma2"] / 2.0) ** 2)
    P_SH = f_rad * p["gamma2"] * 2.0 * a2_sq                     # U2 = 2 N2 |a2|^2
    P_in = I_in * A_CELL
    return P_SH / P_in


# ----------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------
def _label(p):
    return f"S{p['idx']} ({p['lam1_um']*1e3:.0f} nm, $Q_2$={p['Q2']:.0f})"


def fig_chi3_spectra(samples):
    """Re, Im chi3_casc and R vs Domega/gamma2 for all 4 samples."""
    dop = np.linspace(-6, 6, 1201)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3))
    for i, p in samples.items():
        c = np.array([chi3_casc(p, d) for d in dop])
        col = SAMPLE_COLORS[i]
        axes[0].plot(dop, c.real * 1e21, color=col, label=_label(p))
        axes[1].plot(dop, c.imag * 1e21, color=col)
        axes[2].plot(dop, np.abs(c) / CHI3_NATIVE, color=col)
    for ax in axes:
        ax.axvline(0, color="k", lw=0.5)
        for d in (-DOP, DOP):
            ax.axvline(d, color="0.6", lw=0.8, ls=":")
        ax.set_xlabel(r"$\Delta\omega/\gamma_2$")
    axes[0].axhline(0, color="k", lw=0.5)
    axes[0].set_ylabel(r"$\mathrm{Re}\,\chi^{(3)}_{\rm casc}\ [10^{-21}\,\mathrm{m^2/V^2}]$")
    axes[0].set_title("Self-focusing $\\leftrightarrow$ defocusing (sign flip)")
    axes[1].set_ylabel(r"$\mathrm{Im}\,\chi^{(3)}_{\rm casc}\ [10^{-21}\,\mathrm{m^2/V^2}]$")
    axes[1].set_title("Two-photon-like loss")
    axes[2].axhline(1.0, color="0.4", lw=1.0, ls="--", label="native Kerr ($R=1$)")
    axes[2].set_ylabel(r"$R^{(1\mathrm{R})}=|\chi^{(3)}_{\rm casc}|/\chi^{(3)}_{\rm native}$")
    axes[2].set_title("Enhancement over native Kerr")
    axes[0].legend(fontsize=8.5, loc="upper right")
    axes[2].legend(fontsize=8.5, loc="upper right")
    fig.suptitle(r"Cascaded $\chi^{(3)}$ spectra (dotted = conservative $\Delta\omega=\pm2\gamma_2$ operating points)",
                 y=1.02, fontsize=12)
    _save(fig, "sample_chi3_spectra")


def fig_eta_and_ranking(samples):
    """eta_geom(Domega) from the overlap data, and R / sum-rule bar charts."""
    dop = np.linspace(-6, 6, 1201)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3))
    for i, p in samples.items():
        eta = np.array([eta_geom(p, d) for d in dop])
        axes[0].plot(dop, eta, color=SAMPLE_COLORS[i], label=_label(p))
    axes[0].axhline(0.1, color="0.4", lw=1.0, ls="--", label="note assumed $\\eta\\sim0.1$")
    for d in (-DOP, DOP):
        axes[0].axvline(d, color="0.6", lw=0.8, ls=":")
    axes[0].set_xlabel(r"$\Delta\omega/\gamma_2$")
    axes[0].set_ylabel(r"$\eta_{\rm geom}^{(1\mathrm{R})}=(|{\rm overlap}|/{\rm area})^2$")
    axes[0].set_title("Data-derived geometric overlap")
    axes[0].legend(fontsize=8.5)

    idx = list(samples)
    xs = np.arange(len(idx))
    Rp = [abs(chi3_casc(samples[i], DOP)) / CHI3_NATIVE for i in idx]
    Rm = [abs(chi3_casc(samples[i], -DOP)) / CHI3_NATIVE for i in idx]
    cols = [SAMPLE_COLORS[i] for i in idx]
    axes[1].bar(xs - 0.2, Rp, 0.4, color=cols, label=r"$\Delta\omega=+2\gamma_2$")
    axes[1].bar(xs + 0.2, Rm, 0.4, color=cols, alpha=0.5, label=r"$\Delta\omega=-2\gamma_2$")
    axes[1].axhline(1.0, color="0.4", lw=1.0, ls="--")
    axes[1].set_xticks(xs); axes[1].set_xticklabels([f"S{i}" for i in idx])
    axes[1].set_ylabel(r"$R^{(1\mathrm{R})}$")
    axes[1].set_title("Enhancement at operating point")
    axes[1].legend(fontsize=8.5)
    for x, r in zip(xs, Rp):
        axes[1].text(x - 0.2, r + 0.04, f"{r:.2f}", ha="center", va="bottom", fontsize=8)

    sr = [abs(chi3_casc(samples[i], DOP)) * samples[i]["g_THz"] * 1e12 for i in idx]
    axes[2].bar(xs, sr, 0.6, color=cols)
    axes[2].axhline(CHI3_NATIVE * 1e11, color="0.4", lw=1.0, ls="--",
                    label=r"native $\times$ 100 GHz")
    axes[2].axhline(BULK_PPLN_SUMRULE, color="k", lw=1.0, ls="-.",
                    label="bulk PPLN (note)")
    axes[2].set_yscale("log")
    axes[2].set_xticks(xs); axes[2].set_xticklabels([f"S{i}" for i in idx])
    axes[2].set_ylabel(r"$|\chi^{(3)}_{\rm casc}|\cdot\Delta f\ [\mathrm{m^2\,Hz/V^2}]$")
    axes[2].set_title("Bandwidth-integrated strength")
    axes[2].legend(fontsize=8.5)
    _save(fig, "sample_eta_and_ranking")


def fig_feasibility(samples):
    """phi_NL and induced loss vs incident intensity: CW/ns vs ps-pulse regimes."""
    I = np.logspace(8, 16.3, 600)                    # 1e4 .. ~2e12 W/cm^2
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    for i, p in samples.items():
        col = SAMPLE_COLORS[i]
        axes[0].loglog(I * 1e-4, np.abs([phi_NL(p, DOP, ii) for ii in I]), color=col, label=_label(p))
        axes[1].loglog(I * 1e-4, np.abs([induced_loss(p, DOP, ii) for ii in I]), color=col)
    for ax in axes:
        # two usable-intensity ceilings
        ax.axvspan(DAMAGE_BAND[0] * 1e-4, DAMAGE_BAND[1] * 1e-4, color="#d62728", alpha=0.12)
        ax.axvspan(DAMAGE_PS[0] * 1e-4, DAMAGE_PS[1] * 1e-4, color="#1f77b4", alpha=0.12)
        y0 = ax.get_ylim()[0]
        ax.text(np.sqrt(DAMAGE_BAND[0] * DAMAGE_BAND[1]) * 1e-4, y0 * 2.5, "CW/ns\ndamage",
                fontsize=7.5, color="#d62728", ha="center", va="bottom")
        ax.text(np.sqrt(DAMAGE_PS[0] * DAMAGE_PS[1]) * 1e-4, y0 * 2.5, "ps-pulse\nceiling",
                fontsize=7.5, color="#1f77b4", ha="center", va="bottom")
        ax.set_xlabel(r"peak intensity $I_{\rm in}\ [\mathrm{W/cm^2}]$")
    for y, lab, c in [(0.1, "device level (0.1 rad)", "0.3"),
                      (1e-3, "$1$ mrad (measurable)", "0.45"),
                      (1e-4, "interferometric floor", "0.6")]:
        axes[0].axhline(y, color=c, lw=0.9, ls="--" if y == 0.1 else ":")
        axes[0].text(2e4, y * 1.3, lab, fontsize=7.5, color=c)
    axes[0].set_ylabel(r"$|\varphi_{\rm NL}^{(1\mathrm{R})}|$ [rad]")
    axes[0].set_title(r"Single-pass Kerr phase ($\Delta\omega=+2\gamma_2$)")
    axes[1].set_ylabel(r"single-pass $|\Delta T/T|$")
    axes[1].set_title(r"Cascaded saturable absorption ($\Delta\omega=+2\gamma_2$)")
    axes[0].legend(fontsize=8.0, loc="lower right")
    fig.suptitle(r"Pulsed (ps) access $\sim\!10^2$--$10^3\times$ higher peak intensity than CW/ns "
                 r"(adiabatic if $\tau_p\gtrsim1/\gamma_2\sim$ ps)", y=1.02, fontsize=10.5)
    _save(fig, "sample_feasibility")


def fig_fwm_shg(samples):
    """FWM and SHG conversion efficiency vs intensity."""
    I = np.logspace(8, 16.3, 600)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))
    for i, p in samples.items():
        col = SAMPLE_COLORS[i]
        axes[0].loglog(I * 1e-4, [eta_FWM(p, DOP, ii) for ii in I], color=col, label=_label(p))
        axes[1].loglog(I * 1e-4, [eta_SHG(p, DOP, ii) for ii in I], color=col)
    for ax in axes:
        ax.axvspan(DAMAGE_BAND[0] * 1e-4, DAMAGE_BAND[1] * 1e-4, color="#d62728", alpha=0.12)
        ax.axvspan(DAMAGE_PS[0] * 1e-4, DAMAGE_PS[1] * 1e-4, color="#1f77b4", alpha=0.12)
        ax.set_xlabel(r"peak intensity [$\mathrm{W/cm^2}$]")
    axes[0].set_ylabel(r"$\eta_{\rm FWM}=|A_i/A_s|^2$")
    axes[0].set_title(r"Broadband cascaded FWM ($\Delta\omega=+2\gamma_2$)")
    axes[1].set_ylabel(r"$\eta_{\rm SHG}=P_{\rm SH}/P_{\rm in}$")
    axes[1].set_title(r"SHG conversion (baseline: $C_{\rm bridge}=1$, $\gamma_{\rm rad}=\gamma_2/2$)")
    axes[0].legend(fontsize=8.5, loc="lower right")
    _save(fig, "sample_fwm_shg")


def _save(fig, name):
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"{name}.{ext}")
    print(f"wrote figures/samples/{name}.{{png,pdf}}")
    plt.close(fig)


# ----------------------------------------------------------------------
def summary(samples):
    print("\n=== OPERATING-POINT SUMMARY  (Domega = +2 gamma2, C_bridge=1) ===")
    hdr = ("S", "eta", "|chi3|", "R", "Df[GHz]", "|c3|Df", "phi@1GW", "phi@10GW", "loss@1GW", "FWM@1GW", "SHG@1GW")
    print("%-2s %7s %9s %6s %8s %9s %9s %9s %9s %9s %9s" % hdr)
    for i, p in samples.items():
        c3 = chi3_casc(p, DOP)
        print("%-2d %7.4f %9.2e %6.2f %8.1f %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e" % (
            i, eta_geom(p, DOP), abs(c3), abs(c3) / CHI3_NATIVE, p["g_THz"] * 1e3,
            abs(c3) * p["g_THz"] * 1e12,
            abs(phi_NL(p, DOP, I_1GW)), abs(phi_NL(p, DOP, I_10GW)),
            abs(induced_loss(p, DOP, I_1GW)),
            eta_FWM(p, DOP, I_1GW), eta_SHG(p, DOP, I_1GW)))
    print("(phi/loss/FWM/SHG at 1 or 10 GW/cm^2 incident; SHG = order-of-magnitude baseline.)")


def main():
    samples = {i: build_sample(i) for i in (1, 2, 3, 4)}
    fig_chi3_spectra(samples)
    fig_eta_and_ranking(samples)
    fig_feasibility(samples)
    fig_fwm_shg(samples)
    summary(samples)


if __name__ == "__main__":
    main()
