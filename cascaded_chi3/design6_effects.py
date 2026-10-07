"""SH-only cascaded-chi3 observable suite for Design 6 (TiO2/LiNbO3 hybrid).

Consumes the verified Design-6 parameters (design6_params.build_design6) and the
data-derived eta_geom(Domega) (volume-normalised COMSOL overlap), then evaluates the
same observables that survive in the singly-resonant case as for samples 1-4
(notes/Cascaded_chi3_SH_only_note.tex, sec:chi3-1R .. sec:applications-1R):

  * chi3_casc(Domega)  Re (self-focus/defocus sign flip), Im (2-photon-like loss), R^(1R)
  * sum rule |chi3_casc| * Df
  * single-pass per-photon Kerr phase phi_NL(I)  + damage-threshold check
  * single-pass induced loss (saturable absorber) dT/T(I)
  * broadband cascaded FWM efficiency eta_FWM(I)
  * SHG conversion efficiency eta_SHG(I)   (bridge + coupling caveated)
  * optimal detuning scan  (data-derived eta over the single-pole window)
  * time-domain pulsed TCMT: phase reduction r_phi(tau_p) + sign-flip survival

All the CW physics functions are P-KEYED and reused verbatim from sample_effects; only
eta_SHG is redefined here to use Design 6's supercell area p["A_cell"] (= 2 p^2) instead
of the samples' (400 nm)^2 module constant.  The pulsed lineshape machinery is reused
from sample_pulsed.  Samples 3 and 4 (the earlier best) are overlaid faintly for context.

    conda run -n ml_base python design6_effects.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from style import apply_style
apply_style()

from design6_params import build_design6, CHI2, CHI3_NATIVE, C_BRIDGE
from sample_params import build_sample, chi3_casc, eta_geom
from analytic import alpha3_eff_casc_full, gamma_2pa_eff
from sample_effects import (E_loc_sq, phi_NL, induced_loss, eta_FWM, EPS0,
                            I_1GW, I_10GW, DAMAGE_BAND, DAMAGE_PS,
                            BULK_PPLN_SUMRULE, F_RAD)
from sample_pulsed import reduction_factors, solve_a2, a2_ss, chi3hat

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE / "figures" / "design6"

D6_COLOR = "#9467bd"                                 # purple
REF_COLORS = {3: "#2ca02c", 4: "#d62728"}            # S3, S4 faint reference
DOP = 2.0


def _label(p):
    return f"Design 6 ({p['lam1_um']*1e3:.0f} nm, $Q_2$={p['Q2']:.0f})"


def _save(fig, name):
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"{name}.{ext}")
    print(f"wrote figures/design6/{name}.{{png,pdf}}")
    plt.close(fig)


# ----------------------------------------------------------------------
# Design-6 SHG efficiency: identical to sample_effects.eta_SHG but with the
# supercell area p["A_cell"] and V_cell = LN volume p["volume"].
# ----------------------------------------------------------------------
def eta_SHG(p, dop, I_in, T_geom=1.0, f_rad=F_RAD):
    eta = eta_geom(p, dop)
    V_cell = p["volume"]                              # = A_cell * L_meta = V_LN
    M2 = eta * V_cell * CHI2 ** 2 / (EPS0 * p["n2"] ** 2)       # |M~|^2 for N2 = 1
    Els = E_loc_sq(p, I_in, T_geom)
    beta_abs = EPS0 * p["omega1"] * Els * np.sqrt(M2)          # beta^ = eps0 w1 E^2 M~ / N2
    domega = dop * p["gamma2"]
    a2_sq = beta_abs ** 2 / (domega ** 2 + (p["gamma2"] / 2.0) ** 2)
    P_SH = f_rad * p["gamma2"] * 2.0 * a2_sq                    # U2 = 2 N2 |a2|^2 (Boyd)
    P_in = I_in * p["A_cell"]
    return P_SH / P_in


# ----------------------------------------------------------------------
def verify(p):
    """Independent correctness cross-checks (mirror sample_optimal.verify for one design)."""
    print("=" * 70)
    print("VERIFICATION CROSS-CHECKS  (Design 6)")
    print("=" * 70)

    # (i) prefactor-free lineshape: Re/Im chi3 == -2*(Domega/gamma2) at fixed eta
    print("\n(i) Re/Im(chi3) == -2*(Domega/gamma2)  [must hold exactly]")
    ok = True
    for dtil in (0.7, 1.0, 2.0, 3.5):
        c = chi3_casc(p, dtil, eta=0.05)
        if abs(c.real / c.imag - (-2 * dtil)) > 1e-9 * abs(2 * dtil):
            ok = False
            print("   MISMATCH dtil=%.2f: Re/Im=%.6f exp %.6f" % (dtil, c.real / c.imag, -2 * dtil))
    print("   -> %s" % ("PASS" if ok else "FAIL"))

    # (ii) match to analytic.py adiabatic dimensionless forms (constant eta)
    print("\n(ii) chi3 lineshape == analytic.py forms (alpha3_eff, gamma_2pa)")
    eta = 0.18
    K = 2.0 * p["omega1"] * CHI2 ** 2 * eta / (3.0 * p["n2"] ** 2)
    worst = 0.0
    for dtil in np.linspace(-5, 5, 21):
        if abs(dtil) < 1e-6:
            continue
        c = chi3_casc(p, dtil, eta=eta)
        re_exp = (K / p["gamma2"]) * alpha3_eff_casc_full(dtil)
        im_exp = (K / p["gamma2"]) * gamma_2pa_eff(dtil)
        worst = max(worst,
                    abs(c.real - re_exp) / abs(re_exp) if re_exp else 0,
                    abs(c.imag - im_exp) / abs(im_exp) if im_exp else 0)
    print("   max relative deviation over dtil in [-5,5]: %.2e  -> %s"
          % (worst, "PASS" if worst < 1e-12 else "FAIL"))

    # (iii) units (by inspection)
    print("\n(iii) units: omega1[1/s]*chi2^2[m^2/V^2]/(n2^2*denom[1/s]) = m^2/V^2  -> PASS")

    # (iv) hand spot value at +2 gamma2
    e = eta_geom(p, DOP)
    denom = abs(DOP * p["gamma2"] + 1j * p["gamma2"] / 2)
    hand = 2 * p["omega1"] * CHI2 ** 2 * e / (3 * p["n2"] ** 2 * denom)
    code = abs(chi3_casc(p, DOP))
    print("\n(iv) |chi3| @ +2gamma2: eta=%.4f hand=%.3e code=%.3e R=%.2f -> %s"
          % (e, hand, code, code / CHI3_NATIVE, "PASS" if abs(hand - code) < 1e-30 else "FAIL"))


# ----------------------------------------------------------------------
def fig_chi3_spectrum(p, refs):
    """Re, Im chi3_casc and R vs Domega/gamma2 (Design 6, with S3/S4 faint reference)."""
    dop = np.linspace(-6, 6, 1201)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3))
    c = np.array([chi3_casc(p, d) for d in dop])
    for i, pr in refs.items():
        cr = np.array([chi3_casc(pr, d) for d in dop])
        axes[0].plot(dop, cr.real * 1e21, color=REF_COLORS[i], lw=1.0, alpha=0.35)
        axes[1].plot(dop, cr.imag * 1e21, color=REF_COLORS[i], lw=1.0, alpha=0.35)
        axes[2].plot(dop, np.abs(cr) / CHI3_NATIVE, color=REF_COLORS[i], lw=1.0, alpha=0.35,
                     label=f"S{i} ref")
    axes[0].plot(dop, c.real * 1e21, color=D6_COLOR, lw=2.2, label=_label(p))
    axes[1].plot(dop, c.imag * 1e21, color=D6_COLOR, lw=2.2)
    axes[2].plot(dop, np.abs(c) / CHI3_NATIVE, color=D6_COLOR, lw=2.2, label="Design 6")
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
    fig.suptitle(r"Design 6 cascaded $\chi^{(3)}$ spectra (dotted = conservative $\Delta\omega=\pm2\gamma_2$)",
                 y=1.02, fontsize=12)
    _save(fig, "design6_chi3_spectra")


def fig_eta(p, refs):
    """eta_geom(Domega) from the volume-normalised overlap; Design 6's high overlap vs samples."""
    dop = np.linspace(-6, 6, 1201)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    axes[0].plot(dop, [eta_geom(p, d) for d in dop], color=D6_COLOR, lw=2.2, label=_label(p))
    for i, pr in refs.items():
        axes[0].plot(dop, [eta_geom(pr, d) for d in dop], color=REF_COLORS[i], lw=1.0,
                     alpha=0.4, label=f"S{i} ref")
    axes[0].axhline(0.1, color="0.4", lw=1.0, ls="--", label="note assumed $\\eta\\sim0.1$")
    for d in (-DOP, DOP):
        axes[0].axvline(d, color="0.6", lw=0.8, ls=":")
    axes[0].set_xlabel(r"$\Delta\omega/\gamma_2$")
    axes[0].set_ylabel(r"$\eta_{\rm geom}^{(1\mathrm{R})}=(|{\rm overlap}|/V_{\rm LN})^2$")
    axes[0].set_title("Data-derived overlap (volume-normalised)")
    axes[0].legend(fontsize=8.5)

    # dimensionless overlap + mode-symmetry raw vs Domega
    f = p["freq_arr"]
    dw = (f - p["f0"]) / p["g_THz"]
    axes[1].plot(dw, p["norm_ov"], "o-", color=D6_COLOR, ms=3, lw=1.4, label=r"$|{\rm overlap}|/V$")
    axes[1].plot(dw, p["norm_sym"], "s-", color="0.5", ms=3, lw=1.0, label=r"$|{\rm symmetry}|/V$")
    for d, lab in ((-DOP, "0.41"), (DOP, "0.44")):
        axes[1].axvline(d, color="0.6", lw=0.8, ls=":")
    axes[1].annotate("pptx: 0.41 / 0.44", xy=(0, 0.5), fontsize=8.5, color="0.3", ha="center")
    axes[1].set_xlim(-6, 6)
    axes[1].set_xlabel(r"$\Delta\omega/\gamma_2$")
    axes[1].set_ylabel("dimensionless overlap")
    axes[1].set_title(r"Overlap + mode-symmetry (Fano dip at $\Delta\omega=0$)")
    axes[1].legend(fontsize=8.5)
    _save(fig, "design6_eta_overlap")


def fig_feasibility(p):
    """phi_NL and induced loss vs incident intensity, with damage bands."""
    I = np.logspace(8, 16.3, 600)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    axes[0].loglog(I * 1e-4, np.abs([phi_NL(p, DOP, ii) for ii in I]), color=D6_COLOR, lw=2.0,
                   label=_label(p))
    axes[1].loglog(I * 1e-4, np.abs([induced_loss(p, DOP, ii) for ii in I]), color=D6_COLOR, lw=2.0)
    for ax in axes:
        ax.axvspan(DAMAGE_BAND[0] * 1e-4, DAMAGE_BAND[1] * 1e-4, color="#d62728", alpha=0.12)
        ax.axvspan(DAMAGE_PS[0] * 1e-4, DAMAGE_PS[1] * 1e-4, color="#1f77b4", alpha=0.12)
        y0 = ax.get_ylim()[0]
        ax.text(np.sqrt(DAMAGE_BAND[0] * DAMAGE_BAND[1]) * 1e-4, y0 * 2.5, "CW/ns\ndamage",
                fontsize=7.5, color="#d62728", ha="center", va="bottom")
        ax.text(np.sqrt(DAMAGE_PS[0] * DAMAGE_PS[1]) * 1e-4, y0 * 2.5, "ps-pulse\nceiling",
                fontsize=7.5, color="#1f77b4", ha="center", va="bottom")
        ax.set_xlabel(r"peak intensity $I_{\rm in}\ [\mathrm{W/cm^2}]$")
    for y, lab, cc in [(0.1, "device level (0.1 rad)", "0.3"),
                       (1e-3, "$1$ mrad (measurable)", "0.45"),
                       (1e-4, "interferometric floor", "0.6")]:
        axes[0].axhline(y, color=cc, lw=0.9, ls="--" if y == 0.1 else ":")
        axes[0].text(2e4, y * 1.3, lab, fontsize=7.5, color=cc)
    axes[0].set_ylabel(r"$|\varphi_{\rm NL}^{(1\mathrm{R})}|$ [rad]")
    axes[0].set_title(r"Single-pass Kerr phase ($\Delta\omega=+2\gamma_2$, $L=150$ nm)")
    axes[1].set_ylabel(r"single-pass $|\Delta T/T|$")
    axes[1].set_title(r"Cascaded saturable absorption ($\Delta\omega=+2\gamma_2$)")
    axes[0].legend(fontsize=8.5, loc="lower right")
    fig.suptitle(r"Design 6: single-pass phase is $\mu$rad-level (SH-only $\Rightarrow$ no CW-Kerr device)",
                 y=1.02, fontsize=10.5)
    _save(fig, "design6_feasibility")


def fig_fwm_shg(p):
    I = np.logspace(8, 16.3, 600)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))
    axes[0].loglog(I * 1e-4, [eta_FWM(p, DOP, ii) for ii in I], color=D6_COLOR, lw=2.0, label=_label(p))
    axes[1].loglog(I * 1e-4, [eta_SHG(p, DOP, ii) for ii in I], color=D6_COLOR, lw=2.0)
    for ax in axes:
        ax.axvspan(DAMAGE_BAND[0] * 1e-4, DAMAGE_BAND[1] * 1e-4, color="#d62728", alpha=0.12)
        ax.axvspan(DAMAGE_PS[0] * 1e-4, DAMAGE_PS[1] * 1e-4, color="#1f77b4", alpha=0.12)
        ax.set_xlabel(r"peak intensity [$\mathrm{W/cm^2}$]")
    axes[0].set_ylabel(r"$\eta_{\rm FWM}=|A_i/A_s|^2$")
    axes[0].set_title(r"Broadband cascaded FWM ($\Delta\omega=+2\gamma_2$)")
    axes[1].set_ylabel(r"$\eta_{\rm SHG}=P_{\rm SH}/P_{\rm in}$")
    axes[1].set_title(r"SHG conversion (baseline $C_{\rm bridge}=1$, $\gamma_{\rm rad}=\gamma_2/2$)")
    axes[0].legend(fontsize=8.5, loc="lower right")
    _save(fig, "design6_fwm_shg")


def fig_optimal(p, refs, window=(1.0, 5.0)):
    """R vs Domega with the optimum marked (data-derived eta over the single-pole window)."""
    dtil = np.r_[np.linspace(-window[1], -window[0], 700), np.linspace(window[0], window[1], 700)]
    absc = np.array([abs(chi3_casc(p, d)) for d in dtil])
    pos = dtil > 0
    neg = dtil < 0
    jp = np.argmax(absc[pos]); dp = dtil[pos][jp]; vp = absc[pos][jp]
    jn = np.argmax(absc[neg]); dn = dtil[neg][jn]; vn = absc[neg][jn]
    ref = abs(chi3_casc(p, DOP))
    print("\n=== OPTIMAL DETUNING (Design 6) ===")
    print("  defocus (Dw>0): max R=%.2f @ Dw/g=%.2f    focus (Dw<0): max R=%.2f @ Dw/g=%.2f"
          % (vp / CHI3_NATIVE, dp, vn / CHI3_NATIVE, dn))
    print("  ref @ +2g: R=%.2f   -> gain at optimum x%.2f  (eta_opt=%.3f @ Dw=%.1f)"
          % (ref / CHI3_NATIVE, vp / ref, eta_geom(p, dp), dp))

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    dd = np.linspace(-6, 6, 1201)
    ax.plot(dd, [abs(chi3_casc(p, d)) / CHI3_NATIVE for d in dd], color=D6_COLOR, lw=2.2, label=_label(p))
    for i, pr in refs.items():
        ax.plot(dd, [abs(chi3_casc(pr, d)) / CHI3_NATIVE for d in dd], color=REF_COLORS[i],
                lw=1.0, alpha=0.4, label=f"S{i} ref")
    ax.plot([dp], [vp / CHI3_NATIVE], "o", color=D6_COLOR, ms=8, mfc="white", mew=1.8)
    ax.plot([dn], [vn / CHI3_NATIVE], "o", color=D6_COLOR, ms=8, mfc="white", mew=1.8)
    ax.axhline(1.0, color="0.4", lw=1.0, ls="--", label="native Kerr ($R=1$)")
    for d in (-DOP, DOP):
        ax.axvline(d, color="0.6", lw=0.8, ls=":")
    ax.axvspan(-1, 1, color="0.85", alpha=0.5)
    ax.text(0, ax.get_ylim()[1] * 0.9, "single-pole\nbreakdown", ha="center", fontsize=8, color="0.4")
    ax.set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax.set_ylabel(r"$R^{(1\mathrm{R})}=|\chi^{(3)}_{\rm casc}|/\chi^{(3)}_{\rm native}$")
    ax.set_title("Design 6 optimal detuning (circles); dotted = conservative $\\pm2\\gamma_2$")
    ax.legend(fontsize=8.5, loc="upper right")
    _save(fig, "design6_optimal_detuning")


def fig_pulsed(p, refs, delta=DOP):
    """Pulsed phase reduction r_phi(tau_p) and sign-flip survival for Design 6.

    Design 6's low Q_2 (=924) gives a SHORT SH lifetime 1/gamma2 ~ 0.43 ps, so the
    quasi-static cascade holds even at tau_p = 1 ps -- unlike the high-Q samples.
    """
    inv_g_ps = 1e12 / p["gamma2"]
    tp = np.linspace(0.3, 6.0, 120) * 1e-12
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))

    # left: r_phi vs tau_p, Design 6 vs S3/S4
    r6 = np.array([reduction_factors(delta, p["gamma2"] * w)[0] for w in tp])
    axes[0].plot(tp * 1e12, r6, color=D6_COLOR, lw=2.2,
                 label=fr"Design 6 ($1/\gamma_2$={inv_g_ps:.2f} ps)")
    for i, pr in refs.items():
        rr = np.array([reduction_factors(delta, pr["gamma2"] * w)[0] for w in tp])
        axes[0].plot(tp * 1e12, rr, color=REF_COLORS[i], lw=1.0, alpha=0.5,
                     label=fr"S{i} ($1/\gamma_2$={1e12/pr['gamma2']:.2f} ps)")
    axes[0].axvspan(1.0, 3.0, color="#1f77b4", alpha=0.10)
    axes[0].text(2.0, 0.12, "ps-laser\n(1--3 ps)", ha="center", fontsize=8.5, color="#1f77b4")
    axes[0].axhline(1.0, color="0.4", lw=1.0, ls="--", label="CW / quasi-static")
    axes[0].set_ylim(0, 1.1)
    axes[0].set_xlabel(r"pulse duration $\tau_p$ (intensity FWHM) [ps]")
    axes[0].set_ylabel(r"$r_\varphi=\varphi_{\rm pulsed}/\varphi_{\rm quasi\text{-}static}$")
    axes[0].set_title(r"Low $Q_2$ $\Rightarrow$ fast SH $\Rightarrow$ no pulse penalty")
    axes[0].legend(fontsize=8.0, loc="lower right")
    r6_123 = np.interp([1, 2, 3], tp * 1e12, r6)

    # right: sign flip survives pulsing (I^2-weighted effective Re chi3 vs detuning)
    dd = np.linspace(-6, 6, 241)
    cw = np.real(chi3hat(dd)); norm = np.max(np.abs(cw))
    axes[1].plot(dd, cw / norm, color="0.4", lw=2, ls="--", label="CW (quasi-static)")
    for tpi, col in zip((1.0, 3.0), ("#d62728", "#2ca02c")):
        W = p["gamma2"] * tpi * 1e-12
        eff = [reduction_factors(d, W)[3] for d in dd]
        axes[1].plot(dd, np.array(eff) / norm, color=col, lw=1.8, label=fr"pulsed $\tau_p={tpi:.0f}$ ps")
    axes[1].axhline(0, color="k", lw=0.5); axes[1].axvline(0, color="k", lw=0.5)
    axes[1].set_xlabel(r"$\Delta\omega/\gamma_2$")
    axes[1].set_ylabel(r"$\mathrm{Re}\,\chi^{(3)}_{\rm eff}$ (norm. to CW peak)")
    axes[1].set_title("Sign flip survives pulsing (Design 6)")
    axes[1].legend(fontsize=8.5, loc="upper right")
    _save(fig, "design6_pulsed")
    return r6_123


# ----------------------------------------------------------------------
def summary(p, r6_123):
    c3 = chi3_casc(p, DOP)
    print("\n=== DESIGN 6 OPERATING-POINT SUMMARY  (Domega=+2 gamma2, C_bridge=1) ===")
    print("  eta_geom          = %.4f" % eta_geom(p, DOP))
    print("  |chi3_casc|       = %.3e m^2/V^2   R^(1R) = %.2f" % (abs(c3), abs(c3) / CHI3_NATIVE))
    print("  Re/Im chi3        = %+.3e / %+.3e   (sign-flip + 2-photon loss)" % (c3.real, c3.imag))
    print("  Df = gamma_THz    = %.1f GHz        |chi3|*Df = %.3e m^2 Hz/V^2"
          % (p["Df_GHz"], abs(c3) * p["Df_GHz"] * 1e9))
    print("  phi_NL @1GW/@10GW = %.2e / %.2e rad (%.0f nm single pass)"
          % (abs(phi_NL(p, DOP, I_1GW)), abs(phi_NL(p, DOP, I_10GW)), p["L_meta"] * 1e9))
    print("  |dT/T| @1GW       = %.2e            eta_FWM @1GW = %.2e" %
          (abs(induced_loss(p, DOP, I_1GW)), eta_FWM(p, DOP, I_1GW)))
    print("  eta_SHG @1GW      = %.2e  (order-of-mag baseline)" % eta_SHG(p, DOP, I_1GW))
    print("  pulsed r_phi @ tau_p=1/2/3 ps = %.3f / %.3f / %.3f  (1/gamma2=%.2f ps -> quasi-static holds)"
          % (r6_123[0], r6_123[1], r6_123[2], 1e12 / p["gamma2"]))


def main():
    p = build_design6()
    refs = {3: build_sample(3), 4: build_sample(4)}
    verify(p)
    fig_chi3_spectrum(p, refs)
    fig_eta(p, refs)
    fig_feasibility(p)
    fig_fwm_shg(p)
    fig_optimal(p, refs)
    r6_123 = fig_pulsed(p, refs)
    summary(p, r6_123)


if __name__ == "__main__":
    main()
